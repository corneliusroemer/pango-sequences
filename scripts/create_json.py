"""
Creates summary.json including the following for each sequence:
    - Unaliased
    - Parent
    - Nuc substitutions
    - AA substitutions
    - Nuc deletions
    - AA deletions
    - New nuc mutations vs parent
    - New AA mutations vs parent
    - New nuc deletions vs parent
    - New AA deletions vs parent
    - Reversion nuc mutations vs parent
    - Reversion AA mutations vs parent
    - Undeletion nuc mutations vs parent
    - Undeletion AA mutations vs parent
    # - Designation date
    # - Sequence last updated

Usage:
python scripts/create_json.py \
    --nextclade-tsv {input.nextclade.tsv} \
    --previous-summary {input.previous_summary} \
    --output {output}
"""

from copy import deepcopy
import datetime
import json
import typer
import pandas as pd
from pango_aliasor.aliasor import Aliasor
import requests
from natsort import natsort_keygen
from deepdiff import DeepDiff


def del_all(mapping, to_remove):
    """Remove list of elements from mapping."""
    for key in to_remove:
        mapping.pop(key, None)


def split_comma(input: str) -> list[str]:
    # Split comma-separated string into list
    # However, instead of return list with empty string, return empty list
    retval = input.split(",")
    return retval if retval != [""] else []


def create_summary(
    nextclade_tsv: str = typer.Option(
        "build/nextclade.tsv", "--nextclade-tsv", help="Nextclade TSV file"
    ),
    # previous_summary: str = typer.Option(
    #     "data/pango-consensus-sequences_summary.json",
    #     "--previous-summary",
    #     help="Previous summary file",
    # ),
    # designation_dates:  str = typer.Option(..., "--designation-dates", help="Designation dates file"),
    output: str = typer.Option(
        "data/pango-consensus-sequences_summary.json",
        "--output",
        help="Output file",
    ),
):
    ALIAS_FILE = "build/alias_key.json"
    ALIAS_KEY_JSON_URL = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json"
    DESIGNATION_DATES_URL = "https://raw.githubusercontent.com/corneliusroemer/pango-designation-dates/main/data/lineage_designation_date.csv"
    with open(ALIAS_FILE, "w") as f:
        f.write(requests.get(ALIAS_KEY_JSON_URL).text)
    aliasor = Aliasor(ALIAS_FILE)

    df = pd.read_csv(
        nextclade_tsv,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        index_col="seqName",
    )

    df["unaliased"] = df.apply(lambda x: aliasor.uncompress(x.name), axis=1)

    df.sort_values(by="unaliased", key=natsort_keygen(), inplace=True)

    r = {}

    # Add data from Nextclade
    for seq_name, row in df.iterrows():
        r[seq_name] = {
            "lineage": seq_name,
            "unaliased": row["unaliased"],
            "parent": None,
            "children": [],
            "nextstrainClade": row["clade_nextstrain"],
            "nucSubstitutions": split_comma(row["substitutions"]),
            "aaSubstitutions": split_comma(row["aaSubstitutions"]),
            "nucDeletions": split_comma(row["deletions"]),
            "aaDeletions": split_comma(row["aaDeletions"]),
            "nucSubstitutionsNew": None,
            "aaSubstitutionsNew": None,
            "nucDeletionsNew": None,
            "aaDeletionsNew": None,
            "nucSubstitutionsReverted": None,
            "aaSubstitutionsReverted": None,
            "nucDeletionsReverted": None,
            "aaDeletionsReverted": None,
            "frameShifts": split_comma(row["frameShifts"]),
            "designationDate": None,
            # "sequenceLastUpdated": None,
            # "entryLastUpdated": None,
        }

    aliasor.parent = lambda x: aliasor.compress(
        ".".join(aliasor.uncompress(x).split(".")[:-1])
    )

    for key, val in r.items():
        parent = None
        candidate = key
        while parent not in r and parent != "":
            parent = aliasor.parent(candidate)
            candidate = parent

        r[key]["parent"] = parent
        if parent in r:
            r[parent]["children"].append(key)
            r[parent]["children"].sort(key=natsort_keygen())

    for key, val in r.items():
        for field in [
            "nucSubstitutions",
            "aaSubstitutions",
            "nucDeletions",
            "aaDeletions",
        ]:
            parent = val["parent"]
            child_val = dict.fromkeys(val[field])
            if parent == "":
                parent_val = {}
            else:
                parent_val = dict.fromkeys(r[parent][field])

            new = deepcopy(child_val)

            del_all(new, parent_val.keys())
            del_all(parent_val, child_val.keys())

            new = list(new.keys())
            reversions = list(parent_val.keys())

            r[key][field + "New"] = new
            r[key][field + "Reverted"] = reversions

    def mut_to_dict(mut: str) -> (str, str, str):
        gene = ""
        if ":" in mut:
            # Amino acid
            gene, mut = mut.split(":")
            gene = gene + ":"
        ref = mut[0]
        alt = mut[-1]
        pos = gene + mut[1:-1]
        return pos, ref, alt
    
    def format_mut(pos: str, ref: str, alt: str) -> str:
        if ":" in pos:
            gene, pos = pos.split(":")
            gene = gene + ":"
        else:
            gene = ""
        return gene + ref + pos + alt

    # Get mutations relative to parent
    for key, val in r.items():
        for field in [
            "nucSubstitutions",
            "aaSubstitutions",
        ]:
            field_relative = field + "Relative"
            parent = val["parent"]
            parent_muts = r.get(parent,{}).get(field, [])
            child_muts = val[field]
            parent_dict = {}
            child_dict = {}
            for mut in parent_muts:
                pos, ref, alt = mut_to_dict(mut)
                parent_dict[pos] = (ref,alt)
            for mut in child_muts:
                pos, ref, alt = mut_to_dict(mut)
                child_dict[pos] = (ref,alt)
            r[key][field_relative] = [] 
            for pos, mut in child_dict.items():
                ref, alt = mut
                parent_ref, parent_alt = parent_dict.get(pos, (None, None))
                if parent_alt is not None:
                    if parent_alt != alt:
                        # Parent and child have different mutations at this position
                        r[key][field_relative].append(format_mut(pos, parent_alt, alt))
                else:
                    r[key][field_relative].append(format_mut(pos, ref, alt))
            for pos, mut in parent_dict.items():
                ref, alt = mut
                alt = child_dict.get(pos, (None, None))
                if pos not in child_dict:
                    # Reversion to reference
                    r[key][field_relative].append(format_mut(pos, alt, ref))
        

    # Add data from designation dates
    designation_dates = pd.read_csv(
        DESIGNATION_DATES_URL,
        dtype=str,
        keep_default_na=False,
        index_col="lineage",
    )

    for key, val in r.items():
        if key in designation_dates.index:
            r[key]["designationDate"] = designation_dates.loc[
                key, "designation_date"
            ]

    # # Add date of last change
    # if previous_summary not in ["", None]:
    #     with open(previous_summary, "r") as f:
    #         q = json.load(f)
    #     for key, val in r.items():
    #         r[key]["entryLastUpdated"] = q[key]["entryLastUpdated"]
    #         r[key]["sequenceLastUpdated"] = q[key]["sequenceLastUpdated"]
    #         if (
    #             val["nucSubstitutions"] != q[key]["nucSubstitutions"]
    #             or val["nucDeletions"] != q[key]["nucDeletions"]
    #         ):
    #             r[key]["sequenceLastUpdated"] = datetime.datetime.now().strftime(
    #                 "%Y-%m-%dT%H:%M:%SZ"
    #             )
    #         if DeepDiff(val, q[key], ignore_order=True):
    #             r[key]["entryLastUpdated"] = datetime.datetime.now().strftime(
    #                 "%Y-%m-%dT%H:%M:%SZ"
    #             )

    # Write out
    with open(output, "w") as f:
        json.dump(r, f, indent=1)


if __name__ == "__main__":
    typer.run(create_summary)
