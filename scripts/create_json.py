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
            "unaliased": row["unaliased"],
            "parent": None,
            "nextstrainClade": row["clade_nextstrain"],
            "nucSubstitutions": row["substitutions"].split(","),
            "aaSubstitutions": row["aaSubstitutions"].split(","),
            "nucDeletions": row["deletions"].split(","),
            "aaDeletions": row["aaDeletions"].split(","),
            "nucSubstitutionsNew": None,
            "aaSubstitutionsNew": None,
            "nucDeletionsNew": None,
            "aaDeletionsNew": None,
            "nucSubstitutionsReverted": None,
            "aaSubstitutionsReverted": None,
            "nucDeletionsReverted": None,
            "aaDeletionsReverted": None,
            "frameShifts": row["frameShifts"].split(","),
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
