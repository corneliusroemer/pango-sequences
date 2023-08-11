# %%
from itertools import chain
import json


xbb = {
    "nucSubstitutions": [
        "A1G",
        "C44T",
        "C241T",
        "A405G",
        "T670G",
        "C2790T",
        "C3037T",
        "G4184A",
        "C4321T",
        "C9344T",
        "A9424G",
        "C9534T",
        "C9866T",
        "C10029T",
        "C10198T",
        "G10447A",
        "C10449A",
        "C12880T",
        "C14408T",
        "G15451A",
        "C15714T",
        "C15738T",
        "T15939C",
        "T16342C",
        "C17410T",
        "T17859C",
        "A18163G",
        "A19326G",
        "C19955T",
        "A20055G",
        "C21618T",
        "T21810C",
        "G21987A",
        "C22000A",
        "C22109G",
        "T22200A",
        "G22577C",
        "G22578A",
        "G22599C",
        "C22664A",
        "C22674T",
        "T22679C",
        "C22686T",
        "A22688G",
        "G22775A",
        "A22786C",
        "G22813T",
        "T22882G",
        "G22895C",
        "T22896C",
        "G22898A",
        "T22942G",
        "G22992A",
        "C22995A",
        "A23013C",
        "T23019C",
        "T23031C",
        "A23055G",
        "A23063T",
        "T23075C",
        "A23403G",
        "C23525T",
        "T23599G",
        "C23604A",
        "C23854A",
        "G23948T",
        "A24424T",
        "T24469A",
        "C25000T",
        "C25416T",
        "C25584T",
        "C26060T",
        "C26270T",
        "A26275G",
        "C26577G",
        "G26709A",
        "C26858T",
        "A27259C",
        "G27382C",
        "A27383T",
        "T27384C",
        "C27807T",
        "A28271T",
        "C28311T",
        "G28881A",
        "G28882A",
        "G28883C",
        "A29510C",
    ],
    "aaSubstitutions": [
        "E:T9I",
        "E:T11A",
        "M:Q19E",
        "M:A63T",
        "N:P13L",
        "N:R203K",
        "N:G204R",
        "N:S413R",
        "ORF1a:K47R",
        "ORF1a:S135R",
        "ORF1a:T842I",
        "ORF1a:G1307S",
        "ORF1a:L3027F",
        "ORF1a:T3090I",
        "ORF1a:L3201F",
        "ORF1a:T3255I",
        "ORF1a:P3395H",
        "ORF1b:P314L",
        "ORF1b:G662S",
        "ORF1b:S959P",
        "ORF1b:R1315C",
        "ORF1b:I1566V",
        "ORF1b:T2163I",
        "ORF3a:T223I",
        "ORF6:D61L",
        "ORF9b:P10S",
        "S:T19I",
        "S:A27S",
        "S:V83A",
        "S:G142D",
        "S:H146Q",
        "S:Q183E",
        "S:V213E",
        "S:G339H",
        "S:R346T",
        "S:L368I",
        "S:S371F",
        "S:S373P",
        "S:S375F",
        "S:T376A",
        "S:D405N",
        "S:R408S",
        "S:K417N",
        "S:N440K",
        "S:V445P",
        "S:G446S",
        "S:N460K",
        "S:S477N",
        "S:T478K",
        "S:E484A",
        "S:F486S",
        "S:F490S",
        "S:Q498R",
        "S:N501Y",
        "S:Y505H",
        "S:D614G",
        "S:H655Y",
        "S:N679K",
        "S:P681H",
        "S:N764K",
        "S:D796Y",
        "S:Q954H",
        "S:N969K",
    ],
    "nucDeletions": [
        "11288-11296",
        "21633-21641",
        "21992-21994",
        "28362-28370",
        "29734-29758",
    ],
    "aaDeletions": [
        "N:E31-",
        "N:R32-",
        "N:S33-",
        "ORF1a:S3675-",
        "ORF1a:G3676-",
        "ORF1a:F3677-",
        "ORF9b:E27-",
        "ORF9b:N28-",
        "ORF9b:A29-",
        "S:L24-",
        "S:P25-",
        "S:P26-",
        "S:Y144-",
    ],
}

GENEMAP = {
    "ORF1a": {
        "start": 266,
        "end": 13468,
    },
    "ORF1b": {
        "start": 13468,
        "end": 21555,
    },
    "S": {
        "start": 21563,
        "end": 25384,
    },
    "ORF3a": {
        "start": 25393,
        "end": 26220,
    },
    "E": {
        "start": 26245,
        "end": 26472,
    },
    "M": {
        "start": 26523,
        "end": 27191,
    },
    "ORF6": {
        "start": 27202,
        "end": 27387,
    },
    "ORF7a": {
        "start": 27394,
        "end": 27759,
    },
    "ORF7b": {
        "start": 27756,
        "end": 27887,
    },
    "ORF8": {
        "start": 27894,
        "end": 28259,
    },
    "N": {
        "start": 28274,
        "end": 29533,
    },
    "ORF9b": {
        "start": 28284,
        "end": 28577,
    },
}


# Function to take amino acid mutation and return corresponding codon
# Takes string "ORF9b:T5I" and returns [28296, 28297, 28298]
def nucs_from_aa(aa_mut):
    gene, mut = aa_mut.split(":")
    aa_pos = int(mut[1:-1])
    nuc_pos = GENEMAP[gene]["start"] + 3 * (aa_pos - 1)
    nucs = [nuc_pos, nuc_pos + 1, nuc_pos + 2]
    return nucs


def nuc_pos_from_nuc(nuc: str) -> int:
    # Takes nuc mut and returns corresponding nuc pos
    # nuc_pos_from_nuc("T28144C") -> 28144
    return int(nuc[1:-1])


def responsible_nucs(aa: str, nuc_subs: list[str]) -> list[str]:
    # Takes aa mut and returns list of nuc muts in corresponding codon
    # responsible_nucs("ORF8:L84S") -> ["T28144C"]
    nucs = nucs_from_aa(aa)
    retval = [nuc for nuc in nuc_subs if nuc_pos_from_nuc(nuc) in nucs]
    return retval


# # %%
# Some aa substitutions are caused by nuc deletions
# Need to also check if a deletion may be responsible for the aa substitution
# Let's treat deletions as ranges, as they are usually single events
# An aa may be affected if the deletion range overlaps with the aa range


# def nuc_del_overlap(nuc_pos: int, del_range: str) -> bool:
#     """
#     Check if a nuc position overlaps with a deletion range
#     """
#     if not del_range:
#         print("Empty deletion range")
#         raise ValueError
#     if "-" not in del_range:
#         del_start = del_end = int(del_range)
#     else:
#         del_start, del_end = map(int, del_range.split("-"))
#     return del_start <= nuc_pos <= del_end


# def responsible_dels(aa: str, dels: List[str]) -> List[str]:
#     """
#     Given an amino acid substitution, find the deletions that may be responsible for it
#     """
#     # Find nuc positions that may be responsible for the aa substitution
#     nucs = nucs_from_aa(aa)

#     retval = [
#         del_range
#         for del_range in dels
#         if any(nuc_del_overlap(nuc, del_range) for nuc in nucs)
#     ]

#     return retval


# %%
# Attribute AA deletions to nuc deletions


# def nuc_del_causing_aa_del(aa_del: str, nuc_dels: List[str]) -> List[str]:
#     """
#     Given an amino acid deletion, find the deletions that may be responsible for it
#     """
#     # Find nuc positions that may be responsible for the aa substitution
#     nucs = nucs_from_aa(aa_del)

#     retval = [
#         del_range
#         for del_range in nuc_dels
#         if any(nuc_del_overlap(nuc, del_range) for nuc in nucs)
#     ]

#     return retval


# %%
# Make dict of nuc substitutions, where key is nuc pos int and value is dict with ref/alt keys and values


def nuc_sub_dict(nuc_subs: list[str]) -> dict[int, dict[str, str]]:
    """
    Make dict of nuc substitutions, where key is nuc pos int and value is dict with ref/alt keys and values
    """
    retval = {}
    for nuc_sub in nuc_subs:
        nuc_pos = int(nuc_pos_from_nuc(nuc_sub))
        ref, alt = nuc_sub[0], nuc_sub[-1]
        retval[nuc_pos] = {"ref": ref, "alt": alt}
    return retval


# %%


# Make list of AA changes, where each item is a dict with keys "gene", "pos", "ref", "alt", and associated "nucPos"
# Treat substitutions and deletions on equal footing
# Hence needs input: aaSubstitutions, aaDeletions, nucSubstitutions, nucDeletions
def aa_changes(
    aa_subs: list[str],
    aa_dels: list[str],
    nuc_dict: dict[int, dict[str, str]],
) -> dict[str, dict[int, dict]]:
    retval = {}
    for aa_sub in chain(aa_subs, aa_dels):
        gene, mut = aa_sub.split(":")
        aa_pos = int(mut[1:-1])
        codon_start = GENEMAP[gene]["start"] + 3 * (aa_pos - 1)
        nucs = [codon_start, codon_start + 1, codon_start + 2]
        # Check which of the nuc subs and nuc del positions overlap the codon
        nuc_pos = set()
        for pos in nucs:
            if pos in nuc_dict:
                nuc_pos.add(pos)
        ref, alt = mut[0], mut[-1]
        if gene not in retval:
            retval[gene] = {}
        retval[gene][aa_pos] = {
            "ref": ref,
            "alt": alt,
            "nucPos": list(nuc_pos),
        }
    return retval


# %%
def massage_frameshifts(
    frameshifts: list[str],
    nuc_dict: dict[int, dict[str, str]],
) -> dict[str, dict[int, dict]]:
    retval = {}
    for fs in frameshifts:
        gene, mut = fs.split(":")
        if "-" in mut:
            start, end = map(int, mut.split("-"))
        else:
            start = end = int(mut)
        for side, aa_pos in ("start", start), ("end", end):
            codon_start = GENEMAP[gene]["start"] + 3 * (aa_pos - 1)
            nucs = [codon_start, codon_start + 1, codon_start + 2]
            # Check which of the nuc subs and nuc del positions overlap the codon
            nuc_pos = set()
            for pos in nucs:
                if pos in nuc_dict:
                    nuc_pos.add(pos)
            if gene not in retval:
                retval[gene] = {}
            retval[gene][aa_pos] = {
                "ref": "",
                "alt": "frameShiftStart"
                if side == "start"
                else "frameShiftEnd",
                "nucPos": list(nuc_pos),
            }
    return retval


# %%
# Put deletion ranges itemized or not based on whether the range appears in the aa_changes
# If it does, then it is itemized
# If it does not, then it is not itemized


# This is code to decide whether a deletion range can be aggregated/itemized or not
# Currently we don't use deletion ranges so we don't need this
# def delection_dicts(
#     nuc_dels: list[str], aa_changes: dict[str, dict[int, dict]]
# ) -> (list[dict], dict[int, dict]):
#     range_list = []
#     range_dict = {}
#     list_of_nuc_pos = set()
#     for gene in aa_changes:
#         for aa_pos in aa_changes[gene]:
#             list_of_nuc_pos.update(aa_changes[gene][aa_pos]["nucPos"])
#     for nuc_del in nuc_dels:
#         # Get start and end inclusive of deletion range
#         if "-" not in nuc_del:
#             start = end = int(nuc_del)
#         else:
#             start, end = map(int, nuc_del.split("-"))
#         itemize = list_of_nuc_pos.intersection(range(start, end + 1)) != set()
#         if itemize:
#             for pos in range(start, end + 1):
#                 range_dict[pos] = {"alt": "-"}
#         else:
#             range_list.append({"start": start, "end": end})
#     return range_list, range_dict


# ranges, items = delection_dicts(
#     xbb["nucDeletions"],
#     aa_changes(
#         xbb["aaSubstitutions"],
#         xbb["aaDeletions"],
#         xbb["nucSubstitutions"],
#         xbb["nucDeletions"],
#     ),
# )

# %%
import json

# d: dict = json.load(open("data/pango-consensus-sequences_summary.json"))
d: dict = json.load(open("test.json"))

# %%
# Frameshifts can be split into two mutations: start and end
# Then treated like any other AA mutations


# %%
# Add new fields for each lineage

# Columns to add from d to dout
pass_through = [
    "lineage",
    "unaliased",
    "parent",
    "children",
    "nextstrainClade",
    "frameShifts",
    "designationDate",
]

dout = {
    lineage: {k: v for k, v in d[lineage].items() if k in pass_through}
    for lineage in d
}
print(dout)
# %%
for lineage, val in d.items():
    dout[lineage]["mutations"] = {}
    for ref, data in val["mutations"].items():
        nucs = {
            **nuc_sub_dict(data["nucSubstitutions"]),
            **nuc_sub_dict(data["nucDeletions"]),
        }
        df = {
            "nuc": nucs,
            "aa": aa_changes(
                data["aaSubstitutions"],
                data["aaDeletions"],
                nucs,
            ),
            "frameshifts": massage_frameshifts(val["frameShifts"], nucs),
        }
        # Add constructed dict to dout
        dout[lineage]["mutations"][ref] = df

# %%

# Add some dummy annotations
# Two types of annotation:
# - "nuc" for nucleotide substitutions
# - "aa" for amino acid substitutions
# Can embed these right in the "nuc" and "aa" dicts

# Add dummy annotation to BA.2.75
annotations = {
    "BQ.1": {
        "byNuc": {
            "26858": "Reverted relative to 21L",
            "27259": "Reverted relative to 21L",
            "27382": "Reverted relative to 21L",
            "27383": "Reverted relative to 21L",
            "27384": "Reverted relative to 21L",
            "28311": "This mutation causes ORG9b:P10S, but the adjacent mutation changes to P10F",
            "28312": "This mutation on top of 28311, change ORf9b:P10S to P10F",
        },
        # keys here can be changed to whatever matches most easily to the main file
        "byAA": {
            "S:24": "Deletion removes last 2 bases from AA 24, entirety of 25/26, first base of 27",
            "S:25": "Deletion removes last 2 bases from AA 24, entirety of 25/26, first base of 27",
            "S:26": "Deletion removes last 2 bases from AA 24, entirety of 25/26, first base of 27",
            "S:27": "This change is a consequence of the prior deletion",
            "S:493": "Note mutation S:Q493R (in 21L) is reverted here",
            "ORF6:61": "Note mutation in ORF6:D61L (in 21L) is reverted here",
            # this mutation also affects ORF9b but for a note, just listing it once should be enough..?!
            "N:31": "Deletion in-frame for ORF9b; in N, deletes last 2 bases of AA 30, entirety of 31/32, first base of 33",
            "N:32": "Deletion in-frame for ORF9b; in N, deletes last 2 bases of AA 30, entirety of 31/32, first base of 33",
            "N:33": "Deletion in-frame for ORF9b; in N, deletes last 2 bases of AA 30, entirety of 31/32, first base of 33",
        },
    }
}

for lineage, data in annotations.items():
    for ref, mut_dicts in dout[lineage]["mutations"].items():
        print(lineage, ref)
        for pos in data["byNuc"]:
            if int(pos) in mut_dicts["nuc"]:
                dout[lineage]["mutations"][ref]["nuc"][int(pos)][
                    "annotation"
                ] = data["byNuc"][pos]
            else:
                print(f"Warning: {pos} not in {lineage} nuc substitutions")
        for pos in data["byAA"]:
            gene, aa_pos = pos.split(":")
            aa_pos = int(aa_pos)
            if gene in mut_dicts["aa"] and aa_pos in mut_dicts["aa"][gene]:
                dout[lineage]["mutations"][ref]["aa"][gene][aa_pos][
                    "annotation"
                ] = data["byAA"][pos]
            else:
                print(f"Warning: {pos} not in {lineage} aa substitutions")

# %%
json.dump(
    dout,
    open("covariants.json", "w"),
    indent=1,
)
# %%
