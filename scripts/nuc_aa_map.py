# - Mapping of nucleotide and amino acid substitutions for improved display on coVariants.org
#   - One field for each direction
#   - `nucToAa`: For each nuc that is substituted, deleted or inserted, insert corresponding aa changes, if any
#   - `aaToNuc`: For each aa that is substituted, deleted or inserted, insert corresponding nuc changes

# For each aa mutation, find the corresponding codon, then pull all nuc mutations for that codon

# %%

# Get map from simple genemap config


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


# %%
# print(nucs_from_aa("ORF9b:T5I"))
# print(nucs_from_aa("S:D614G"))
# print(nucs_from_aa("S:A403P"))
# # %%

# test_dict = {
#     "nucSubstitutions": ["C8782T", "C17747T", "A17858G", "C18060T", "T28144C"],
#     "aaSubstitutions": ["ORF1b:P1427L", "ORF1b:Y1464C", "ORF8:L84S"],
# }

# %%


def nuc_pos_from_nuc(nuc: str) -> int:
    # Takes nuc mut and returns corresponding nuc pos
    # nuc_pos_from_nuc("T28144C") -> 28144
    return int(nuc[1:-1])


import json
from typing import List


def responsible_nucs(aa: str, nuc_subs: List[str]) -> List[str]:
    # Takes aa mut and returns list of nuc muts in corresponding codon
    # responsible_nucs("ORF8:L84S") -> ["T28144C"]
    nucs = nucs_from_aa(aa)
    retval = [nuc for nuc in nuc_subs if nuc_pos_from_nuc(nuc) in nucs]
    return retval


# responsible_nucs("ORF8:L84S", test_dict["nucSubstitutions"])

# # %%
# for aa in test_dict["aaSubstitutions"]:
#     print(aa, responsible_nucs(aa, test_dict["nucSubstitutions"]))

# # %%
# xbb = {
#     "nucSubstitutions": [
#         "A1G",
#         "C44T",
#         "C241T",
#         "A405G",
#         "T670G",
#         "C2790T",
#         "C3037T",
#         "G4184A",
#         "C4321T",
#         "C9344T",
#         "A9424G",
#         "C9534T",
#         "C9866T",
#         "C10029T",
#         "C10198T",
#         "G10447A",
#         "C10449A",
#         "C12880T",
#         "C14408T",
#         "G15451A",
#         "C15714T",
#         "C15738T",
#         "T15939C",
#         "T16342C",
#         "C17410T",
#         "T17859C",
#         "A18163G",
#         "A19326G",
#         "C19955T",
#         "A20055G",
#         "C21618T",
#         "T21810C",
#         "G21987A",
#         "C22000A",
#         "C22109G",
#         "T22200A",
#         "G22577C",
#         "G22578A",
#         "G22599C",
#         "C22664A",
#         "C22674T",
#         "T22679C",
#         "C22686T",
#         "A22688G",
#         "G22775A",
#         "A22786C",
#         "G22813T",
#         "T22882G",
#         "G22895C",
#         "T22896C",
#         "G22898A",
#         "T22942G",
#         "G22992A",
#         "C22995A",
#         "A23013C",
#         "T23019C",
#         "T23031C",
#         "A23055G",
#         "A23063T",
#         "T23075C",
#         "A23403G",
#         "C23525T",
#         "T23599G",
#         "C23604A",
#         "C23854A",
#         "G23948T",
#         "A24424T",
#         "T24469A",
#         "C25000T",
#         "C25416T",
#         "C25584T",
#         "C26060T",
#         "C26270T",
#         "A26275G",
#         "C26577G",
#         "G26709A",
#         "C26858T",
#         "A27259C",
#         "G27382C",
#         "A27383T",
#         "T27384C",
#         "C27807T",
#         "A28271T",
#         "C28311T",
#         "G28881A",
#         "G28882A",
#         "G28883C",
#         "A29510C",
#     ],
#     "aaSubstitutions": [
#         "E:T9I",
#         "E:T11A",
#         "M:Q19E",
#         "M:A63T",
#         "N:P13L",
#         "N:R203K",
#         "N:G204R",
#         "N:S413R",
#         "ORF1a:K47R",
#         "ORF1a:S135R",
#         "ORF1a:T842I",
#         "ORF1a:G1307S",
#         "ORF1a:L3027F",
#         "ORF1a:T3090I",
#         "ORF1a:L3201F",
#         "ORF1a:T3255I",
#         "ORF1a:P3395H",
#         "ORF1b:P314L",
#         "ORF1b:G662S",
#         "ORF1b:S959P",
#         "ORF1b:R1315C",
#         "ORF1b:I1566V",
#         "ORF1b:T2163I",
#         "ORF3a:T223I",
#         "ORF6:D61L",
#         "ORF9b:P10S",
#         "S:T19I",
#         "S:A27S",
#         "S:V83A",
#         "S:G142D",
#         "S:H146Q",
#         "S:Q183E",
#         "S:V213E",
#         "S:G339H",
#         "S:R346T",
#         "S:L368I",
#         "S:S371F",
#         "S:S373P",
#         "S:S375F",
#         "S:T376A",
#         "S:D405N",
#         "S:R408S",
#         "S:K417N",
#         "S:N440K",
#         "S:V445P",
#         "S:G446S",
#         "S:N460K",
#         "S:S477N",
#         "S:T478K",
#         "S:E484A",
#         "S:F486S",
#         "S:F490S",
#         "S:Q498R",
#         "S:N501Y",
#         "S:Y505H",
#         "S:D614G",
#         "S:H655Y",
#         "S:N679K",
#         "S:P681H",
#         "S:N764K",
#         "S:D796Y",
#         "S:Q954H",
#         "S:N969K",
#     ],
#     "nucDeletions": [
#         "11288-11296",
#         "21633-21641",
#         "21992-21994",
#         "28362-28370",
#         "29734-29758",
#     ],
#     "aaDeletions": [
#         "N:E31-",
#         "N:R32-",
#         "N:S33-",
#         "ORF1a:S3675-",
#         "ORF1a:G3676-",
#         "ORF1a:F3677-",
#         "ORF9b:E27-",
#         "ORF9b:N28-",
#         "ORF9b:A29-",
#         "S:L24-",
#         "S:P25-",
#         "S:P26-",
#         "S:Y144-",
#     ],
# }

# for aa in xbb["aaSubstitutions"]:
#     print(aa, responsible_nucs(aa, xbb["nucSubstitutions"]))
# %%

# Some aa substitutions are caused by nuc deletions
# Need to also check if a deletion may be responsible for the aa substitution
# Let's treat deletions as ranges, as they are usually single events
# An aa may be affected if the deletion range overlaps with the aa range


def nuc_del_overlap(nuc_pos: int, del_range: str) -> bool:
    """
    Check if a nuc position overlaps with a deletion range
    """
    if not del_range:
        print("Empty deletion range")
        raise ValueError
    if "-" not in del_range:
        del_start = del_end = int(del_range)
    else:
        del_start, del_end = map(int, del_range.split("-"))
    return del_start <= nuc_pos <= del_end


def responsible_dels(aa: str, dels: List[str]) -> List[str]:
    """
    Given an amino acid substitution, find the deletions that may be responsible for it
    """
    # Find nuc positions that may be responsible for the aa substitution
    nucs = nucs_from_aa(aa)

    retval = [
        del_range
        for del_range in dels
        if any(nuc_del_overlap(nuc, del_range) for nuc in nucs)
    ]

    return retval

# for aa in xbb["aaSubstitutions"]:
#     print(aa, responsible_dels(aa, xbb["nucDeletions"]))

# %%
# Attribute AA deletions to nuc deletions

def nuc_del_causing_aa_del(aa_del: str, nuc_dels: List[str]) -> List[str]:
    """
    Given an amino acid deletion, find the deletions that may be responsible for it
    """
    # Find nuc positions that may be responsible for the aa substitution
    nucs = nucs_from_aa(aa_del)

    retval = [
        del_range
        for del_range in nuc_dels
        if any(nuc_del_overlap(nuc, del_range) for nuc in nucs)
    ]

    return retval

# for aa in xbb["aaDeletions"]:
#     print(aa, nuc_del_causing_aa_del(aa, xbb["nucDeletions"]))
# # %%

# # Produce structured output
# # For each aaSubstition, have a dict of responsibleNucSubstitutions and responsibleNucDeletions

# aaSubstitutionsToNuc = {}
# for aa in xbb["aaSubstitutions"]:
#     aaSubstitutionsToNuc[aa] = {
#         "responsibleNucSubstitutions": responsible_nucs(aa, xbb["nucSubstitutions"]),
#         "responsibleNucDeletions": responsible_dels(aa, xbb["nucDeletions"]),
#     }

# aaDeletionsToNuc = {}
# for aa in xbb["aaDeletions"]:
#     aaDeletionsToNuc[aa] = {
#         "responsibleNucDeletions": nuc_del_causing_aa_del(aa, xbb["nucDeletions"]),
#     }

# print(aaSubstitutionsToNuc)
# print(aaDeletionsToNuc)
    
# # %%

# # Can also invert: for each nucSubstitution, show the aaSubstitutions that it may be responsible for

# nucSubstitutionsToAa = { nuc: set() for nuc in xbb["nucSubstitutions"] }

# for aa in aaSubstitutionsToNuc:
#     for nuc in aaSubstitutionsToNuc[aa]["responsibleNucSubstitutions"]:
#         nucSubstitutionsToAa[nuc].add(aa) 

# print(nucSubstitutionsToAa)

# # %%
# nucDeletionsToAa = { nuc: set() for nuc in xbb["nucDeletions"] }

# for aa in aaSubstitutionsToNuc:
#     for nuc in aaSubstitutionsToNuc[aa]["responsibleNucDeletions"]:
#         nucDeletionsToAa[nuc].add(aa)

# for aa in aaDeletionsToNuc:
#     for nuc in aaDeletionsToNuc[aa]["responsibleNucDeletions"]:
#         nucDeletionsToAa[nuc].add(aa)

# print(nucDeletionsToAa)

# %%

# Need special treatment for changes from parent other than reference
# Basically need to find the differences now to the parent
# This gives a new set of nucSubstitutions, nucDeletions, aaSubstitutions, aaDeletions
# In principle, nothing changes, causality should still hold

# Discussed with Emma, we shouldn't just have nucSubstitutionsNew (relative to Pango parent)
# We should also do nucSubstitutionsSinceClade (relative to clade parent)

# These differences create the new category of "reversion", if a mutation is present in parent
# but not (or different) in child

# These are noteworthy evolutionary events, but require special treatment

# Read in json
# Add new fields for each lineage
# Write out json
import json
d: dict = json.load(open("data/pango-consensus-sequences_summary.json"))

#%%
# Add new fields for each lineage

for lineage, data in d.items():
    data["aaSubstitutionsToNuc"] = {}
    for aa in data["aaSubstitutions"]:
        data["aaSubstitutionsToNuc"][aa] = {
            "responsibleNucSubstitutions": responsible_nucs(aa, data["nucSubstitutions"]),
            "responsibleNucDeletions": responsible_dels(aa, data["nucDeletions"]),
        }

    data["aaDeletionsToNuc"] = {}
    for aa in data["aaDeletions"]:
        data["aaDeletionsToNuc"][aa] = {
            "responsibleNucDeletions": nuc_del_causing_aa_del(aa, data["nucDeletions"]),
        }
    

    data["nucSubstitutionsToAa"] = { nuc: set() for nuc in data["nucSubstitutions"] }
    for aa in data["aaSubstitutionsToNuc"]:
        for nuc in data["aaSubstitutionsToNuc"][aa]["responsibleNucSubstitutions"]:
            data["nucSubstitutionsToAa"][nuc].add(aa) 

    data["nucDeletionsToAa"] = { nuc: set() for nuc in data["nucDeletions"] }

    for aa in data["aaSubstitutionsToNuc"]:
        for nuc in data["aaSubstitutionsToNuc"][aa]["responsibleNucDeletions"]:
            data["nucDeletionsToAa"][nuc].add(aa)
        
    for aa in data["aaDeletionsToNuc"]:
        for nuc in data["aaDeletionsToNuc"][aa]["responsibleNucDeletions"]:
            data["nucDeletionsToAa"][nuc].add(aa)

# From https://stackoverflow.com/a/8230505/7483211
class SetEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        return json.JSONEncoder.default(self, obj)

json.dump(d, open("data/pango-consensus-sequences_summary_new.json", "w"), indent=1, cls=SetEncoder)
# %%
