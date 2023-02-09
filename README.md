# Consensus sequences for each Pango lineage

This repository contains semi-automatically generated prototype sequences for each Pango lineage. These sequences are not real sequences in databases but are algorithmically constructed consensus sequences that try to represent the common ancestor sequence of that lineage. They are based on the sequences designated in the cov-lineages/pango-designation repository. There is some manual curation involved to overwrite erroneous sites - errors happen when a lineage's designated sequences have dropout or reversions. The algorithm used to create these sequences has a high threshold to allow reversions, almost all sequences need to be reverted, otherwise it's assumed the reversions are an artefact.

Please be aware that due to the semi-automatic generation of these synthetic sequences, they can contain errors.

The data in this repo is automatically updated every day.

For lineages that are derived from BA.2*, BA.4* and BA.5*, there has been significant curation/overwriting, but for previous lineages the amount of curation is very limited. So do expect errors.

If you find errors, please open an issue here. The same is true if you have ideas what other data about lineages you would like to be included here.

The sequences contained here are the ones used in Nextclade reference trees and produced by code contained in the nextclade_data_workflows/sars-cov-2 repository.

## Contents

The repository contains:

- Zstd compressed fasta files of the prototype nucletide sequences and amino acid sequences of the translations of common genes. They are available in the `data/` folder of this repository. You can easily decompress them after downloading with `zstd -d data/pango-consensus-sequences_nuc.fasta.zst`. You can pick a single sequence of interest using `seqkit grep -r -p "B.1.1.7" pango_lineages.fasta`.:
  - Nucleotide sequences are at: `data/pango-consensus-sequences_nuc.fasta.zst`
  - Translations are at: `data/pango-consensus-sequences_S.fasta.zst` etc.
- A JSON file dictionary with Pango lineage names as keys, e.g. `BQ.1` and dicts with the following values per sequence:
  - `lineage`: Same as key, so that one can treat dict as list and not need to use the key
  - `unaliased`: Unaliased pango lineage name, e.g. `B.1.1.529.5.3.1.1.1.1`
  - `parent`: The direct parent lineage, e.g. `BE.1.1.1` (or the first parent for which there is a consensus sequence if the direct parent doesn't have one)
  - `children`: Array of child lineages that are present in the dataset, e.g. `["BQ.1.1", "BQ.1.2"]`
  - `nextstrainClade`: The Nextstrain clade this lineage belongs to, e.g. `22E`
  - `nucSubstitutions`: Array of nucleotide substitution strings, e.g. `["C241T","C8782T"]`
  - `aaSubstitutions`: Array of amino acid substitution strings, e.g. `["S:L452R", ORF1b:P1427"]`
  - `nucDeletions`: Array of nucleotide deletion strings, e.g. `["21983-21991"]`
  - `aaDeletions`: Array of amino acid deletion strings, e.g. `["S:Y144-"]`
  - `nucSubstitutionsNew`: Like `nucSubstitutions`, but only includes substitutions that are new with respect to the parent lineage (can be considered as _defining_ mutations of this lineage). Same format as `nucSubstitutions`.
  - `aaSubstitutionsNew`: Like `nucSubstitutionsNew` but for amino acid substitutions. Same format as `aaSubstitutions`.
  - `nucDeletionsNew`: Like `nucSubstitutions`, but for nucleotide deletions. Same format as `nucDeletions`.
  - `aaDeletionsNew`: Like `nucDelertionsNew`, but for amino acid deletions. Same format as `aaDeletions`.
  - `nucSubstitutionsReverted`: Like `nucSubstitutions`, but only includes substitutions that are present in the parent but not in this lineage
  - `aaSubstitutionsReverted`: Like `nucSubstitutionsReverted`, but for amino acid substitutions
  - `nucDeletionsReverted`: Like `nucSubstitutionsReverted`, but for nucleotide deletions. Same format as `nucDeletions`.
  - `aaDeletionsReverted`: Like `nucDelertionsReverted`, but for amino acid deletions. Same format as `aaDeletions`.
  - `frameShifts`: Array of frame shift strings, e.g. `["S:142-143"]`
  - `designationDate`: Designation date of this lineage, e.g. `2022-12-01` or `null` if not available
  These files are useful if you want to check which mutations are contained and don't want to align the fasta files themselves. They are also useful to see what changed between different versions of the sequences, as they can be diffed.

## Caveats

- Not all designated lineages are included in this repo. Some (in particular early) pandemic lineages have too few or no designated sequences available for the algorithm to be confident.
- Consensus sequences are only published if there are more than 3 sequences of a lineage available as open data (Genbank/RKI/COG-UK)
- Insertions are currently never included in consensus sequences. They may or may not be added in the future.
- Currently, when one position is mutated twice, e.g. `486F->S->L` this gives rise to both a reversion and a new substitution. This might be changed in the future.

If you need to be sure that a sequence is correct, e.g. when you're creating a Spike protein for an experiment, please double check using the pango designation issue (if such an issue exists) and the annotation on the Usher tree - which is independently curated by @AngieHinrichs. Also, see <https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/utils/otto/sarscov2phylo/pango.clade-mutations.tsv> for the paths extracted from the Usher tree for each lineage.
