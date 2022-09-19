# Consensus sequences for each Pango lineage

This repository contains semi-automatically generated prototype sequences for each Pango lineage. These sequences are not real sequences in databases but are algorithmically constructed consensus sequences that try to represent the common ancestor sequence of that lineage. They are based on the sequences designated in the cov-lineages/pango-designation repository. There is some manual curation involved to overwrite erroneous sites - errors happen when a lineage's designated sequences have dropout or reversions. The algorithm used to create these sequences has a high threshold to allow reversions, almost all sequences need to be reverted, otherwise it's assumed the reversions are an artefact.

Please be aware that due to the semi-automatic generation of these synthetic sequences, they can contain errors.

For lineages that are part of BA.2*, BA.4* and BA.5* there has been significant curation/overwriting, but for previous lineages the amount of curation is very limited. So do expect errors.

If you find errors, please open an issue here and I will fix it for everyone.

The sequences contained here are the ones used in Nextclade reference trees and produced by code contained in the nextclade_data_workflows repository.

The repository contains:

- A zst compressed fasta file with all sequences, you can easily decompress with `zstdcat pango_lineages.fasta.zst > pango_lineages.fasta`. And you can pick a sequence of interest using `seqkit grep -r -p "B.1.1.7" pango_lineages.fasta`.
- A tsv file with the following columns: `lineage_name`, `unaliased_lineage_name`, `substitutions`, `aa_substitutions`, `deletions`, `new_substitutions`, `new_aa_substitutions`, `new_deletions`, `new_reversions`, `new_aa_reversions`, `new_undeletions`. These files are useful if you want to check which mutations are contained and don't want to align the fasta files themselves. They are also useful to see what changed between different versions of the sequences, as they can be diffed. The `new_` columns refer to the differences with respect to the parent lineage.

If you need to be sure that a sequence is correct, e.g. when you're creatin a Spike protein for an experiment, please double check using the pango designation issue (if such an issue exists) and the annotation on the Usher tree - which is independently curated by @AngieHinrichs. Also, see <https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/utils/otto/sarscov2phylo/pango.clade-mutations.tsv> for the paths extracted from the Usher tree for each lineage.
