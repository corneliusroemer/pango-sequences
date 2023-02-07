# Dev docs

## Run SARS-CoV-2 nextclade_workflow

```bash
gh repo clone neherlab/nextclade_data_workflows ~/code
cd ~/code/nextclade_data_workflows/sars-cov-2
snakemake --profile profiles/clades all
```

Make sure this runs every day using a cronjob on scicore.

## Updating data on the nightly branch

```bash
gh repo clone corneliusroemer/pango-sequences ~/code
cd ~/code/pango-sequences
git checkout nightly
cp ~/code/nextclade_data_workflows/sars-cov-2/pre-processed/synthetic.fasta data/pango_consensus_sequences.fasta
snakemake -c1
git add --all
git commit -m "Nightly data update: $(date)"
```

## Fixing wrong bases, frameshifts etc

The file to do overwrites in is at: <https://github.com/neherlab/nextclade_data_workflows/blob/master/sars-cov-2/profiles/clades/lineage_overwrite.tsv>

