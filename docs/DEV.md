# Dev docs

## Run SARS-CoV-2 nextclade_workflow

```bash
gh repo clone neherlab/nextclade_data_workflows ~/code
git checkout pango-consensus-only
cd ~/code/nextclade_data_workflows/sars-cov-2
snakemake --profile profiles/clades all
```

Make sure this runs every day using a cronjob on scicore.

## Updating data in this repo on server

```bash
gh repo clone corneliusroemer/pango-sequences ~/code
cd ~/code/pango-sequences
snakemake -c5
```

## Local development without deploying

```bash
snakemake -c10 -F --ri -p --config local=True test=True
```

`local=True` will download sequences from server
`test=True` will not commit and publish changes

## Fixing wrong bases, frameshifts etc

The file to do overwrites in is at: <https://github.com/neherlab/nextclade_data_workflows/blob/master/sars-cov-2/profiles/clades/lineage_overwrite.tsv>
