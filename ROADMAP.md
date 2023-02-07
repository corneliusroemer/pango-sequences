# Roadmap

- Auto publish unchecked sequences `nightly` branch
- At least once a week, go and review sequences, do overwrite and release to `main` branch
- Execute TSV generation on `nightly` sequences

## Workflow

- Cronjob publishes nuc sequences to AWS s3 after full run
  - Sort alphabetically for stability
  - Runs nextalign to make translations and also uploads
  - Path: ? (TODO: get data.nextstrain.org path)
  - Temporary path:
    - s3://nextstrain-data/pango-consensus-sequences_nuc.fasta
    - s3://nextstrain-data/pango-consensus-sequences_S.fasta
    - ...
- Github repo or scicore takes these and post processes them:
  - Implement the following with snakemake and python scripts
  - Zstd compress and commit to data folder with same path as above
  - Produce JSON analysis file for better diffs, with indent of just 1 for smaller file size
    - Sorted alphabetically
    - Dict with key `B.1.1.7`
    - Properties:
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
      - Designation date
      - Sequence last updated


## Future work

- Sort topologically by unaliased name
- Sort mutations numerically
- Comparison to Usher
- Add insertions