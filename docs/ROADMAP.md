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
    - Potential future additions:
      - Pango notes
      - Earliest designated sequences included in first designation and earliest few by submission date, including lab
      - Breakpoints and donors of recombinants

## Future work

- Sort topologically by unaliased name
- Sort mutations numerically
- Comparison to Usher
- Add insertions
- Commits of changes of summary
- How to deal with hypothetical parents like B.1.640 that don't have consensus sequences?
  - Give them real parent even if higher up
- What to do with pango sequences that don't have a consensus sequence? (there are ~100)
- Deal with mutations that are not reverted like S:F486S -> S:F486P
  - List genotypes not substitutions in nucMutations?
  - And make reversions contain from to?
  - Maybe have extra fields for further substitutions?
