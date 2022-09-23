genes = [
    "S",
    "ORF1a",
    "ORF1b",
    "ORF3a",
    "E",
    "M",
    "N",
    "ORF6",
    "ORF7a",
    "ORF7b",
    "ORF8",
    "ORF9b",
]


rule all:
    input:
        expand("data/translations/pango_consensus_gene_{gene}.fasta.zst", gene=genes),
        "build/nextclade_subset.tsv",
        "build/nextclade_subset_unaliased.tsv",
        "build/nextclade_subset_unaliased_relative.tsv",


rule run_nextclade:
    input:
        "data/pango_consensus_sequences.fasta.zst",
    output:
        tsv="build/nextclade.tsv",
    shell:
        "nextclade run --in-order -d sars-cov-2 --output-all build {input}"


rule compress_translations:
    input:
        "build/nextclade_gene_{gene}.translation.fasta",
    output:
        "data/translations/pango_consensus_gene_{gene}.fasta.zst",
    shell:
        "zstd -T0 -19 {input} -o {output}"


rule subset_tsv:
    input:
        "build/nextclade.tsv",
    output:
        "build/nextclade_subset.tsv",
    params:
        fields="seqName,clade,totalSubstitutions,totalDeletions,totalInsertions,"
        + "totalFrameShifts,totalAminoacidSubstitutions,totalAminoacidDeletions,"
        + "frameShifts,aaDeletions,deletions,substitutions,aaSubstitutions",
    shell:
        "tsv-select -H -f {params.fields} {input} > {output}"


# Use pango_aliasor package
rule add_unaliased_name:
    input:
        "build/nextclade_subset.tsv",
    output:
        "build/nextclade_subset_unaliased.tsv",
    shell:
        "python scripts/add_unaliased.py {input} {output}"


# Add mutations that are different from the parent lineage
# Add as private nuc, private del, private reversions
# private aa, private aadel, private aareversions
# Make a script that infers these
rule add_relative_mutations:
    input:
        "build/nextclade_subset_unaliased.tsv",
    output:
        "build/nextclade_subset_unaliased_relative.tsv",
    shell:
        "python scripts/add_relative_mutations.py {input} {output}"
