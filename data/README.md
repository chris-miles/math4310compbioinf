# Course Data

## phiX174 reference sequence

[phix174-NC_001422.1.fasta](phix174-NC_001422.1.fasta) is the complete 5,386-base Escherichia phage phiX174 reference genome, accession **NC_001422.1**, downloaded from [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1) on 2026-09-22. It is small enough for word counting, synthetic-read assembly experiments, and simple sequence-model checks on a laptop. It is a reference sequence, not observed sequencing reads or a CpG-island benchmark.

The FASTA preserves the downloaded sequence, header, and line wrapping; the extra blank line at the end was removed. Source: [NCBI EFetch](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text). File SHA-256: `a0d9bfd990c62c0908392b962bb65b574466ad9fba3777681e40f7fb3debf3a2`. The download was checked for one record, accession, length, and an A/C/G/T-only alphabet.

The biological genome is circular, but the FASTA writes it from one chosen starting position. The introductory examples count linear substrings and omit words crossing the end/start boundary. For 3-mers, that gives 5,384 total occurrences; a circular count would give 5,386. Neither count is the number of distinct words.

NCBI makes molecular sequence data available for scientific use; see its [data policies](https://www.ncbi.nlm.nih.gov/home/about/policies/). Credit the accession and source when reusing this file. This folder contains no human-subject data.

## Phage G4 reference sequence

[g4-NC_001420.2.fasta](g4-NC_001420.2.fasta) is the complete 5,577-base genome of Escherichia phage G4, accession **NC_001420.2**, downloaded from [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/nuccore/NC_001420.2) on 2026-10-24. G4 and phiX174 are related single-stranded DNA phages of E. coli, so their genes give real homologous pairs for global and local alignment. The record notes that it is identical to GenBank V00657 and is still a provisional RefSeq. Gene coordinates quoted in the problem sets come from this record's annotation; for example, gene D is at 1976–2434 and gene F at 2600–3883.

The FASTA preserves the downloaded sequence, header, and line wrapping; the extra blank line at the end was removed. Source: [NCBI EFetch](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001420.2&rettype=fasta&retmode=text). File SHA-256: `f43b70c2fa4fe809db4feb74f5faa27062eeb76c168319a65433c148f8c99bbb`. The download was checked for one record, accession, length, and an A/C/G/T-only alphabet. Like phiX174, the genome is circular and its FASTA starts at one chosen position; the two files do not start at corresponding places.

## Working with course data

Examples read local files and do not download during rendering. Tiny artificial strings stay inline in the notes. Synthetic reads generated from this reference must be labeled synthetic, with the sampling and error assumptions stated.
