# Harmonized annotation schema

- `sample_annotation.tsv`: every study- and sample-level field a harmonized dataset should provide,
  with its type and controlled vocabulary.
- `vocabularies.tsv`: allowed values of the local enums.

External vocabularies: **OncoTree** / **NCI Thesaurus** for disease, **HGNC** symbols for genes
(features in exported matrices use HGNC symbols; Entrez IDs may follow in parentheses, e.g.
`MTOR (2475)`).

This is a draft to be revised during the jamboree's data inventory.
