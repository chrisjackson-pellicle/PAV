# BLAST Database: order_genomes_blastdb

## Database Information
- **Database Name**: order_genomes_blastdb
- **Type**: Nucleotide database
- **Source**: extracted_features_location_extract.fasta
- **Sequences**: 90,632 features
- **Total Bases**: 61,509,485
- **Longest Sequence**: 8,052 bases
- **Creation Date**: August 18, 2025

## Database Files
- `order_genomes_blastdb.ndb` - Database index
- `order_genomes_blastdb.nhr` - Header file
- `order_genomes_blastdb.nin` - Index file
- `order_genomes_blastdb.njs` - Journal file
- `order_genomes_blastdb.not` - Note file
- `order_genomes_blastdb.nsq` - Sequence file
- `order_genomes_blastdb.ntf` - Taxonomy file
- `order_genomes_blastdb.nto` - Taxonomy offset file

## Source Data
This database was created from `extracted_features_location_extract.fasta`, which contains functional features extracted from GenBank files in the order_genomes directory. The extraction process uses BioPython's `location.extract()` method for improved sequence extraction, especially for complex features with multiple exons or features spanning sequence joins.

The extraction process excludes:
- source features
- LSC (Large Single Copy) regions
- SSC (Small Single Copy) regions
- repeat_region features
- misc_feature features
- intron features
- gene features (only functional features like CDS, tRNA, rRNA are included)
- inverted repeat features

## Sequence Extraction Improvements
This version uses BioPython's `location.extract()` method instead of manual concatenation, which ensures:
- Correct exon ordering for complex features
- Proper handling of features spanning sequence joins
- Accurate sequence extraction for compound locations
- Better handling of features with different strands

## Usage
To use this database for BLAST searches:

```bash
blastn -db plastid_annotation_validator/data/order_genomes_blastdb/order_genomes_blastdb -query your_sequence.fasta -out results.txt
```

## FASTA Header Format
The sequences in this database have headers in the format:
`Order_Family_genus_Species_accession_feature_type_[gene_name]`

For example:
`Acorales_Acoraceae_Acorus_tatarinowii_NC_045294.1_CDS_[rps12]`
