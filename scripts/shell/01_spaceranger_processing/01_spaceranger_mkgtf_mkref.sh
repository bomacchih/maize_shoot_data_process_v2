# Prepare a reference for use with 10x analysis software. Requires a GTF and FASTA.
spaceranger mkgtf Zm-B73-Reference-NAM-5.0.gtf Zm-B73-Reference-NAM-5.0_filtered.gtf \
                   --attribute=gene_biotype:protein_coding \
                   --attribute=gene_biotype:lncRNA \
                   --attribute=gene_biotype:antisense

spaceranger mkref --genome B73_V5 --fasta AGPv5_B73.fa --genes Zm-B73-Reference-NAM-5.0_filtered.gtf