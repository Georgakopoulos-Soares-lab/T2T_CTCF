#!/bin/bash

agat_convert_sp_gxf2gxf.pl --gff ../data/chm13/chm13.draft_v2.0.gene_annotation.gff3 --out ../data/chm13/chm13_standardized.gff3
agat_sp_add_introns.pl --gff ../data/chm13/chm13_standardized.gff3 --out ../data/chm13/chm13_gencode_full.gff3