#!/usr/bin/env bash
# VERSION: 0.6.0

# Slightly modified version of mixcr v3.0.14, using a custom
# repseqio/library/human/TRA_TRD.json file which does not combine
# TRA and TRD V genes, except with hybrid genes such as TRAV29DV5 
export MIXCR="${HOME}/mixcr_mod/mixcr"

# migec v1.2.9
export MIGEC=$(which migec)

export SPECIES="human"

for p in parallel $MIGEC $MIXCR; do
    command -v $p >/dev/null 2>&1 || { echo >&2 "error: cannot find $p in path."; exit 1;  }
done

test -f procfile || { echo >&2 "error: procfile missing in current directory."; exit 1; }

set -ex

# Before running, need to break samplesheet into parts
./src/split_samplesheet.py -o multi RA_Rodriguez.tsv

# Assembling UMI-collapsed fastqs with migec 
ls multi | parallel --bar -P procfile 'mkdir -p checkout/{/.} && $MIGEC CheckoutBatch -cute --skip-undef multi/{} checkout/{/.}'
ls multi | parallel -P procfile --bar 'mkdir -p histogram/{/.} && $MIGEC Histogram checkout/{/.} histogram/{/.}'
ls multi | parallel -j4 --bar 'mkdir -p assembly/{/.} && $MIGEC AssembleBatch -c checkout/{/.} histogram/{/.} assembly/{/.}'

# Creating the exportClones mixcr output format, following the guidelines
# for "Assemble full TCR/Ig receptor sequences" as described here:
# https://mixcr.readthedocs.io/en/develop/assembleContigs.html
mkdir -p alignments
awk 'FNR>1{print $1","$5","$6}' assembly/*/assemble.log.txt | parallel --bar -j1 --colsep ',' '$MIXCR align -OreadsLayout=Collinear --species $SPECIES -j alignments/{1}.json {2} {3} alignments/{1}.vdjca || echo {/.} >> FAILED1'
ls alignments | grep vdjca | parallel -P procfile --bar '$MIXCR assemble -OcloneClusteringParameters=null --write-alignments -OseparateByV=true -OseparateByJ=true -OseparateByC=true -j alignments/{/.}.json alignments/{/.}.vdjca alignments/{/.}.clna || echo {/.} >> FAILED2'
ls alignments | grep vdjca | parallel -P procfile --bar '$MIXCR assembleContigs -j alignments/{/.}.json alignments/{/.}.clna alignments/{/.}.clns || echo {/.} >> FAILED3'
ls alignments | grep vdjca | parallel -P procfile --bar '$MIXCR exportClones --filter-out-of-frames --filter-stops alignments/{/.}.clns alignments/{/.}.tsv || echo {/.} >> FAILED4'
ls alignments | grep tsv | parallel --bar './split_mixcr_file.py alignments/{} -o exported_clones'

# Getting maximum length sequences using exportAlignments. This file format is needed
# to calculate SHM (somatic hypermutation) for each sequence, as well as the isotype
# indexes such as CSR. 
fd vdjca alignments | parallel --bar '$MIXCR exportAlignments -p fullImputed {} alignments/{/.}.tsv'
fd tsv alignments | parallel --bar './mixcr_imputed_hotfix.py -r human_reference.tsv {} -o {}'
fd tsv alignments | parallel --bar './split_exported_alignments_file.py {} -o exported_alignments_imputed'

# Assembling into a simpler format, using the V, J, and C genes, and collapsing
# sequences from CDR1 -> FR4. Primer sequences are in FR1, so all sequences
# should cover this region. 
mkdir full_assembled
fd IGH exported_alignments_imputed | parallel --bar './convert_exported_alignments_to_assembled.py {} --region-start CDR1 --is-igh -o full_assembled/{/.}.csv --species $SPECIES'
fd "IGK|IGL|TRA|TRB|TRD|TRG" exported_alignments_imputed | parallel --bar './convert_exported_alignments_to_assembled.py {} --region-start CDR1 -o full_assembled/{/.}.csv --species $SPECIES'

# shm_calculator is a custom program that realigns the full sequence to the IMGT V reference sequences
# (using the already called V gene to reduce computation) and reports the number of mismatches and alignment length
cd full_assembled
ls | rg IGH | parallel --bar 'shm_calculator -f IMGT_202113-2/IGHV.fasta -i7 -v1 --has-headers {} > {/.}.tmp && xsv cat columns {} {/.}.tmp | sponge {} && rm {/.}.tmp'
ls | rg IGK | parallel --bar 'shm_calculator -f IMGT_202113-2/IGKV.fasta -i7 -v1 --has-headers {} > {/.}.tmp && xsv cat columns {} {/.}.tmp | sponge {} && rm {/.}.tmp'
ls | rg IGL | parallel --bar 'shm_calculator -f IMGT_202113-2/IGLV.fasta -i7 -v1 --has-headers {} > {/.}.tmp && xsv cat columns {} {/.}.tmp | sponge {} && rm {/.}.tmp'
cd ..
