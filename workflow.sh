#!/bin/bash
# workflow.sh
# Purpose: Process VCF and extract DP, SNP/INDEL, mutation type, and combine into a single TSV with proper headers

zcat /data-shared/vcf_examples/luscinia_vars.vcf.gz | grep -v '^##' > no-headers.vcf
IN=no-headers.vcf
<$IN cut -f1-6 > cols1-6.tsv

<$IN awk '{
    if(match($0,/DP=([^;]+)/))
        print substr($0,RSTART+3,RLENGTH-3)
    else
        print "NA"
}' > dp.tsv

<$IN awk '{if($0 ~ /INDEL/) print "INDEL"; else print "SNP"}' > col-type.tsv
<$IN awk -F'\t' '{
    if(length($4)==1 && length($5)==1){
        pair=$4$5
        if(pair=="AG" || pair=="GA" || pair=="CT" || pair=="TC")
            print "Transition"
        else
            print "Transversion"
    } else print "NA"
}' > col-transition.tsv

paste cols1-6.tsv dp.tsv col-type.tsv col-transition.tsv > cols-all.tsv

# Replace first line (header) in cols-all.tsv
sed -i '1s/.*/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tDP\tSNP.INDEL\tMutation_type/' cols-all.tsv

