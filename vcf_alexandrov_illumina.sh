#!/bin/bash

# activate conda env
eval "$(conda shell.bash hook)"
conda activate hrd

# get input vcf files
cd vcf_files

mkdir ../unziped

# remove spaces in vcf file names
for f in *\ *; do mv "$f" "${f// /_}"; done

# unzip the files and remove the .gz files after
for f in *.gz; do
  echo ${f}
  STEM=$(basename "${f}" .gz)
  gunzip -c "${f}" > ../unziped/"${STEM}"
  rm ${STEM}.gz
done

# get names of vcf files to convert 
VAR=$(ls ../unziped)

for i in ${VAR}
  do
    echo "converting vcf: ${i}"
    awk '$1 !~/^##/' ../unziped/${i} > ${i}_vcf_removed_header.txt
    awk -F'\t' -vcols=#CHROM,POS,FILTER,REF,ALT '(NR==1){n=split(cols,cs,",");for(c=1;c<=n;c++){for(i=1;i<=NF;i++)if($(i)==cs[c])ci[c]=i}}{for(i=1;i<=n;i++)printf "%s" FS,$(ci[i]);printf "\n"}' ${i}_vcf_removed_header.txt  > ${i}.txt
    awk '$5 !~ /^[.]/' ${i}.txt > ${i}_filtered.txt
    awk '$3 == "PASS"' ${i}_filtered.txt > ${i}
    
  done
  
rm *.txt
rm -r ../unziped

# start python script
#cd ..
#python3 alexandrov_sig.py
