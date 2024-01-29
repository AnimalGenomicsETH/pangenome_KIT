
if [[ $1 == S* ]] ;
then
  samtools view -h --reference /PacBio_CLR/${1}.cram 6:69999582-70223136 | samtools fasta | seqtk seq -L 10000 > ${1}.fa
else
  samtools view -h --reference /REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa /PacBio_CLR/${1}.mm2.bam 6:69999582-70223136 | samtools fasta | seqtk seq -L 10000 > ${1}.fa
fi

minimap2 -ax map-pb -r 2k -P --end-bonus 1000 -Y ${1}.fa SxG.repeat.fa | samtools sort --write-index -o ${1}.bam
samtools view -e 'rlen>13000' ${1}.bam | cut -f 3 | sort | uniq -c | sort -k1,1nr | awk -v S=${1} -v L=$(cat ${1}.fa | wc -l)  '{++a[$1]} END {for (k in a) {print S,L,k,a[k]}}' > ${1}.count
