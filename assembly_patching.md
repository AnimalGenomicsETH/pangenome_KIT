

Align the target sequence to get the cut region

```
minimap2 -xasm5 SIM_2.fa blast.seq > 6.paf
```

We now cut either side of the matched region (add a good chunk either side to be safe)

```
samtools faidx SIM_2.fa SIM_2#2#6:1-76110239 > base.fa
samtools faidx SIM_2.fa SIM_2#2#6:76324556-1000000000 >> base.fa
```

Manually need to join those two contigs and separate by some "N"s, so we have a fake "gap" to fill

Extract reads from the region of interest, and then assemble with extra repeat care (`-D 10`)

```
samtools view -h -D HP:<(echo -e "2\n0") --reference ARS-UCD1.2_Btau5.0.1Y.fa SxGCHEF120170207244.mm2.cram 6:69099582-71123136 | samtools fasta -0 SIM_2.SIM_1.fa.gz
hifiasm --primary -t 8 -o SIM -r 3 -a 5 -n 5 -D 10 -l 0 SxG.SIM.fa.gz
```

Now use ragtag to patch in the extra copy into that gap

```
ragtag.py patch base.fa SIM_1.fa --aligner minimap2 --mm2-params '-cx asm5' -t 4 -f 5000 --remove-small -s 10000 -i 100000000000
```
