grep -Ff <(for i in publicSamples/fastq/*1.fixed.fastq.gz; do echo -n "$i "; zcat $i | head -n 1; done | grep -v ":" | grep -oP "\w+(?=_)") sample_information.csv | awk '{print $1}'
