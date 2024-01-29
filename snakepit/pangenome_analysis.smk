import subprocess

def rename_region(pangenome,fasta,just_name=False):

    if just_name:
        chromosome = pangenome
    elif pangenome in config['graphs']:
        region = config['graphs'][pangenome]['region']
        chromosome = region['ARS'].split(':')[0]
    elif pangenome in config['gggenes']:
        region = config['gggenes'][pangenome]['region']
        chromosome = region['ARS'].split(':')[0]
    try:
        panSN = subprocess.check_output(f"zcat {fasta} | head -n 1 | sed 's/>//'",shell=True,encoding="utf-8",stderr=subprocess.STDOUT).strip()
    except subprocess.CalledProcessError:
        return '' #file doesn't exist yet, but will when we need this
    if 'gzip' in panSN:
        return '<TBD>' # panSN doesn't really exist yet
    if just_name:
        return panSN
    return region['ARS'].replace(chromosome,panSN,1)

def convert_region_to_chromosome(wildcards):
    key = 'gggenes' if wildcards.pangenome in config['gggenes'] else 'graphs'
    return config[key][wildcards.pangenome]['region']['ARS'].split(':')[0]

rule odgi_extract:
    input:
        og = lambda wildcards: expand(rules.pggb_construct.output['og'],chromosome=convert_region_to_chromosome(wildcards)),
        fasta = lambda wildcards: expand(rules.panSN_spec.output[0],chromosome=convert_region_to_chromosome(wildcards))
    output:
        '{pangenome}/pggb.subgraph.{renamed}.gfa'
    params:
        region = lambda wildcards, input: rename_region(wildcards.pangenome,input.fasta),
        rename = lambda wildcards: "sed 's/:[0-9]\+-[0-9]\+//'" if wildcards.renamed == 'renamed' else "cat"
    threads: 2
    resources: 
        mem_mb = 2500
    shell:
        '''
        odgi extract -t {threads} -i {input.og} -r "{params.region}" -o /dev/stdout |\
        odgi sort -t {threads} -p "wc" -i - -o - |\
        odgi view -t {threads} -i - -g |\
        {params.rename} > {output}
        '''

rule node_sizes:
    input:
        expand(rules.odgi_extract.output,renamed='renamed',allow_missing=True)
    output:
        '{pangenome}/subgraph.nodes.csv'
    localrule: True
    shell:
        '''
        awk '$1=="S" {{print $2,length($3)}}' {input} > {output}
        '''

rule odgi_sort:
    input:
        og = rules.pggb_construct.output['og']
    output:
        gfa = 'graphs/{chromosome}.pggb.sorted.gfa',
        og = 'graphs/{chromosome}.pggb.sorted.og'
    threads: 1
    resources:
        mem_mb = 20000
    shell:
        '''
        odgi sort -t {threads} -i {input} -p "wc" -o /dev/stdout | tee {output.og} | odgi view -i - -g > {output.gfa}
        '''

rule find_bubbles:
    input:
        rules.odgi_sort.output['gfa']
    output:
        'graphs/{chromosome}.graph_bubbles.list'
    resources:
        mem_mb = 35000
    shell:
        '''
        gfatk SSC --size 5 {input} | sort -k1,1n | sed 's/^0/1/' > {output}
        '''

rule jaccard_index:
    input:
        graph = rules.odgi_sort.output['og'],
        bubbles = rules.find_bubbles.output,
        fasta = rules.panSN_spec.output[0]
    output:
        'graphs/{chromosome}.graph_bubbles.jaccard.gz'
    params:
        reference_name = lambda wildcards, input: rename_region(wildcards.chromosome,input.fasta,True),
        graph = lambda wildcards, input: Path(input.graph).resolve(),
        _output = lambda wildcards, output: Path(output[0]).resolve(),
        padding = 1000
    threads: 8
    resources:
        mem_mb = 2500
    shell:
        '''
        set +e #some similarities still segfault, but we can ignore them
        echo "path sample1 sample2 jaccard" | pigz > {output}

        odgi position -t {threads} -i {input.graph} -r {params.reference_name} -G <(cat {input.bubbles} | tr ' ' '\\n') | awk -F',' 'NR>1 {{print $1,$4}}' |\

        awk 'FNR==NR {{ a[$1] = $2; next }} {{ print a[$1],a[$2]}}' - {input.bubbles} | awk 'NF==2&&!/NA/' |\
        awk -v OFS='\\t' '$1>0&&$2>0&&$2>$1 {{print "{params.reference_name}",$1-{params.padding},$2+{params.padding};next}} {{print R,$2-{params.padding},$1+{params.padding} }}' | awk -v prev1=0 -v prev2=0 -v OFS='\\t' '$2>0&&$3>0&&$3>$2&&$3>prev2&&prev2!=0 {{print $1,prev1,prev2}} {{prev1=$2;prev2=$3}}'  > $TMPDIR/positions.bed

        cd $TMPDIR
        odgi extract -t {threads} -i {params.graph} -b $TMPDIR/positions.bed -E -d 10000 -s

        for g in *.og
        do
          odgi similarity -t {threads} -i $g | sed 's/:[0-9]\+-[0-9]\+//g' | awk -v S=${{g}} 'NR>1 && $2>$1 {{print S,$1,$2,$6}}'
        done | pigz -p {threads} >> {params._output}
        '''

checkpoint bedtools_makewindows:
    input:
        fai = 'graphs/{chromosome}.fa.gz.fai',
    output:
        directory('graphs/{chromosome}_windows_{size,\d+}')
    params:
        reference_path = 'ARS_UCD1.2#0#{chromosome}',
        split_every = 1000
    localrule: True
    shell:
        '''
        mkdir -p {output}
        bedtools makewindows -w {wildcards.size} -g <(echo -e "{params.reference_path}\\t"$(awk '$1=="{params.reference_path}" {{print $2}}' {input.fai})) |\
        split --lines {params.split_every} -a 4 --additional-suffix .bed -d - {output}/
        '''

rule odgi_extract_windows:
    input:
        og = rules.pggb_construct.output['og'],
        bed = 'graphs/{chromosome}_windows_{size}/{window}.bed'
    output:
        'graphs/{chromosome}_windows_{size}/{window}.jaccard.gz'
    params:
        og = lambda wildcards, input: Path(input['og']).resolve(),
        bed = lambda wildcards, input: Path(input['bed']).resolve(),
        _output = lambda wildcards, output: Path(output[0]).resolve()
    threads: 6
    resources:
        mem_mb = 6000,
        walltime = '4h'
    shell:
        '''
        cd $TMPDIR
        odgi extract -i {params.og} -b {params.bed} -t {threads} -s
        for g in *.og
        do
          odgi similarity -t {threads} -i $g | sed 's/:[0-9]\+-[0-9]\+//g' | awk -v S=${{g}} 'NR>1 && $2>$1 {{print S,$1,$2,$6}}'
        done | pigz -p {threads} > {params._output}
        '''

def aggregate_windows(wildcards):
    checkpoint_output = checkpoints.bedtools_makewindows.get(**wildcards).output[0]
    return expand('graphs/{{chromosome}}_windows_{{size}}/{window}.jaccard.gz', window=glob_wildcards(PurePath(checkpoint_output).joinpath('{window}.bed')).window)

rule gather:
    input:
        aggregate_windows
    output:
        'graphs/{chromosome}.windows.{size}.jaccard.gz'
    localrule: True
    shell:
        '''
        echo "path sample1 sample2 jaccard" | pigz > {output}
        cat {input} >> {output}
        '''

rule summarise_jaccard:
    input:
        rules.gather.output
    output:
        'graphs/{chromosome}.summary.{size}.jaccard.csv'
    params:
        min_paths = 153
    threads: 4
    resources:
        mem_mb = 5000,
        walltime = '10m'
    run:
        import polars as pl
        df = (pl.read_csv(input[0],separator=' ',dtypes=[pl.Utf8,pl.Utf8,pl.Utf8,pl.Float32])
            .with_columns([pl.lit(wildcards.chromosome).alias('Chromosome'),pl.col("path").str.extract(r":(\d+)-").cast(pl.UInt32, strict=False).alias('Start')])
            .drop('path')
         )
        sample_classes = {"white":[f'{S}#{wildcards.chromosome}' for S in ('SIM_2','HER','SIM_3','SIM_1')],
                     "colored":[f'{S}#{wildcards.chromosome}' for S in ["BSW_1","BSW_2","BSW_3","EVO","RGV",
                                                              "BSW_4","BSW_5","BSW_6","BSW_7","OBV_1",
                                                              "OBV_2","OBV_3","OBV_4","BV_1","BV_2",
                                                              "BV_3","BV_4","BV_5","BV_6","HIG"]]}
        invert_samples = {v:K for K,samples in sample_classes.items() for v in samples}

        df = df.with_columns([(pl.col('sample1').replace(invert_samples)==pl.col('sample2').replace(invert_samples)).alias('grouping')]).drop(['sample1','sample2'])
        (df.group_by(['Chromosome','Start','grouping'])
                     .agg(pl.col('jaccard').mean())
                     .join(df.group_by(['Chromosome','Start']).agg(pl.count()),on=['Chromosome','Start'])
                     .filter(pl.col('count')>=params.min_paths)
                     .pivot(index=['Chromosome','Start','count'],columns='grouping',values='jaccard')
                     .with_columns([(pl.col('true')/pl.col('false')).alias('Jaccard similarity ratio')])
                     .drop(['true','false'])
                     .write_csv(output[0]))

rule odgi_procbed:
    input:
        graph = lambda wildcards: expand(rules.odgi_extract.output,pangenome=wildcards.region,renamed='raw'),
        bed = lambda wildcards: config['gggenes'][wildcards.region]['repeats']
    output:
        repeats = 'gggenes/{region}.bed',
        og = 'gggenes/{region}.og'
    params:
        chromosome = lambda wildcards: config['gggenes'][wildcards.region]['region']['ARS'].split(':')[0],
        path = lambda wildcards: PurePath(config['gggenes'][wildcards.region]['repeats']).with_suffix('').with_suffix('').name
    localrule: True
    shell:
        '''
        awk -v OFS='\\t' '$1=={params.chromosome} {{print "{params.path}#{params.chromosome}",$2,$3,$4"_"NR,$6,$5,$7,$8}}' {input.bed} | grep -vE "(Simple_repeat|Low_complexity)" | sort -k 1,1V -k2,2n | tee {output.repeats} |
        odgi procbed -i {input.graph} -b /dev/stdin |\
        odgi inject -i {input.graph} -b /dev/stdin -o {output.og}
        '''
        
rule odgi_untangle:
    input:
        bed = rules.odgi_procbed.output['repeats'],
        og = rules.odgi_procbed.output['og']
    output:
        'gggenes/{region}.tangle',
        'gggenes/{region}.repeats.tsv'
    params:
        jaccard = config['gggenes'].get('jaccard',0.5)
    localrule: True
    shell:
        '''
        odgi untangle -j {params.jaccard} -R <(odgi paths -i {input.og} -L | grep -v "#") -i {input.og} -g | tee {output[0]} |\
        grep -E $(odgi paths -i {input.og} -L | awk 'BEGIN {{printf "^mol|"}} /#/ {{printf $1"|"}}' | sed 's/|$//') |\
        awk -v OFS='\\t' 'BEGIN {{b["gene"]="repeat_class"; c["gene"]="repeat_family"}} NR==FNR {{a[$4]=$5; b[$4]=$7; c[$4]=$8; next}} {{if (a[$2]=="-") {{$5=0}}; print $1,$2,b[$2],c[$2],$3,$4,$5 }}' {input.bed} - |\
        sed -E 's/_[0-9]+//' |\
        sed -E 's/#[0-9]+:[0-9]+-[0-9]+//' > {output[1]}
        '''
