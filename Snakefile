import os, glob, collections

# ---------
# functions
# ---------

def count_complete_replicates_per_sample(bam_dir, sample, nrep):
    # Sometimes all files are available only for Rep1 but not for Rep2. Here, the replicates are iteratively checked for completeness.
    # If for example Rep2 fails the check but Rep1 passed, nrep will be adjusted to 1 and the sample will be returned so that enhancers can be called for the complete Rep1.
    n_complete_reps = 0
    for i, rep in [(i,f'Rep{i}') for i in range(1,nrep+1)]:
        for feature in ('ATAC-seq', 'H3K27ac', 'H3K4me1', 'H3K4me3'):
            bam = glob.glob('%s/%s_%s_%s*.rmdup.bam' %(bam_dir, feature, sample, rep))
            if not len(bam) == 1:
                print('%s of sample %s will be skipped. ' %(rep,sample))
                print('There were %s instead of 1 unique files found with the pattern %s/%s_%s_%s*.rmdup.bam\n' %(len(bam), bam_dir, feature, sample, rep))
                n_complete_reps -= 1
                break
        n_complete_reps += 1
    return sample, n_complete_reps

def get_samples(samples, bam_dir):
    # function to return samples in the bam_dir and distinguish those that are available in 1 replicate from those that are available in multiple replicates
    # for every sample, check if all the bam files are there and skip the sample (or adjust the number of complete replicates) if files are missing.
    print('_'*100)
    if samples == 'all':
        atac_only_message = '\nATAC_ONLY mode: Only calling ATAC-seq peaks.' if atac_only else ''
        print(atac_only_message + '\nFetching sample names from %s:' %bam_dir)
        atac_files = glob.glob('%s/ATAC-seq*%s*Rep*.rmdup.bam' %(bam_dir,build)) # build: make sure only samples matching the passed build are taken in case there are data from multiple species in the bam_dir
    else:
        atac_files = sum([glob.glob('%s/ATAC-seq_%s_Rep*.rmdup.bam' %(bam_dir,sample)) for sample in samples.split(',')], [])
    print('\nBam files found: (ATAC-seq)\n%s' %'\n'.join(['%s' %sample for sample in atac_files]))
    atac_samples = ['_'.join(os.path.basename(x).split('_')[:5]) for x in atac_files]
    print('\nUnique samples found: (ATAC-seq)\n%s\n' %'\n'.join(['%s' %sample.replace('ATAC-seq_','') for sample in set(atac_samples)]))
    atac_samples_nrep_counter = collections.Counter(atac_samples)
    if atac_only:
        print('Valid samples:\n' + '\n'.join(sum([['%s Rep%s' %(sample.replace('ATAC-seq_',''), rep) for rep in range(1,nreps+1)] for sample, nreps in atac_samples_nrep_counter.items() if not nreps == 0], [])) + '\n' + '_'*100 + '\n')
        return {sample.replace('ATAC-seq_',''): count for sample, count in atac_samples_nrep_counter.items() if count > 0}
    samples_nrep_dict = {sample: n_complete_reps for sample, n_complete_reps in [count_complete_replicates_per_sample(bam_dir, sample.replace('ATAC-seq_', ''), count) for sample, count in atac_samples_nrep_counter.items()] if n_complete_reps > 0}
    print('Complete samples (ATAC-seq and Histone Modifications):\n' + '\n'.join(sum([['%s Rep%s' %(sample, rep) for rep in range(1,n_complete_reps+1)] for sample, n_complete_reps in samples_nrep_dict.items()], [])) + '\n' + '_'*100 + '\n')
    return samples_nrep_dict

def get_targets(sample, nrep, merge_reps=None, transcripts_file=None, atac_only=False):
    if atac_only:
        rep_peaks = ['%s/ATAC-seq_%s_Rep%s.%s_peaks.bed' %(atac_peaks_dir, sample, i, peak_caller) for i in range(1,nrep+1)]
        combined_peaks = []
        if nrep > 1:
            if peak_caller == 'macs2':
                combined_peaks = ['%s/ATAC-seq_%s.macs2_peaks.idr.bed' %(atac_peaks_dir, sample)]
            elif peak_caller == 'genrich':
                combined_peaks = ['%s/ATAC-seq_%s.genrich_peaks.bed' %(atac_peaks_dir, sample)]
        return rep_peaks + combined_peaks
    else:
        if nrep == 1:
            enhancers_rep = ['%s/%s_Rep1.enhancers.filtered.bed' %(outdir, sample)]
            targets = enhancers_rep # + [re.sub('.bed', '.stats.pdf', enhancers)]
        else:
            if merge_reps == 'intersect_predictions':
                enhancers_intersect = ['%s/%s.enhancers.intersect_predictions.bed' %(outdir, sample)]
                enhancers_rep = ['%s/%s_Rep%i.enhancers.filtered.bed' %(outdir, sample, i) for i in range(1,(nrep+1))]
                targets = sum([enhancers_intersect, enhancers_rep], [])
                # targets += [re.sub('.bed', '.stats.pdf', enhancers_intersect)]
            elif merge_reps == 'merge_bams':
                enhancers = ['%s/%s_merged.enhancers.filtered.bed' %(outdir, sample)]
                merged_bam_files = ['%s/%s_%s_merged.cpm.bw' %(bw_dir, feature, sample) for feature in ('H3K27ac','H3K4me1','H3K4me3')]
                targets = sum([enhancers, merged_bam_files], [])
                # targets += [re.sub('.bed', '.stats.pdf', enhancers)]
    return targets

def get_replicate_bam(feature, sample, rep, index=False):
    # rep has to be passed again separately in case of the merge_bams rule where the wildcards' rep is "_merged", but we want to find the Rep1 and Rep2 input files
    if rep == 'merged':
        bam = '%s/%s_%s_merged.rmdup.bam' %(bam_dir, feature, sample)
        return bam if not index else bam + '.csi'
    else:
        bamfiles = glob.glob('%s/%s_%s_%s*.rmdup.bam' %(bam_dir, feature, sample, rep))
        if not len(bamfiles) == 1:
            raise FileNotFoundError('There were %s instead of 1 unique files found with the pattern %s/%s_%s_%s*.rmdup.bam' %(len(bamfiles), bam_dir, feature, sample, rep))
        return bamfiles[0] if not index else bamfiles[0] + '.csi'

def get_input_for_enhancer_filtering(wc):
    if (merge_reps == 'intersect_predictions') or (samples_nrep_dict[wc.sample] == 1):
        atac_peaks = '%s/ATAC-seq_%s_%s.%s_peaks_%s.bed' %(atac_peaks_dir, wc.sample, wc.rep, peak_caller, pval)
    elif (merge_reps == 'merge_bams') and (samples_nrep_dict[wc.sample] > 1):
        if peak_caller == 'macs2':
            atac_peaks = '%s/ATAC-seq_%s.macs2_peaks_%s.idr.bed' %(atac_peaks_dir, wc.sample, pval)
        elif peak_caller == 'genrich': # genrich can handle multiple replicates, making IDR obsolete
            atac_peaks = '%s/ATAC-seq_%s.genrich_peaks_%s.bed' %(atac_peaks_dir, wc.sample, pval)
    crup_predictions = '%s/%s_%s.enhancers.crup.bed' %(crup_dir, wc.sample, wc.rep)
    promoters = '%s/%s_%s.promoters.bed' %(promoter_dir, wc.sample, wc.rep)
    transcripts = transcripts_file
    return [atac_peaks, crup_predictions, promoters, transcripts]

# convert config dict entries to variables, e.g. a = d['a']. (bam_dir, samples, build, transcripts, merge_reps, peak_caller, crup_cutoff, outdir, atac_only)
for k,v in config.items():
    exec(k + '=v')

# convert passed directory / file paths to absolute paths
bam_dir = os.path.abspath(config['bam_dir'])
bw_dir = bam_dir.replace('bam','bigwig')
outdir = os.path.abspath(config['outdir'])
transcripts_file = os.path.abspath(transcripts) if not transcripts == 'UCSC' else '%s/transcripts_%s_UCSC.bed' %(outdir, build)

# variables
crup_dir = '%s/crup' %outdir
atac_peaks_dir = '%s/%s' %(outdir, peak_caller)
promoter_dir = '%s/ehmm' %outdir
fai = '/project/MDL_ChIPseq/data/genome/fasta/%s.fa.fai' %build
sizes = '/project/MDL_ChIPseq/data/genome/assembly/%s.sizes' %build
genome_size = sum([int(x.strip().split('\t')[1]) for x in open(sizes).readlines()])

# identify samples: either a comma-separated string of sample names or "all" for all samples in the bam_dir that come in >1 replicates
samples_nrep_dict = get_samples(samples, bam_dir)

# define targets
targets = sum([get_targets(sample, nrep, merge_reps, transcripts_file, atac_only) for sample, nrep in samples_nrep_dict.items()], [])

# -----
# rules
# -----

rule all:
    input:
        targets
    params:
        peak_caller = peak_caller
    shell:
        'prun mdl cite crup ehmm {params.peak_caller}'

rule get_transcripts:
    output:
        '%s/transcripts_%s_UCSC.bed' %(outdir, build)
    params:
        build = build,
        outdir = outdir
    shell:
        'prun mdl get_transcripts.R {params}'
    
rule intersect_replicate_predictions:
    # Take the union of replicate-specific filtered enhancers
    input:
        lambda wc: ['%s/%s_Rep%s.enhancers.filtered.bed' %(wc.outdir, wc.sample, rep) for rep in range(1,(samples_nrep_dict[wc.sample]+1))]
    output:
        '{outdir}/{sample,((?!merged).)*}.enhancers.intersect_predictions.bed'
    params:
        nrep = lambda wc: samples_nrep_dict[wc.sample]
    shell:
        'prun mdl intersect_replicate_enhancers.R {output} {params.nrep} {input}'

rule plot_stats:
    input:
        '{outdir}/{sample}.bed'
    output:
        '{outdir}/{sample}.stats.pdf'
    params:
        nrep = lambda wc: samples_nrep_dict[wc.sample],
        transcripts = transcripts,
        merge_reps = merge_reps,
        peak_caller = peak_caller,
        build = build
    shell:
        'plot_stats.R {output} {params}'
    
rule filter_enhancers:
    input:
        get_input_for_enhancer_filtering
    output:
        '{outdir}/{sample}_{rep}.enhancers.filtered.bed'
    params:
        build = build
    shell:
        'prun mdl enhancer_filter.R {output} {input} {params}'

rule genome_size_to_bed:
    input:
        ancient(sizes)
    output:
        '{dir}/{build}.sizes.bed'
    shell:
        'cat {input} | awk -F"\t" -v OFS="\t" "{{\$2=int(\$2/100)*100}} \$1 !~ /Un|random|chrM/ {{print \$1,100,\$2}}" > {output}'
        
rule call_promoters:
    input:
        bam_dir = ancient(bam_dir),
        sizes = ancient(sizes),
        regions = '{promoter_dir}/%s.sizes.bed' %build,
        atac = lambda wc: get_replicate_bam('ATAC-seq', wc.sample, wc.rep),
        atac_idx = lambda wc: get_replicate_bam('ATAC-seq', wc.sample, wc.rep, index=True),
        h3k27ac = lambda wc: get_replicate_bam('H3K27ac', wc.sample, wc.rep),
        h3k27ac_idx = lambda wc: get_replicate_bam('H3K27ac', wc.sample, wc.rep, index=True),
        h3k4me1 = lambda wc: get_replicate_bam('H3K4me1', wc.sample, wc.rep),
        h3k4me1_idx = lambda wc: get_replicate_bam('H3K4me1', wc.sample, wc.rep, index=True),
        h3k4me3 = lambda wc: get_replicate_bam('H3K4me3', wc.sample, wc.rep),
        h3k4me3_idx = lambda wc: get_replicate_bam('H3K4me3', wc.sample, wc.rep, index=True)
    output:
        bed = '{promoter_dir}/{sample}_{rep,[a-zA-Z0-9]*}.promoters.bed',
        outdir = directory('{promoter_dir}/{sample}_{rep,[a-zA-Z0-9]*}')
    threads:
        min(workflow.cores, 20)
    shell:
        '''
        mkdir -p {output.outdir}
        prun mdl ehmm applyModel --regions {input.regions} --mark ATAC-seq:{input.atac} --mark H3K27AC:{input.h3k27ac} --mark H3K4ME1:{input.h3k4me1} --mark H3K4ME3:{input.h3k4me3} --genomeSize {input.sizes} --outdir {output.outdir} -n {threads}
        ln -sr {output.outdir}/promoterRegions.bed {output.bed}
        '''

rule idr:
    input:
        lambda wc: expand('%s/ATAC-seq_%s_Rep{i}.%s_peaks.narrowPeak' %(wc.atac_peaks_dir, wc.sample, wc.peak_caller), i=range(1,(samples_nrep_dict[wc.sample]+1)))
    output:
        '{atac_peaks_dir}/ATAC-seq_{sample}.{peak_caller}_peaks.idr'
    shell:
        'idr --idr-threshold 0.05 --input-file-type narrowPeak -s {input} -o {output}'

rule idr_to_bed:
    input:
        '{atac_peaks_dir}/ATAC-seq_{sample}.{peak_caller}_peaks.idr'
    output:
        '{atac_peaks_dir}/ATAC-seq_{sample}.{peak_caller}_peaks.idr.bed'
    shell:
        'sort -k1,1 -k2,2n {input} | awk -F "\\t" -v OFS="\\t" "{{print \$1,\$2,\$3,\\"idr_\\"NR,\$5,\$6}}" > {output}'

rule name_sort:
    input:
        '{sample}.rmdup.bam'
    output:
        '{sample}.rmdup.nsort.bam'
    threads:
        min(workflow.cores, 30)
    shell:
        'samtools sort -n -@ {threads} -o {output} {input}'
        
rule call_single_replicate_atac_peaks_genrich:
    input:
        lambda wc: get_replicate_bam('ATAC-seq', wc.sample, wc.rep).replace('.rmdup.bam','.rmdup.nsort.bam')
    output:
        '{atac_peaks_dir}/ATAC-seq_{sample}_{rep,((?!merged).)*}.genrich_peaks_%s.narrowPeak' %pval
    params:
        pval = pval
    shell:
        'prun mdl Genrich -j -y -p {params.pval} -t {input} -o {output}' # -j: ATAC-seq mode, -y: keep unpaired reads (some old ATAC are single-end), permissive p-value (default 0.01)

rule call_multi_replicate_atac_peaks_genrich:
    # genrich can handle multiple replicates, making IDR obsolete.
    input:
        lambda wc: [get_replicate_bam('ATAC-seq', wc.sample, 'Rep%i' %i).replace('.rmdup.bam','.rmdup.nsort.bam') for i in range(1,(samples_nrep_dict[wc.sample]+1))]
    output:
        '{atac_peaks_dir}/ATAC-seq_{sample}.genrich_peaks_%s.narrowPeak' %pval
    params:
        pval = pval
    shell:
        'prun mdl Genrich -j -y -p {params.pval} -t "{input}" -o {output}' # -j: ATAC-seq mode, -y: keep unpaired reads (some old ATAC are single-end), permissive p-value (default 0.01)
        
rule call_atac_peaks_macs2:
    input:
        sizes = ancient(sizes),
        bam = lambda wc: get_replicate_bam('ATAC-seq', wc.sample, wc.rep)
    output:
        '{atac_peaks_dir}/ATAC-seq_{sample}_{rep}.macs2_peaks_%s.narrowPeak' %pval,
    params:
        genome_size = genome_size,
        name = 'ATAC-seq_{sample}_{rep}.macs2',
        pval = pval
    shell:
        'macs2 callpeak -t {input.bam} --outdir {wildcards.atac_peaks_dir} -n {params.name} -g {params.genome_size} -q {params.pval}'

rule narrowPeak_to_bed:
    input:
        '{sample}.{peak_caller}_peaks_%s.narrowPeak' %pval
    output:
        '{sample}.{peak_caller}_peaks_%s.bed' %pval
    shell:
        'awk -F "\\t" -v OFS="\\t" "{{print \$1,\$2,\$3,\$4,\$5,\$6}}" {input} > {output}'
        
rule bedgraph_to_bed:
    # crup's predictions may contain negative coordinates (removed with awk) and such that extend chromosome sizes (removed with bedtools slop)
    input:
        bedgraph = '{crup_dir}/{sample}_{rep}/singleEnh.bedGraph',
        sizes = ancient(sizes)
    output:
        '{crup_dir}/{sample}_{rep}.enhancers.crup.bed'
    shell:
        'cat {input.bedgraph} | awk -F"\t" -v OFS="\t" "{{if (\$2<0) \$2=0}} {{print \$0}}" | bedtools slop -i /dev/stdin -b 0 -g {input.sizes} > {output}'
    
rule call_crup_enhancers:
    input:
        data = '{crup_dir}/{sample}_{rep}/data_matrix.rds',
        sizes = ancient(sizes)
    output:
        '{crup_dir}/{sample}_{rep}/singleEnh.bedGraph'
    log:
        '{crup_dir}/{sample}_{rep}/crup.log'
    threads:
        min(workflow.cores, 10)
    params:
        crup_cutoff = crup_cutoff
    shell:
        '''
        mkdir -p {wildcards.crup_dir}/{wildcards.sample}_{wildcards.rep}
        prun mdl CRUP.R -P -x {threads} -m {input.data} -u {params.crup_cutoff} -o {wildcards.crup_dir}.{wildcards.sample}_{wildcards.rep} > {log}
        '''

rule normalize_crup_data:
    input:
        '{data_dir}/data_table.txt'
    output:
        '{data_dir}/data_matrix.rds'
    log:
        '{data_dir}/data_matrix.log'
    params:
        fai = fai
    threads:
        min(workflow.cores, 10)
    shell:
        '''
        sequencing_type=$(prun mdl get_sequencing_type.py {input} {threads})
        prun mdl CRUP.R -N -I -x {threads} -g {params.fai} -s $sequencing_type -i {input} > {log}
        mv {wildcards.data_dir}/data_table.data_matrix.rds {output}
        '''

rule create_crup_data_table:
    input:
        h3k27ac = lambda wc: get_replicate_bam('H3K27ac', wc.sample, wc.rep),
        h3k4me1 = lambda wc: get_replicate_bam('H3K4me1', wc.sample, wc.rep),
        h3k4me3 = lambda wc: get_replicate_bam('H3K4me3', wc.sample, wc.rep)        
    output:
        tbl = '{crup_dir}/{sample}_{rep,[a-zA-Z0-9]*}/data_table.txt',
        outdir = directory('{crup_dir}/{sample}_{rep,[a-zA-Z0-9]*}')
    shell:
        '''
        mkdir -p {output.outdir}
        echo -e "feature\tbam_file\nH3K27ac\t{input.h3k27ac}\nH3K4me1\t{input.h3k4me1}\nH3K4me3\t{input.h3k4me3}" > {output.tbl}
        '''

rule merge_bams:
    input:
        lambda wc: [get_replicate_bam(wc.feature, wc.sample, 'Rep%i' %i) for i in range(1,(samples_nrep_dict[wc.sample]+1))]
    output:
        '{bam_dir}/{feature,((?![/_]).)*}_{sample}_merged.rmdup.bam'
    params:
        bw_dir = bw_dir
    threads:
        min(workflow.cores, 20)
    shell:
        'samtools merge -@ {threads} {output} {input}'

rule bw_atac:
    input:
        bam = '%s/{sample}_merged.rmdup.bam' %bam_dir,
        csi = '%s/{sample}_merged.rmdup.bam.csi' %bam_dir
    output:
        '%s/{sample,ATAC.*}_merged.cpm.bw' %bw_dir
    log:
        '%s/{sample}_merged.bamCoverage.log' %bam_dir
    threads:
        min(workflow.cores, 20)
    shell:
        'bamCoverage --binSize 10 --normalizeUsing CPM --minMappingQuality 30 -p {threads} -b {input.bam} -o {output} > {log}'

rule bw:
    input:
        bam = '%s/{sample}_merged.rmdup.bam' %bam_dir,
        csi = '%s/{sample}_merged.rmdup.bam.csi' %bam_dir
    output:
        '%s/{sample,(?!ATAC).*}_merged.cpm.bw' %bw_dir
    log:
        '%s/{sample}_merged.bamCoverage.log' %bam_dir
    threads:
        min(workflow.cores, 20)
    shell:
        'bamCoverage --binSize 10 --normalizeUsing CPM --centerReads --minMappingQuality 30 -p {threads} -b {input.bam} -o {output} > {log}'

rule index:
    input:
        '{bam_dir}/{sample}_merged.rmdup.bam'
    output:
        '{bam_dir}/{sample}_merged.rmdup.bam.csi'
    threads:
        min(workflow.cores, 20)
    shell:
        'samtools index -c -@ {threads} {input}'
