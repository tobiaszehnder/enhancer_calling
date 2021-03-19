# import modules
import numpy as np, pandas as pd, re

# convert dict entries to variables, e.g. a = d['a']
for k,v in config.items():
    exec(k + '=v')

# define variables
crup_dir = '%s/crup' %outdir
atac_peaks_dir = '%s/%s' %(outdir, peak_caller)
promoter_dir = '%s/ehmm' %outdir
fai = '/project/MDL_ChIPseq/data/genome/fasta/%s.fa.fai' %build
sizes = '/project/MDL_ChIPseq/data/genome/assembly/%s.sizes' %build

# define targets
def get_targets(merge_reps):
    if merge_reps:
        final_enhancers = '%s/%s_merged.enhancers.bed' %(outdir, sample)
        enhancers_rep = '%s/%s_merged.enhancers.filtered.bed' %(outdir, sample)
    else:
        final_enhancers = '%s/%s.enhancers.bed' %(outdir, sample)
        enhancers_rep = ['%s/%s_Rep%s.enhancers.filtered.bed' %(outdir, sample, rep) for rep in np.arange(nrep)+1]
    stats = re.sub('.bed', '.stats.pdf', final_enhancers)
    return final_enhancers, enhancers_rep, stats

rule all:
    input: get_targets(merge_reps)
    params:
        peak_caller = peak_caller
    shell:
        'prun mdl cite crup ehmm {params.peak_caller}'

def get_filtered_enhancers_per_replicate(wc, nrep, merge_reps):
    print(merge_reps)
    if merge_reps == True:
        return '%s/%s_merged.enhancers.filtered.bed' %(wc['outdir'], wc['sample'])
    else:
        return ['%s/%s_Rep%s.enhancers.filtered.bed' %(wc['outdir'], wc['sample'], rep) for rep in np.arange(nrep)+1]
    
rule combine_replicate_predictions:
    # Take the union of replicate-specific filtered enhancers
    input:
        lambda wc: get_filtered_enhancers_per_replicate(wc, nrep, merge_reps)
    output:
        '{outdir}/{sample,((?!merged}.)*}.enhancers.bed'
    params:
        nrep = nrep
    shell:
        'prun mdl enhancer_combine_replicates.R {output} {params.nrep} {input}'

rule link_bam_merged_enhancers:
    # In case the bams were merged, just link the {enhancers}.filtered.bed to {enhancers}.bed
    input:
        '{outdir}/{sample}_merged.enhancers.filtered.bed'
    output:
        '{outdir}/{sample}_merged.enhancers.bed'
    shell:
        'ln -sf {input} {output}'
        
rule plot_stats:
    input:
        '{outdir}/{sample}.enhancers.bed'
    output:
        '{outdir}/{sample}.enhancers.stats.pdf'
    params:
        nrep = nrep,
        transcripts = transcripts,
        merge_reps = merge_reps,
        peak_caller = peak_caller
    shell: 'enhancer_stats.R {output} {params}'

rule filter_enhancers:
    input:
        atac_peaks = '%s/{sample}_{rep}.atac_peaks.bed' %atac_peaks_dir,
        crup_predictions = '%s/{sample}_{rep}.enhancers.crup.bed' %crup_dir,
        promoters = '%s/{sample}_{rep}.promoters.bed' %promoter_dir
    output:
        '{outdir}/{sample}_{rep}.enhancers.filtered.bed'
    params:
        transcripts = transcripts,
        build = build
    shell:
        'prun mdl enhancer_filter.R {output} {input} {params}'

rule call_promoters:
    input:
        bam_dir = bam_dir,
        sizes = sizes
    output:
        '{promoter_dir}/{sample}.promoters.bed'
    shell:
        '''
        mkdir -p {wildcards.promoter_dir}/{wildcards.sample}
        prun mdl call_promoters {wildcards.sample} {input.bam_dir} {input.sizes} {wildcards.promoter_dir}/{wildcards.sample}
        ln -srv {wildcards.promoter_dir}/{wildcards.sample}/promoterRegions.bed {output}
        '''

rule call_atac_peaks:
    input:
        sizes
    output:
        '{atac_peaks_dir}/{sample}_{rep}.atac_peaks.bed'
    params:
        peak_caller = peak_caller,
        pval = pval,
        bam_dir = bam_dir
    shell:
        '''
        mkdir -p {wildcards.atac_peaks_dir}/{wildcards.sample}_{wildcards.rep}
        prun mdl call_atac_peaks {params.peak_caller} {params.pval} {params.bam_dir} {wildcards.sample}_{wildcards.rep} {input} {output}
        '''

rule call_crup_enhancers:
    input:
        '{crup_dir}/{sample}_{rep}/data_matrix.rds'
    output:
        '{crup_dir}/{sample}_{rep}.enhancers.crup.bed'
    log:
        '{crup_dir}/{sample}_{rep}.enhancers.crup.log'
    threads:
        min(workflow.cores, 20)
    params:
        crup_cutoff = crup_cutoff
    shell:
        '''
        mkdir -p {wildcards.crup_dir}/{wildcards.sample}_{wildcards.rep}
        prun mdl CRUP.R -P -x {threads} -m {input} -u {params.crup_cutoff} -o {sample} > {log}
        ln -sr {wildcards.crup_dir}/{wildcards.sample}_{wildcards.rep}/singleEnh.bedGraph {output}
        '''

rule normalize_crup_data:
    input:
        '{data_dir}/data_table.txt'
    output:
        '{data_dir}/data_matrix.rds'
    log:
        '{data_dir}/data_matrix.log'
    params:
        fai = fai,
        sequencing_type = sequencing_type
    threads:
        min(workflow.cores, 20)
    shell:
        '''
        prun mdl CRUP.R -N -I -x {threads} -g {params.fai} -s {params.sequencing_type} -i {input} > {log}
        mv {wildcards.data_dir}/data_table.data_matrix.rds {output}
        '''

rule create_crup_data_table:
    input:
        bam_dir
    output:
        '{data_dir}/data_table.txt'
    shell:
        '''
        mkdir -p {wildcards.data_dir}
        prun mdl create_crup_datatable $(basename {wildcards.data_dir}) {input} {output}
        '''
