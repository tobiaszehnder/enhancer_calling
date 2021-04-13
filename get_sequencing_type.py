#! /usr/bin/env python

def main():
    import sys, os
    _, data_table, nthreads = sys.argv
    
    # fetch the canonical path of the specified filename, eliminating any symbolic links encountered in the path, and determine if it contains "paired-end" or "single-end"
    with open(data_table, 'r') as f:
        lines = f.readlines()
    files = [x for x in sum([x.strip().split('\t') for x in lines], []) if os.path.isfile(x)]
    flag_sums = [int(os.popen('samtools view -@ %s -c -f 1 %s' %(nthreads,fn)).read().strip()) for fn in files]
    if all([flag_sum == 0 for flag_sum in flag_sums]):
        print('single')
        return 'single'
    elif all([flag_sum > 0 for flag_sum in flag_sums]):
        print('paired')
        return 'paired'
    else:
        raise ValueError('Histone modification bam files in %s must be of equal sequencing type (either single-end or paired-end).' %data_table)

if __name__ == '__main__':
    main()
