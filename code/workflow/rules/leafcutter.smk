rule IndexSampleID_ColumnID:
    '''
    output all sample ids (rail_id) along with corresponding column number in 
    the `resources/snaptron-yil/run/full/All_perind.constcounts.gz` file.
    '''
    message: 'Gather Sample Rail IDs'
    input: 
        'resources/snaptron-yil/run/full/All_perind.constcounts.gz'
    output: 
        'resources/All_perind_rail_id_column_index.txt'
    log: 
        'logs/GatherSampleRailIDs.log'
    threads: 1
    resources: cpu = 1, mem_mb = 15000, time = 2100
    shell: 
        '''
        printf "%s\n" $(zcat {input} | head -1 | sed -r 's/chrom //') | awk '{{print $1, NR+1}}' 1> {output} 2> {log}
        
        '''

rule MakeSampleLookupTable:
    message: 'Make lookup tablef or sample ids(rail_id)'
    input: 
        gtex = 'resources/snaptron-yil/scripts/data/samples_gtex.tsv',
        tcga = 'resources/snaptron-yil/scripts/data/samples_tcga.tsv',
        srav = 'resources/snaptron-yil/scripts/data/samples_srav2.tsv',
        all_perind = rules.IndexSampleID_ColumnID.output
    output: 
        'resources/All_perind_rail_id_tissue_type_lookup.tsv'
    threads: 1
    resources: cpu = 1, mem_mb = 15000, time = 2100
    script: '../scripts/MakeSampleLookupTable.R'



Perind_RailID_Tissue_Lookup_df = pd.read_csv(
    config['RailID_Tissue_Lookup'], sep='\t').set_index(
        keys=['rail_id', 'tissue', 'ds'], drop=False)

def getColumnIndexFromCountsFile(wildcards):
    Columns = Perind_RailID_Tissue_Lookup_df.query(
        'ds == @wildcards.data_source & tissue == @wildcards.tissue' 
    )['counts_col_id'].reset_index(drop=True).unique()
    Columns = ['1'] + [str(i) for i in Columns] # add in the first column in *count.gz file, which is intron coordinates
    Columns = ','.join(Columns)
    return Columns

rule GetLeafcutterCounts:
    message: 'Subset leafcutter counts output by sample'
    input: 'resources/snaptron-yil/run/full/All_perind.constcounts.gz'
    # input: 'test.counts.gz'
    output: 'results/SubsetCounts/{data_source}/{tissue}/perind.counts.gz'
    log: 'logs/SubsetCounts/{data_source}/{tissue}.log'
    params: 
        Fields = getColumnIndexFromCountsFile
    threads: 1
    resources: cpu = 1, mem_mb = 20000, time = 2100
    shell: 
        '''
        zcat {input} | cut -d " " -f {params.Fields} | bgzip -c 1> {output} 2> {log}
        '''

rule GetLeafcutterNumerators:
    message: 'Get the numerators and get input ready for leafcutterMD analyses.'
    input: rules.GetLeafcutterCounts.output 
    output: 'results/SubsetCounts/{data_source}/{tissue}/perind_num.counts.gz'
    log: 'logs/SubsetCounts/{data_source}/{tissue}_nums.log'
    threads: 1
    resources: cpu = 1, mem_mb = 25000, time = 2100
    shell: 
        '''
        cat <(zcat {input} | head -1 | sed -r 's/chrom //')   `# remove "chrom " from header` \
            <(zcat {input} | awk 'NR > 1' | sed -r 's/\/[0-9]+//g') |   `# extract only numerators` \
            bgzip -c 1> {output} 2> {log}
        '''


rule FilterLeafcutterNumerators:
    message: 'Filter out clusters that have <=5 total reads'
    input: rules.GetLeafcutterNumerators.output
    output: 'results/SubsetCounts/{data_source}/{tissue}/perind_num_filtered.counts.gz'
    threads: 1
    resources: cpu = 1 , mem_mb = 25000, time = 2100
    shell: 
        '''
        unzipped=$(echo {output} | sed 's/counts\.gz/counts/')
        Rscript workflow/scripts/FilterLeafCutterCountsByCluster.R {input} $unzipped
        cat <(head -1 $unzipped | sed -r 's/chrom //') <(awk 'NR > 1' $unzipped) | bgzip -c > {output}
        rm $unzipped
        '''


rule LeafcutterMD:
    message: 'Use leafcutterMD to call outlier splicing'
    input: rules.FilterLeafcutterNumerators.output
    output: 'results/leafcutterMD/{data_source}/{tissue}/outlier_pVals.txt'
    log: 'logs/LeafcutterMD/{data_source}/{tissue}.log'
    params: 
        output_prefix = 'results/leafcutterMD/{data_source}/{tissue}/outlier'
    threads: 4 
    resources: cpu = 4, mem_mb = 30000, time = 2100
    shell: 
        '''
        workflow/scripts/leafcutter/scripts/leafcutterMD.R \
            -p {threads} \
            --output_prefix {params.output_prefix} \
            {input} &> {log}
        '''