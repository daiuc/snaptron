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


def getColumnIndexFromCountsFile(wildcards):
    Columns = SampleLookupTable.query(
        'ds == @wildcards.data_source & tissue == @wildcards.tissue' 
    )['counts_col_id'].reset_index(drop=True).unique()
    #Columns = ['1'] + [str(i) for i in Columns] # add in the first column in *count.gz file, which is intron coordinates
    #Columns = ','.join(Columns)
    Columns = [str(i) for i in Columns]
    return Columns


rule GetLeafcutterCounts:
    message: 'Subset leafcutter counts output by sample (also by groups)'
    input: 
        'resources/snaptron-yil/run/full/All_perind.constcounts.gz'
    output: 
        d = directory('results/SubsetCounts/{data_source}/{tissue}/counts'),
        flag = touch('results/SubsetCounts/{data_source}/{tissue}/counts/GetLeafcutterCounts.done')
    params: 
        Fields = getColumnIndexFromCountsFile,
        Chunk_Size = 100,
        input = 'resources/snaptron-yil/run/full/All_perind.constcounts.gz',
        output_prefix = 'results/SubsetCounts/{data_source}/{tissue}/counts/group'
    threads: 1
    resources: cpu = 1, mem_mb = 15000, time = 2100
    run:
        import subprocess

        N_chunks = max(1, round(len(params.Fields)/params.Chunk_Size)) # make sure at least 1 chunk
        cut_fields = []
        print(wildcards)
        print('## Splitting ' + str(len(params.Fields)) + ' samples into ' + str(N_chunks) + 
            ' chunks, chunk size = ' + str(params.Chunk_Size) + ' samples.')
        if os.path.exists(str(output.d)):
            print('Path - ' + str(output.d) + ' exists? True')
        else:
            print('Path - ' + str(output.d) + ' exists? False. Create path.')
            os.makedirs(str(output.d), exist_ok=True)

        for r in range(N_chunks):
            if r < N_chunks - 1:
                cut_fields.append([params.Fields[r * params.Chunk_Size + i] for i in range(params.Chunk_Size)])
            else:
                cut_fields.append(params.Fields[r * params.Chunk_Size: len(params.Fields)])
        
        r = 0 # r: run for chunk of samples
        for cf in cut_fields:
            cmd = 'zcat ' + params.input + ' | cut -d " " -f 1,' + ','.join(cf) + \
                ' | bgzip -c 1> ' + params.output_prefix + str(r) + '_perind.counts.gz'
            print('\n## Chunk ' + str(r) + ' : Running command: \n' + cmd)
            try: 
                res = subprocess.run(cmd, shell=True, capture_output=True, check=True)
            except CalledProcessError:
                print("Bash cmd error!")
            if len(res.stderr) > 0:
                print('# !!Failed!! bash command run was unsuccessful! Check code.')
                print(res.stderr)
            else:
                print('# success!')
            
            r = r + 1
        
        print("\n### All Done!")


def getGetLeafcutterNumeratorsInput(wildcards):
    input_files = glob.glob1(
        'results/SubsetCounts/' + wildcards.data_source + '/' + wildcards.tissue + '/counts', 
        'group*_perind.counts.gz')
    input_files = ['results/SubsetCounts/' + wildcards.data_source + '/' + wildcards.tissue + \
        '/counts/' + f for f in input_files]
    return input_files

rule GetLeafcutterNumerators:
    message: 'Get the numerators and get input ready for leafcutterMD analyses.'
    input: 
        counts = getGetLeafcutterNumeratorsInput,
        pull_previous_step = rules.GetLeafcutterCounts.output.flag
    output: 
        d = directory('results/SubsetCounts/{data_source}/{tissue}/nums'),
        flag = touch('results/SubsetCounts/{data_source}/{tissue}/nums/GetLeafcutterNumerators.done')
    log: 'logs/SubsetCounts/{data_source}/{tissue}_nums.log'
    threads: 1
    resources: cpu = 1, mem_mb = 25000, time = 2100
    shell: 
        '''
        mkdir -p {output.d}
        for f in {input.counts}; do
            echo "## Extracting numerators from $f"
            group=$(basename -s _perind.counts.gz $f)
            cat <(zcat $f | head -1 | sed -r 's/chrom //')   `# remove "chrom " from header` \
                <(zcat $f | awk 'NR > 1' | sed -r 's/\/[0-9]+//g') |   `# extract only numerators` \
                bgzip -c 1> {output.d}/${{group}}_perind.nums.gz 2> {log}
        done
        '''


def getFilterLeafcutterNumeratorsInput(wildcards):
    input_files = glob.glob1(
        'results/SubsetCounts/' + wildcards.data_source + '/' + wildcards.tissue + '/nums', 
        'group*_perind.nums.gz')
    input_files = ['results/SubsetCounts/' + wildcards.data_source + '/' + wildcards.tissue + \
        '/nums/' + f for f in input_files]
    return input_files

rule FilterLeafcutterNumerators:
    message: 'Filter out clusters that have <=5 total reads'
    input: 
        nums = getFilterLeafcutterNumeratorsInput,
        pull_previous_step = rules.GetLeafcutterNumerators.output.flag
    output: 
        d = directory('results/SubsetCounts/{data_source}/{tissue}/nums_filtered'),
        flag = touch('results/SubsetCounts/{data_source}/{tissue}/nums_filtered/FilterLeafcutterNumerators.done')
    threads: 1
    resources: cpu = 1 , mem_mb = 25000, time = 2100
    shell: 
        '''
        mkdir -p {output.d}
        for f in {input.nums}; do
            echo -e "\n\n### Filtering ${{f}} ###\n\n"
            group=$(basename -s _perind.nums.gz $f)
            out_f={output.d}/${{group}}_perind.nums.gz

            unzipped=$(echo $out_f | sed 's/nums\.gz/nums/')
            Rscript workflow/scripts/FilterLeafCutterCountsByCluster.R $f $unzipped
            cat <(head -1 $unzipped | sed -r 's/chrom //') \
                <(awk 'NR > 1' $unzipped) | bgzip -c > $out_f
            rm $unzipped
        done
        '''


def getLeafcutterMDInput(wildcards):
    input_files = glob.glob1(
        'results/SubsetCounts/' + wildcards.data_source + '/' + wildcards.tissue + \
        '/nums_filtered', 'group*_perind.nums.gz')
    input_files = ['results/SubsetCounts/' + wildcards.data_source + '/' + 
        wildcards.tissue + '/nums_filtered/' + f for f in input_files]
    return input_files

rule LeafcutterMD:
    message: 'Use leafcutterMD to call outlier splicing'
    input: 
        nums_filtered = getLeafcutterMDInput,
        pull_previous_step = rules.FilterLeafcutterNumerators.output.flag
    output: 
        d = directory('results/leafcutterMD/{data_source}/{tissue}'),
        flag = touch('results/leafcutterMD/{data_source}/{tissue}/LeafcutterMD.done')
    log: 'logs/LeafcutterMD/{data_source}/{tissue}.log'
    threads: 8 
    resources: cpu = 8, mem_mb = 30000, time = 2100
    shell: 
        '''
        mkdir -p {output.d}
        for f in {input.nums_filtered}; do
            echo -e "\n\n### run leafcutterMD on ${{f}} ###\n\n"
            group=$(basename -s _perind.nums.gz $f)
            # out_f={output.d}/${{group}}_perind.nums.gz
            out_prefix={output.d}/${{group}}

            workflow/scripts/leafcutter/scripts/leafcutterMD.R \
                -p {threads} \
                --output_prefix $out_prefix \
                $f &> {log}
        done
        
        txts=({output.d}/*txt)
        for txt in ${{txts[@]}}; do
              bgzip $txt
        done
        '''