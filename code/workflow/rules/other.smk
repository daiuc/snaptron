

# this function only look up data_source, e.g. GTEx
# returns column index of samples in the big resources/snaptron-yil/run/full/All_perind.constcounts.gz file
def getMakeRareSplicingTrainingData_Fields(wildcards):
    Columns = SampleLookupTable.query(
        'ds == @wildcards.data_source' 
    )['counts_col_id'].reset_index(drop=True).unique()
    Columns = [str(i) for i in Columns]
    return ','.join(Columns)

rule MakeRareSplicingTrainingData:
    input: 'resources/snaptron-yil/run/full/All_perind.constcounts.gz'
    # input: 'test.gz'
    output: 'results/raresplicing_training_data/{data_source}_all.tsv.gz'
    params:
        Fields = getMakeRareSplicingTrainingData_Fields,
        Min_Reads = 30
    shell:
        '''
        zcat {input} | cut -d " " -f "1,{params.Fields}" | \
            awk -v Min_Reads={params.Min_Reads} 'BEGIN {{OFS="\t"}}
                {{  
                    if (NR == 1) print $1, "N_used", "ratio_used"
                    else 
                    {{  
                        N_used=0
                        for (i=2; i<=NF; i++) 
                        {{
                            split($i, arr, "/")
                            if (arr[1] > Min_Reads) N_used = N_used + 1
                            else N_used = N_used
                        }}
                        
                        if (N_used > 0) print $1, N_used, N_used/(NF-1)
                    }}
                }}' | bgzip -c > {output} 
        '''
# note ` zcat somefile | head ` will not work with snakemake, because snakemake work in strict bash mode
# commands like so will result a broken pipe, which results in error for snakemake
# thus for test run, make a shorter but full file, and do not use `head`