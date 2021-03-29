def import_hmm_df(hmm_file,plasmid_info,database_name):

    import pandas as pd

    mod_hmm_file = hmm_file.split('.')[:-1]
    mod_hmm_file.append('_mod.out')
    mod_hmm_file = ''.join(mod_hmm_file)

    headers = ['target name', 'target_accession',  'query name', 
           'query_accession',    'E-value_all',  'score_all',  'bias_all',
           'E-value', 'score',  'bias', 'exp', 'reg', 'clu',  'ov', 'env',
           'dom', 'rep', 'inc']

    n_keep = len(headers)

    with open(hmm_file, 'r') as file:
        data = file.readlines()
        for i in range(len(data)):
            data[i] = ' '.join(data[i].split())
            data[i] = ','.join(data[i].split()[0:n_keep])
        file.close()

    new_data = [item for item in data if (item[0] != '#')]
    
    with open(mod_hmm_file,"w") as mod_file:
        for item in new_data:
            mod_file.write("%s\n" % item)

    hmm_df = pd.read_csv(mod_hmm_file, sep = ',', index_col=None,
                       names = headers)
    
    #Get base names of target without frame info
    true_name = [ '_'.join(x.split('_')[:-1]) for x in hmm_df['target name'] ] 
    hmm_df['true_target'] = true_name

    #Sort by results with Evalue below cutoff
    hmm_df = hmm_df[hmm_df['E-value_all'] < 1e-30]

    plasmid_info[database_name] = ''
    for i in range(len(plasmid_info['GenBank-Accn'])):
        target = plasmid_info.loc[i,'GenBank-Accn']
        matching_results = hmm_df[hmm_df['true_target']==target]
        plasmid_info[database_name][i] = matching_results['query name'].tolist()

    return(plasmid_info)


