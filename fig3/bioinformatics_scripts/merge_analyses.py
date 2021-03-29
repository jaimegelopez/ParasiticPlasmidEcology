import pandas as pd
import merge_util

main_dir = '/tigress/DONIA/jglopez/ParasiticEcologyComp/'
element_info_loc = 'aggregate_assembly_reports/element_df.pkl'
hmm_dir = 'hmm_search/'

element_info = pd.read_pickle(main_dir + element_info_loc)

#Read hmm results
element_info =  merge_util.import_hmm_df(main_dir + hmm_dir + 'cas_results.out',element_info,'cas')

element_info.to_csv(path_or_buf = 'merged_element_info.csv',sep = '\t')
