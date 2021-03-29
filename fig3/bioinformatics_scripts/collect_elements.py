

from Bio import SeqIO
import pandas as pd

#Get pickled element dataframe
genome_dir = '/tigress/DONIA/data/NCBI_complete_bacterial_genomes/'
final_element_file = '/tigress/DONIA/jglopez/ParasiticEcologyComp/collect_elements/all_elements.fasta'
element_pkl = '/tigress/DONIA/jglopez/ParasiticEcologyComp/aggregate_assembly_reports/element_df.pkl'
element_df = pd.read_pickle(element_pkl)

#Get all unique origins
unique_origins = set(element_df['Origin'])

dna_prefix = '_genomic.fna'

record_set = [];

#Loop through origins and add all 
for origin in unique_origins:
  
    #Get elements from that origin
    origin_df = element_df.loc[element_df['Origin']==origin,:]

    #Define genomic fasta file of genome
    genome_fasta = genome_dir + origin +'/' + origin + dna_prefix

    #Get genbank ids of element
    genbank_ids = origin_df['GenBank-Accn'].tolist()

    #Loop through sequences and add those matching element ids
    for record in SeqIO.parse(genome_fasta,'fasta'):
        if any(gid in record.description for gid in genbank_ids):       
           record_set.append(record)
    
SeqIO.write(record_set,final_element_file,'fasta')
