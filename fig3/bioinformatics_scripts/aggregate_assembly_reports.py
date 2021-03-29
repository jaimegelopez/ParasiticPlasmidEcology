#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 11:32:53 2019

@author: jglopez
"""

#This scripts analyzes the assembly reports of all the genomes and extract relevant elements

import pandas as pd
import os
from io import StringIO

#Define directories for genomes and results
genome_directory = '/tigress/DONIA/data/NCBI_complete_bacterial_genomes/'        
analysis_directory = '/tigress/DONIA/jglopez/ParasiticEcologyComp/aggregate_assembly_reports/'

os.chdir(genome_directory)

#Get all genome directories
dirs = [d for d in os.listdir(genome_directory) if os.path.isdir(d)]

#Define results df
element_df = pd.DataFrame()

#Minimal size for a chromosome to be considered real
minimal_chromosome = 500000

#Loops through directories
for dir in dirs:
    file = dir + '/' + dir + '_assembly_report.txt'
    
    #Get the part of the assembly report describing molecules
    summary_file = open(genome_directory + file, 'r')
    summary_lines = summary_file.readlines()
    summary_table = [n for n in summary_lines if not n.startswith('#')]
    summary_table = "".join(summary_table)
    summary_table_IO = StringIO(summary_table)          
    summary_df = pd.read_table(summary_table_IO,
                               sep = '\t', index_col = None, header = None)

    organism = summary_lines[1].replace('# Organism name:','').replace('\n','').strip()    
    summary_df['Organism'] = organism    

    bioproject_line = [line for line in summary_lines if '# BioProject:' in line]
    bioproject = bioproject_line[0].replace('# BioProject:','').replace('\n','').strip()
    summary_df['Bioproject'] = bioproject
                       
    #Add in the origin genome to the df
    summary_df['Origin'] = dir

    #Add to dfs if criteria met    
    num_chromosome = (summary_df.iloc[:,3] == 'Chromosome').sum()
    len_chromosome = summary_df.loc[summary_df[3]=='Chromosome',8].max()

    #Add in all elements to one df
    if (num_chromosome > 0) & (len_chromosome > minimal_chromosome):
        element_df = element_df.append(summary_df)
    summary_file.close()

element_df.columns = ['Sequence-Name','Sequence-Role','Assigned-Molecule','Assigned-Molecule-Location/Type','GenBank-Accn',
                     'Relationship','RefSeq-Accn','Assembly-Unit','Sequence-Length','UCSC-style-name','Organism','Bioproject','Origin']

element_df = element_df.reset_index(drop = True)
    
element_df.to_pickle(analysis_directory + 'element_df.pkl')
    
