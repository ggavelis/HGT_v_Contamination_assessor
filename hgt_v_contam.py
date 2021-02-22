#!/usr/bin/python

#Names for input files
annot_directory='/depot/jwisecav/data/ggavelis/all_dinos/8_master_annotations/'

#Names for output files
filtered_directory='/depot/jwisecav/data/ggavelis/all_dinos/9_filtered_master_annotations/'
filtered_suffix='_AI_filtered_annot.tsv'
annot_suffix='_master_annot.tsv'
stats_suffix='_AI_filter_stats.txt'

#For purposes of checking if "prokaryote" genes have a polyA tail:
prokaryote_lineages=['other_Archaea', 'Aenigmarchaeota', 'Asgard_Archaea', 'Crenarchaeota', 'Diapherotrites', 'Euryarchaeota', 'Geoarchaeota', 'Korarchaeota', 'Nanoarchaeota', 'Nanohaloarchaeota', 'Parvarchaeota', 'Thaumarchaeota', 'other_Bacteria', 'Actinobacteria', 'Aquificae', 'Armatimonadetes', 'BacteroChlorobi', 'Caldiserica', 'ChlamydiVerrucomicrobia', 'Chloroflexi', 'Chrysiogenetes', 'Cyanobacteria', 'Deferribacteres', 'DeinocThermus', 'Dictyoglomi', 'Elusimicrobia', 'FibrobacAcidobacteria', 'Firmicutes', 'Fusobacteria', 'Gemmatimonadetes', 'Nitrospinae', 'Nitrospirae', 'Planctomycetes', 'Proteobacteria', 'Spirochaetes', 'Synergistetes', 'Tenericutes', 'Thermodesulfobacteria', 'Thermotogae', 'Viroids', 'Viruses']

import os; import numpy as np; import pandas as pd; #import Bio; from Bio import SeqIO

libs=[\
"Dinophysis_norvegica_DN1_GG~Proteobacteria~Teleaulax", \
"Dinophysis_norvegica_DN4_GG~Proteobacteria~Teleaulax", \
"Amyloodinium_ocellatum0_SRA~other_Embryophyta,Acrogymnospermae,other_Magnoliophyta,monocots,eudicots,Metazoa,Ciliophora~Paraphysomonas", \
"Dinophysis_fortii0_SRA~~Mesodinium,Teleaulax", \
"Durinskia_baltica_MMETSP2~other_Stramenopiles,Bacillariophyta~", \
"Gambierdiscus_australes_MMETSP~Amoebozoa,Euglenozoa~", \
"Green_Dinoflagellate_M0_SRA~BacteroChlorobi,Proteobacteria,Planctomycytes~", \
"Green_Dinoflagellate_T0_SRA~BacteroChlorobi,Proteobacteria~", \
"Gyrodiniellum_shiwhaense0_SRA~other_Stramenopiles,Proteobacteria,BacteroChlorobi~Mesodinium,Cafeteria,Wobblia,Paraphysomonas,Phaeocystis,Pythium~", \
"Kryptoperidinium_foliaceum_MMETSP2_CCMP~Bacillariophyta~", \
"Nusuttodinium_aeruginosum0_SRA~Cryptophyta,Metazoa,other_Embryophyta,Acrogymnospermae,other_Magnoliophyta,monocots,eudicots~", \
"Oxyrrhis_marina0_MM~Haptophyceae,other_Stramenopiles,Bacillariophyta~", \
"Peridinium_bipes0_SRA~Proteobacteria,Sordariomycetes,Planctomycetes,other_Embryophyta,Acrogymnospermae,other_Magnoliophyta,monocots,eudicots~", \
"Ross_Sea_Dinoflagellate0_SRA~Haptophyceae,Proteobacteria,BacteroChlorobi~", \
"Scrippsiella_trochoidea_MMETSP3~Haptophyceae~", \
"Yihiella_yeosuensis0_SRA~Proteobacteria,BacteroChlorobi,Cercozoa,Metazoa~" \
]

for lib in libs:
    
    ###### Setup variables
    set_to_flag, set_to_unflag, set_to_remove = set(), set(), set()
    lin_list, gen_list, len_lin_list, len_gen_list = [],[],0,0
    parameters=lib.split('~'); basename=parameters[0]
    # output files
    infile=annot_directory + basename + annot_suffix; outfile=filtered_directory + basename + filtered_suffix
    outfasta=filtered_directory + basename + '.faa'; fasta_handle=open(outfasta,"w")
    outstats=filtered_directory + basename + stats_suffix; stats_handle=open(outstats,"w")
    lin_list=parameters[1].split(','); len_lin_list=str(len(lin_list)); gen_list=parameters[2].split(',')

    ###### Process annotation file
    if os.path.isfile(outfasta) != True or os.path.isfile(outfile) != True or os.path.isfile(outstats) != True:
        df = pd.read_csv(infile, delimiter='\t') # Read input file into pandas
        df['besthit'] = df['best_hit'].astype(str)
        print('*****' + basename + '*****\n')
        
        #### WHAT % AA have a dino-spliced leader?
        if 'Y' not in df['dinoflagellate_spliced_leader'].values:
            AA_with_SL = 0
        else:
            AA_with_SL = df['dinoflagellate_spliced_leader'].value_counts()['Y']
        total_AA = len(df['dinoflagellate_spliced_leader'])
        percent_AA_with_SL = AA_with_SL / total_AA * 100
        print('% AA with SL: ' + str(percent_AA_with_SL)); stats_handle.write('percent_AA_with_SL:'+str(percent_AA_with_SL)+'\n')
        
        if lin_list[0] != '': ### FILTERING BY LINEAGE
            print('Filtering by lineage'); stats_handle.write('filter_lineage:Y\n\n')
            ##### Filter by subclade_of_best_hit column to exclude hits to contaminate lineages
            df.loc[df['subclade_of_best_hit'] != lin_list[0]][0:1]
            
            for i in lin_list:
                ### Find 'subclade_of_best_hit's to == that cont lineage
                cont_df = df.loc[df['subclade_of_best_hit'] == i]
                for index, row in cont_df.iterrows():
                    set_to_flag.add(row['AA_ID'])
                num_of_that_lineage=len(cont_df)  # How many are there?
                
                ### Are any of them redeemable? (i.e. have a spliced leader?)
                cont_df_with_sl = cont_df.loc[df['dinoflagellate_spliced_leader'] == 'Y']
                for index, row in cont_df_with_sl.iterrows():
                    set_to_unflag.add(row['AA_ID'])
                print('   ' + str(num_of_that_lineage) + ' ' + i + ' proteins. Of which ' + str(len(cont_df_with_sl)) + ' have a dino SL')
                stats_handle.write('lineage:'+i+'\nnum_of_that_lineage:'+str(num_of_that_lineage)+'\nlineage_subset_with_sl:'+str(len(cont_df_with_sl))+'\n')
                
                ### if hits to 'bacteria', redeem/unflag it if it's polyadenylated
                if i in prokaryote_lineages:
                    cont_df_with_polyA = cont_df.loc[df['polyadenylated'] == 'Y']
                    print('   and with ' + str(len(cont_df_with_polyA)) + ' prokaryote proteins that have polyA tails')
                    stats_handle.write('num_with_poly:' + str(len(cont_df_with_polyA))+'\n')
                    for index, row in cont_df_with_polyA.iterrows():
                        #print(row['AA_ID'])
                        set_to_unflag.add(row['AA_ID'])
                stats_handle.write('\n')
        else:
            print('Not filtering by lineage.'); stats_handle.write('filter_lineage:N\n\n')
                    
        if gen_list[0] != '':  ### FILTERING BY GENUS OR SPECIES
            print('\nFiltering by genus/species' ); stats_handle.write('filter_genus_or_species:Y\n')
            for i in gen_list:
                
                #### Find 'best_hit' query names that contain that genus or species match
                cont_df = df.loc[df['best_hit'].str.contains(i, na=False)]
                for index, row in cont_df.iterrows():
                    set_to_flag.add(row['AA_ID'])
                num_of_that_genus=len(cont_df)  # How many are there?
                
                ### Are any of them redeemable? (i.e. have a spliced leader?)
                cont_df_with_sl = cont_df.loc[df['dinoflagellate_spliced_leader'] == 'Y']
                for index, row in cont_df_with_sl.iterrows():
                    set_to_unflag.add(row['AA_ID'])
                print('   ' + str(num_of_that_genus) + ' ' + i + ' proteins. Of which ' + str(len(cont_df_with_sl)) + ' have a dino SL')
                stats_handle.write('genus_or_species:'+i+'\nnum_of_that_taxon:'+str(num_of_that_genus)+'\ntaxon_subset_with_sl:'+str(len(cont_df_with_sl))+'\n\n')
        else:
            print('\nNot filtering by genus/species.'); stats_handle.write('filter_genus_or_species:N\n\n')
        
        ###### WHAT PROTEINS SHOULD WE REMOVE?
        
        set_to_remove = set_to_flag.difference(set_to_unflag)
        print('\n' + str(len(set_to_flag)) + ' proteins flagged as best-hits to taxa thought to be contaminating this library')
        stats_handle.write('count_AA_with_cont_taxa_as_best-hits:'+str(len(set_to_flag))+'\n')
        print('' + str(len(set_to_unflag)) + " were then unflagged based on redeeming sequence features. \n   (Dino-splice leader means contamination is unlikely.)\n   (Polyadenylation means a transcript is likely to be eukaryotic.)\n")
        print(str(len(set_to_remove)) + ' proteins should be filtered out of ' + basename)
        stats_handle.write('count_AA_that_will_be_removed:'+str(len(set_to_remove))+'\n')
        
        ## How many proteins should we keep?
        AA_count = len(df); stats_handle.write('aa_count:'+str(AA_count)+'\n')
        keep_count = AA_count - len(set_to_remove); print(str(keep_count) + ' proteins to keep'); stats_handle.write('keep_count:' + str(keep_count) + '\n')
        
        ### Make a clean dataframe
        clean_df = df[~df.AA_ID.isin(set_to_remove)] # IN pandas speak this means:
        # "clean DF = all rows of original DF where the protein in 'AA_ID' is NOT stored in 'set_to_remove.' ('~' means 'not')"
        # print(str(len(clean_df)) + ' proteins to keep \n\n')
        
        ## percent proteins filtered:
        percent_filtered=str(round(len(set_to_remove) / AA_count * 100)); print(percent_filtered + '% proteins filtered\n\n'); stats_handle.write('percent_filtered:'+percent_filtered+'\n')
        
        ###### WRITE OUTFILES
        clean_df.to_csv(outfile, sep='\t', index=False); print('Written to ' + outfile + '\n') ## Filtered annotations
        ## Filtered fasta
        for index, row in clean_df.iterrows():
            fasta_handle.write('>' + row['AA_ID'] + '\n')
            fasta_handle.write(row['AA_seq'] + '\n')

        ##### Statistics about AI filtering
        
        ###### CLEANUP
        fasta_handle.close(); stats_handle.close()
        
        ##### Write out statistics file
                    
    ######### Possible File Errors #############        
    elif os.path.isfile(infile) != True: print(infile + ' does not exist \n') # If input file doesn't exist
    elif os.path.isfile(outfile) == True: print('outfile already exists for ' + basename) # If output file already exists
    else: print(basename + "analysis failed. \n (Check if problem with input or ouptut files.)") # Other error
print('Done')
