# This Script is for the sub-Project of markus
# Concretely this is for the analyses of the enrichments.
# This containes the functions to perform all the analyses.

# Libraries:

from goatools import obo_parser
import pandas as pd

#############
# Functions #
#############

# Funtion to load the data from the files:

network      = "systematic_PPI_BioGRID"
technique = "GCV-all"
enrichemnt= "BP"
total_run   = 1

technique  = "GCV-O+"
technique  = "triangle"

def Load_Data(network, technique, enrichemnt, total_run, cluster = [43, 63]):
    
    # Load the user options (enrichment for GO terms and distance technique)
    
    NETWORK = network
    TECH = technique
    GO   = enrichemnt
    
    # Load the files (by repetition and cluster):
    # For repetitions and cluster we create the path to the file
    
    # Final dataframe to save all the info from inside of the loop:
    
    Final_df = pd.DataFrame(columns = ["GO_term", "Enriched_in", "Total_K", "Run"])
    
    # Loop start:
    
    for run in range(total_run):
        
        RUN = run
        print(run)
        
        for cluster in range(cluster[0], (cluster[1] + 1)):
            
            NR = cluster
            path = rf"C:\Users\Sergio\Desktop\Markus_Project\{NETWORK}\{TECH}\js_divergence\{GO}\{RUN}_{NR}.csv"
            
            # Load the file using pandas:
            
            Enrichment_df = pd.read_csv(path)
            
            # Checking in which cluster a GO term is enriched:
            
            for GO_term in range(len(Enrichment_df)):
                itera_GO = Enrichment_df.iloc[GO_term]
                # Is enriched at least in one cluster ?
                if sum(itera_GO == True) != 0:
                    # In which clusters is enriched:
                    
                    cluster_in = itera_GO[itera_GO == True].index.values.astype(int)
                    
                    # Create the itera information:
                    # Information:
                    
                    Name_GO          = Enrichment_df.iloc[GO_term,0]
                    Cluster_Enriched = cluster_in
                    Option_K         = NR
                    Repetition       = RUN
                    
                    Itera_Info = pd.DataFrame({"GO_term"     : [Name_GO],
                                               "Enriched_in" : [Cluster_Enriched],
                                               "Total_K"     : [Option_K],
                                               "Run"         : [Repetition]})
                    # Save the info:
                    
                    Final_df = Final_df.append(Itera_Info)
        
        # Do the level Analyses:
        
        save_path = rf"C:\Users\Sergio\Desktop\Markus_Project\{TECH}_{GO}_{NETWORK}.csv"
        obo_file_path = rf'C:\Users\Sergio\Desktop\Markus_Project\go-basic.obo'
        GO_Level_Function(Final_df, obo_file_path, save_path)
                    
# Function to get the GO term levels from the dataframe of enrichments:
# This function uses the obo Parser for the analyses.
               
def GO_Level_Function(Final_df, obo_file_path, save_path):
    
    # Load the obo data to the parser:
    
    go_dag = obo_parser.GODag(obo_file_path)
    
    # Get the levels of the GO terms: 
    
    Final_df_levels = Final_df
    
    # We don´t know the level of these terms (are not in the parser):
     
    missing_GO_terms = set(Final_df_levels["GO_term"]) - set(go_dag)
    
    print("You have a total of " + str(len(missing_GO_terms)) + " without level info")
    
    # Delete the terms to avoid errors:
    
    Final_df_levels = Final_df_levels[~Final_df_levels.GO_term.isin(missing_GO_terms)]
    Final_df_levels_deleted = Final_df_levels[Final_df_levels.GO_term.isin(missing_GO_terms)]
    
    # Get the level of the terms that we can:

    Final_df_levels['Level'] = [go_dag[term].level for term in Final_df_levels['GO_term']]
    
    # Build the final dataFrame with a NaN for those that we don´t have the level:
    # For future studies the users can decide if they include them or not in the analyses.
    
    Final_df_Levels = pd.concat([Final_df_levels, Final_df_levels_deleted])
    Final_df_Levels = Final_df_Levels.sort_values(by=['Total_K'])
    Final_df_Levels = Final_df_Levels.reset_index(drop = True)
    
    # Save the results:
    
    Final_df_Levels.to_csv(save_path, sep = "\t", header = True, index = False)
    
# Funtion to compare the GO terms (dividing by Intersection and Unique)
# Pair-wise based method
# This function is used as the based for the rest of the experiments.
    
def Pair_Wise_GO_Comparison(network, enrichemnt, total_run, cluster = [43, 63]):
    
    # Data from the user:
   
    NETWORK = network  
    GO = enrichemnt

    # Total Techniques:
    
    technique_1 = "GCV-all"
    technique_2 = "GCV-O+"
    technique_3 = "triangle"
    
    # Paths to the files:
    
    path_1 = rf"C:\Users\Sergio\Desktop\Markus_Project\{technique_1}_{GO}_{NETWORK}.csv"
    path_2 = rf"C:\Users\Sergio\Desktop\Markus_Project\{technique_2}_{GO}_{NETWORK}.csv"
    path_3 = rf"C:\Users\Sergio\Desktop\Markus_Project\{technique_3}_{GO}_{NETWORK}.csv"
    
    # Load Data:
    
    technique_1_df = pd.read_csv(path_1, sep = "\t")
    technique_2_df = pd.read_csv(path_2, sep = "\t")
    technique_3_df = pd.read_csv(path_3, sep = "\t")
    
    # Pair-wise and Total Comparisons:
    # We classify each GO term as Unique or Intersection.
    # We add two columns: one to know what we are comparing and another with the info of
    # in which technique is found.
    
    # Data_Frame for the final result:
    
    Final_Result = pd.DataFrame(columns = ["GO_Term", "Level", "Comparison", "In",
                                           "Total_K",  "Run"])
    
    for run in range(total_run):
        
        # Subset of the GO terms enriched in each run:
        
        technique_1_df_rep = technique_1_df[technique_1_df.Run == run]
        technique_2_df_rep = technique_2_df[technique_2_df.Run == run]
        technique_3_df_rep = technique_3_df[technique_3_df.Run == run]
        
        # Subset of GO terms by total number of clusters:
        
        for cluster in range(cluster[0], (cluster[1] + 1)):
            
            technique_1_df_rep_clust = technique_1_df_rep[technique_1_df_rep.Total_K == cluster]
            technique_2_df_rep_clust = technique_2_df_rep[technique_2_df_rep.Total_K == cluster]
            technique_3_df_rep_clust = technique_3_df_rep[technique_3_df_rep.Total_K == cluster]
            
            # Calculate the sets:
            
            set1 = set(technique_1_df_rep_clust.GO_term)
            set2 = set(technique_2_df_rep_clust.GO_term)
            set3 = set(technique_3_df_rep_clust.GO_term)
            
            # Comparisons:
            # 4 (three techniques, three pairwise comparions + all)
        
            for comparison in range(4):
                
                # PairWise:
                if comparison < 3:
                
                    # technique_1 Vs technique_2
                    if comparison == 0:
                    
                        comparison_1     = set1
                        comparison_2     = set2
                        technique_comp_1 = technique_1
                        technique_comp_2 = technique_2
                        data_1 = technique_1_df_rep_clust
                        data_2 = technique_2_df_rep_clust
                 
                    # technique_1 Vs technique_3
                    if comparison == 1:
                    
                        comparison_1     = set1
                        comparison_2     = set3
                        technique_comp_1 = technique_1
                        technique_comp_2 = technique_3
                        data_1 = technique_1_df_rep_clust
                        data_2 = technique_3_df_rep_clust
                
                    # technique_2 Vs technique_3
                    if comparison == 2 :
                    
                        comparison_1     = set2
                        comparison_2     = set3
                        technique_comp_1 = technique_2
                        technique_comp_2 = technique_3
                        data_1 = technique_2_df_rep_clust
                        data_2 = technique_3_df_rep_clust
                    
                    # Create the sets for each of the two techniques:
                    
                    union_tech    = comparison_1.intersection(comparison_2)
                    unique_tech_1 = comparison_1 - union_tech
                    unique_tech_2 = comparison_2 - union_tech
                    
                    # Get the GO terms:
                    
                    union_GO  = data_1[data_1.GO_term.isin(union_tech)][["GO_term", "Level"]]
                    tech_1_GO = data_1[data_1.GO_term.isin(unique_tech_1)][["GO_term", "Level"]]
                    tech_2_GO = data_2[data_2.GO_term.isin(unique_tech_2)][["GO_term", "Level"]]
                    GO_terms   = pd.concat([union_GO, tech_1_GO, tech_2_GO]).reset_index(drop = True)
                    
                    # Create some variables to keep the information:
                    
                    Comparison = pd.Series(technique_comp_1 + "_VS_" + technique_comp_2 ).repeat(len(GO_terms)).reset_index(drop = True)
                    In_Which   = pd.concat([pd.Series("Both").repeat(len(union_GO)), 
                                            pd.Series(technique_comp_1).repeat(len(tech_1_GO)),
                                            pd.Series(technique_comp_2).repeat(len(tech_2_GO))]).reset_index(drop = True)
                    
                    Cluster_it = pd.Series(cluster).repeat(len(GO_terms)).reset_index(drop = True)
                    Run_it = pd.Series(run).repeat(len(GO_terms)).reset_index(drop = True)
                    
                    # Put all together:
                    
                    Result_itera_Comp = pd.DataFrame({"GO_Term"    : GO_terms.GO_term,
                                                      "Level"      : GO_terms.Level,
                                                      "Comparison" : Comparison,
                                                      "In"         : In_Which,
                                                      "Total_K"    : Cluster_it,
                                                      "Run"        : Run_it})
                    # Save the info:
                    
                    Final_Result = Final_Result.append(Result_itera_Comp)
 
                # All Comparisons:
                if comparison == 3:
                    
                    # In this case we consider the intersection between all three sets:
                    
                    union_tech    = set1.intersection(set2, set3)
                    unique_tech_1 = set1 - union_tech
                    unique_tech_2 = set2 - union_tech
                    unique_tech_3 = set3 - union_tech
                    
                    # Take the GO terms per each technique:

                    union_GO  = technique_1_df_rep_clust[technique_1_df_rep_clust.GO_term.isin(union_tech)][["GO_term", "Level"]]
                    tech_1_GO = technique_1_df_rep_clust[technique_1_df_rep_clust.GO_term.isin(unique_tech_1)][["GO_term", "Level"]]
                    tech_2_GO = technique_2_df_rep_clust[technique_2_df_rep_clust.GO_term.isin(unique_tech_2)][["GO_term", "Level"]]
                    tech_3_GO = technique_3_df_rep_clust[technique_3_df_rep_clust.GO_term.isin(unique_tech_3)][["GO_term", "Level"]]
                    
                    # Concat all the GO terms:
                    
                    GO_terms   = pd.concat([union_GO, tech_1_GO, tech_2_GO, tech_3_GO ]).reset_index(drop = True)
                    
                    # Get information about the GO terms:
                    
                    # What are we comparing?
                    
                    Comparison = pd.Series("All").repeat(len(GO_terms)).reset_index(drop = True)
                    
                    # In which technique is unique:
                    
                    In_Which   = pd.concat([pd.Series("Three").repeat(len(union_GO)), 
                                            pd.Series(technique_1).repeat(len(tech_1_GO)),
                                            pd.Series(technique_2).repeat(len(tech_2_GO)),
                                            pd.Series(technique_3).repeat(len(tech_3_GO))]).reset_index(drop = True)
                    
                    Cluster_it = pd.Series(cluster).repeat(len(GO_terms)).reset_index(drop = True)
                    Run_it = pd.Series(run).repeat(len(GO_terms)).reset_index(drop = True)
                    
                    Result_itera_Comp = pd.DataFrame({"GO_Term"    : GO_terms.GO_term,
                                                      "Level"      : GO_terms.Level,
                                                      "Comparison" : Comparison,
                                                      "In"         : In_Which,
                                                      "Total_K"    : Cluster_it,
                                                      "Run"        : Run_it})
                    # Save the info:
                    
                    Final_Result = Final_Result.append(Result_itera_Comp)
    
    # Return the final data frame:
                    
    return(Final_Result) 





















         
                
                

        