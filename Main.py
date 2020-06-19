# MAIN SCRIPT

"""
This script computes all the biological experiments. To run it is necessary to
load the Function_Files script that contains all the functions.
"""
import os
import multiprocessing
from multiprocessing import Pool

# Set work directory:

os.chdir(r"C:\Users\Sergio\Desktop\Markus_Project")
import Functions_Files

######################
# Yeast PPI Network: #
######################

# Arguments:

network      = "systematic_PPI_BioGRID"
technique_1  = "GCV-all"
technique_2  = "GCV-O+"
technique_3  = "triangle"
technique_4  = "GDV"
enrichemnt_1 = "BP"
enrichemnt_2 = "CC"
enrichemnt_3 = "MF"
total_run   = 10

# 1. Prepare the Data: #

# Technique "GCV_all":
# BP, CC, and MF:

sub_command1 = (network, technique_1, enrichemnt_1, 10)
sub_command2 = (network, technique_1, enrichemnt_2, 10)
sub_command3 = (network, technique_1, enrichemnt_3, 10)

# Technique "GCV-O+":
# BP, CC, and MF:

sub_command4 = (network, technique_2, enrichemnt_1, 10)
sub_command5 = (network, technique_2, enrichemnt_2, 10)
sub_command6 = (network, technique_2, enrichemnt_3, 10)

# Technique "triangle":
# BP, CC, and MF:

sub_command7 = (network, technique_3, enrichemnt_1, 10)
sub_command8 = (network, technique_3, enrichemnt_2, 10)
sub_command9 = (network, technique_3, enrichemnt_3, 10)

# Technique "GDV":
# BP, CC, and MF:

sub_command10 = (network, technique_4, enrichemnt_1, 10)
sub_command11 = (network, technique_4, enrichemnt_2, 10)
sub_command12 = (network, technique_4, enrichemnt_3, 10)


# Run the code:

Process = [sub_command1,sub_command2, sub_command3, sub_command4, sub_command5,
           sub_command6, sub_command7,sub_command8,sub_command9, sub_command10,
           sub_command11, sub_command12]

for arguments in Process:
    Functions_Files.Load_Data(*arguments)

# 2. Prepare the PairWise Comparison (for the four techniques):
    
# For BP, CC, and MF:
# For each annotation all the four different techniques are compared:
    
sub_command1 = (network, enrichemnt_1, 10)
sub_command2 = (network, enrichemnt_2, 10)
sub_command3 = (network, enrichemnt_3, 10)

Process = [sub_command1,sub_command2, sub_command3]

for arguments in Process:
    Functions_Files.Pair_Wise_GO_Comparison(*arguments)
    
# 3. Create the plots for the comparisons:


# Technique "GCV_all" VS "GCV-O+":
# BP, CC, and MF:

sub_command1 = (network, technique_1, technique_2, enrichemnt_1, 10)
sub_command2 = (network, technique_1, technique_2, enrichemnt_2, 10)
sub_command3 = (network, technique_1, technique_2, enrichemnt_3, 10)
    
# Technique "GCV_all" VS "triangle":
# BP, CC, and MF:

sub_command4 = (network, technique_1, technique_3, enrichemnt_1, 10)
sub_command5 = (network, technique_1, technique_3, enrichemnt_2, 10)
sub_command6 = (network, technique_1, technique_3, enrichemnt_3, 10)
    
# Technique "GCV_all" VS "GDV":
# BP, CC, and MF:
    
sub_command4 = (network, technique_1, technique_4, enrichemnt_1, 10)
sub_command5 = (network, technique_1, technique_4, enrichemnt_2, 10)
sub_command6 = (network, technique_1, technique_4, enrichemnt_3, 10) 

# Technique "GCV-O+" VS "triangle":
# BP, CC, and MF:  

sub_command7 = (network, technique_2, technique_3, enrichemnt_1, 10)
sub_command8 = (network, technique_2, technique_3, enrichemnt_2, 10)
sub_command9 = (network, technique_2, technique_3, enrichemnt_3, 10) 

# Technique "GCV-O+" VS "GDV":
# BP, CC, and MF:   

sub_command10 = (network, technique_2, technique_4, enrichemnt_1, 10)
sub_command11 = (network, technique_2, technique_4, enrichemnt_2, 10)
sub_command12 = (network, technique_2, technique_4, enrichemnt_3, 10) 

# Technique "triangle" VS "GDV":
# BP, CC, and MF:    
    
sub_command13 = (network, technique_3, technique_4, enrichemnt_1, 10)
sub_command14 = (network, technique_3, technique_4, enrichemnt_2, 10)
sub_command15 = (network, technique_3, technique_4, enrichemnt_3, 10) 
    
Process = [sub_command1,sub_command2, sub_command3, sub_command4, sub_command5,
           sub_command6, sub_command7,sub_command8,sub_command9, sub_command10,
           sub_command11, sub_command12, sub_command13, sub_command14, sub_command15]

for arguments in Process:
    Functions_Files.Main_Plot_Function(*arguments)    
