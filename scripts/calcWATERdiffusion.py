#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import numpy as np

import yaml
import json
import matplotlib.pyplot as plt
import mdtraj
import urllib.request
import seaborn as sns

import MDAnalysis
from MDAnalysis.analysis.waterdynamics import MeanSquareDisplacement as MSD


sys.path.insert(1, '../../Databank/Scripts/BuildDatabank/')
from databankLibrary import download_link, read_mapping_file, read_mapping_file_res, read_mapping_filePAIR, make_positive_angles, lipids_dict, molecules_dict, databank

path = '../../Databank/Data/Simulations/'
db_data = databank(path)
systems = db_data.get_systems()

molecules = []
for key in lipids_dict:
    molecules.append(key)
for key in molecules_dict:
    molecules.append(key)

    
#colors = {'POPC' :'black','POPS':'red','POPE':'blue','POPG':'green'}

#h = []


#for subdir, dirs, files in os.walk(r'../../Data/Simulations/'):
#    for filename in files:
#        filepath = subdir + os.sep + filename
#        if filepath.endswith("README.yaml"):
#            READMEfilepath = subdir + '/README.yaml'
#            with open(READMEfilepath) as yaml_file:

for system in systems:

#                readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
#                #for molname in molecules:
#                doi = readme.get('DOI')
#                trj = readme.get('TRJ')
#                tpr = readme.get('TPR')
#                trj_name = subdir + '/' + readme.get('TRJ')[0][0]
#                tpr_name = subdir + '/' + readme.get('TPR')[0][0]
#                gro_name = subdir + '/conf.gro'
#                trj_url = download_link(doi, trj[0][0])
#                tpr_url = download_link(doi, tpr[0][0])
#                EQtime=float(readme.get('TIMELEFTOUT'))*1000

#                densityFOLDERS = subdir.replace("Simulations","WATERdiffusion")    
#                outfilename = str(densityFOLDERS) + '/WATERlateralMSD.xvg'


    # Extracting information from README.yaml file
    subdir = system['path']
    READMEfilepath = subdir + '/README.yaml'
    doi = system['DOI']
    trj = system.get('TRJ')
    tpr = system.get('TPR')
    trj_name = subdir + '/' + system.get('TRJ')[0][0]
    tpr_name = subdir + '/' + system.get('TPR')[0][0]
    gro_name = subdir + '/conf.gro'
    trj_url = download_link(doi, trj[0][0])
    tpr_url = download_link(doi, tpr[0][0])
    EQtime=float(system.get('TIMELEFTOUT'))*1000

    print(subdir)

    # Create folders where water diffusion results will be saved
    # This follows the organization of simulations in the original databank
    outputFOLDERS = subdir.replace("../../Databank/Data/Simulations/","../Data/WATERdiffusion/")    
    outfilename = outputFOLDERS + '/WATERlateralMSD.xvg'
    
    #CheckOutPutFile = outputFOLDERS + '/WATERlateralMSD.xvg'
    if os.path.isfile(outfilename): # and os.path.getsize(CheckOutPutFile) > 0:
        print('Result already found')
        continue

    os.system('mkdir -p ' + outputFOLDERS)
    os.system('cp ' + READMEfilepath + ' ' + outputFOLDERS)
    
    print('Analyzing in ' + subdir)
                        
    #Download tpr and xtc files to same directory where dictionary and data are located


    if (not os.path.isfile(tpr_name)):
        response = urllib.request.urlretrieve(tpr_url, tpr_name)
                        
    if (not os.path.isfile(trj_name)):
        response = urllib.request.urlretrieve(trj_url, trj_name)

    # Creating a gro file (this works only for Gromacs simulations)
    if (not os.path.isfile(gro_name)):
        os.system('echo System | gmx trjconv -f ' + trj_name + ' -s ' + tpr_name + ' -dump 0 -o ' + gro_name)

    # Creating a trajectory with molecules whole  (this works only for Gromacs simulations) 
    # Note that this leaves out the equilibration time
    xtcwhole=subdir + '/whole.xtc'
    if (not os.path.isfile(xtcwhole)):
        os.system('echo System | gmx trjconv -f ' + trj_name + ' -s ' + tpr_name + ' -o ' + xtcwhole + ' -pbc mol -b ' + str(EQtime))

        
        
    #u = MDAnalysis.Universe(gro_name, trj_name)
    #water_mapping_file = '../BuildDatabank/mapping_files/'+readme['MAPPING_DICT']['SOL']
    water_mapping_file = '../../Databank/Scripts/BuildDatabank/mapping_files/' + system['COMPOSITION']['SOL']['MAPPING']
    waterO = read_mapping_file(water_mapping_file,'M_O_M')
    #waterRES = readme.get('SOL')
                    
    #os.system('echo System | gmx trjconv -f ' + xtc + ' -s ' + tpr + ' -o ' + xtcwhole + ' -pbc mol -b ' + str(EQtime))
    
    ndxSTR = 'keep 0 \na ' + waterO  + ' \nkeep 1 \nq'
    #print(ndxSTR)
    os.system('echo "' + ndxSTR + '" | gmx make_ndx -f ' + tpr_name + ' -o ' + outputFOLDERS + '/SOLindex.ndx') 
    os.system('gmx msd -f ' + xtcwhole + ' -s ' + tpr_name + ' -n  ' + outputFOLDERS + '/SOLindex.ndx -o ' + outfilename + ' -trestart 1000 -lateral z')
                    
