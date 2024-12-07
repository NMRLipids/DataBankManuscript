import os
import sys
import numpy as np

import yaml
import json
import matplotlib.pyplot as plt
import mdtraj
import urllib.request
import seaborn as sns

import MDAnalysis as mda
#from MDAnalysis.analysis.waterdynamics import MeanSquareDisplacement as MSD

sys.path.insert(1, '../../Databank/Scripts/BuildDatabank/')
from databankLibrary import * # download_link, read_mapping_file, read_mapping_file_res, read_mapping_filePAIR, make_positive_angles, lipids_dict, molecules_dict, databank

path = '../../Databank/Data/Simulations/'
db_data = databank(path)
systems = db_data.get_systems()

molecules = []
for key in lipids_dict:
    molecules.append(key)
for key in molecules_dict:
    molecules.append(key)


cwd = os.getcwd()
# Loop over simulations
for system in systems:

    # Extracting information from README.yaml file

    if 'gromacs' not in system['SOFTWARE']:
        continue

    try:
        if 'WARNINGS' in system and 'AMBIGUOUS_ATOMNAMES' in system['WARNINGS']:
            continue

        if 'WARNINGS' in system and 'ORIENTATION' in system['WARNINGS']:
            continue

        if 'WARNINGS' in system and 'GROMACS_VERSION' in system['WARNINGS'] and system['WARNINGS']['GROMACS_VERSION'] == 'gromacs3':
            continue
    except:
        pass

        
    subdir = '../../Databank/Data/Simulations/' + system['path']
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
    outputFOLDERS = subdir.replace("../../Databank/Data/Simulations/","../Data/MD-PERMEATION/")    
    #outfilename = str(densityFOLDERS) + '/WATERlateralMSD.xvg'

    CheckOutPutFile = outputFOLDERS + '/Counting_events.txt'
    print(CheckOutPutFile)
    if os.path.isfile(CheckOutPutFile): # and os.path.getsize(CheckOutPutFile) > 0:
        print('Result already found')
        continue
    

    
    os.system('mkdir -p ' + outputFOLDERS)
    os.system('cp ' + READMEfilepath + ' ' + outputFOLDERS)

    #print(outfilename)                    
    #if os.path.isfile(outfilename):
    #    print('Result found from ' + subdir)
    #    continue

    print('Analyzing in ' + subdir)
                        
    # Downloading the tpr and xtc files locally using the information from the README.yaml file
    if (not os.path.isfile(tpr_name)):
        response = urllib.request.urlretrieve(tpr_url, tpr_name)
                        
    if (not os.path.isfile(trj_name)):
        response = urllib.request.urlretrieve(trj_url, trj_name)

    # Creating a gro file (this works only for Gromacs simulations)
    if (not os.path.isfile(gro_name)):
        os.system('echo System | gmx trjconv -f ' + trj_name + ' -s ' + tpr_name + ' -dump 0 -o ' + gro_name)

    # Creating a trajectory with molecules whole  (this works only for Gromacs simulations) 
    # Note that this leaves out the equilibration time
    xtcwhole=subdir + '/centered.xtc'
    if (not os.path.isfile(xtcwhole)):
        #continue
        os.system('echo System | gmx trjconv -f ' + trj_name + ' -s ' + tpr_name + ' -o ' + xtcwhole + ' -pbc mol -b ' + str(EQtime))
                        
    # Reads the name of water oxygen atom using the information in mapping file and README.yaml 
    water_mapping_file = '../../Databank/Scripts/BuildDatabank/mapping_files/' + system['COMPOSITION']['SOL']['MAPPING']
    waterO = read_mapping_file(water_mapping_file,'M_O_M')

    # Reads the phosphate atom name
    for lipid in system['COMPOSITION']:
        if ('CHOL' in lipid) or ('CER' in lipid) or ('DHMDMAB' in lipid) or ('DMTAP' in lipid) or ('DOG' in lipid) or ('TOCL' in lipid) or ('GB3' in lipid) or ('GM1' in lipid) or ('TLCL_0H' in lipid):
            continue
        elif lipid in lipids_dict:
            lipid_mapping_file = '../../Databank/Scripts/BuildDatabank/mapping_files/' + system['COMPOSITION'][lipid]['MAPPING']
            print(lipid_mapping_file)
            P = read_mapping_file(lipid_mapping_file,'M_G3P2_M')

    # Creating OW.gro and P.gro this works only for Gromacs simulations) 
    ndxSTR = 'keep 0 \na ' + waterO  + ' \nkeep 1 \na' + P +  ' \nq'
    if (not os.path.isfile(outputFOLDERS + '/OWPindex.ndx')):
        os.system('echo "' + ndxSTR + '" | gmx make_ndx -f ' + tpr_name + ' -o ' + outputFOLDERS + '/OWPindex.ndx')
    if (not os.path.isfile(outputFOLDERS + '/OW.gro')):
        os.system('echo "' + waterO + '" | gmx trjconv -f ' + xtcwhole + ' -s ' + tpr_name + ' -n  ' + outputFOLDERS + '/OWPindex.ndx ' + ' -o ' + outputFOLDERS + '/OW.gro')
    if (not os.path.isfile(outputFOLDERS + '/P.gro')):
        os.system('echo "' + P + '" | gmx trjconv -f ' + xtcwhole + ' -s ' + tpr_name + ' -n  ' + outputFOLDERS + '/OWPindex.ndx ' + ' -o ' + outputFOLDERS + '/P.gro')

    # Reading estimated center, simulation time, precision, and the number of the first water residue

    # center of mass
    try:
        u = mda.Universe(gro_name,xtcwhole)
    except:
        print('MDanalysis failed')
        continue
    
    lipids = getLipids(system)
    c = u.select_atoms(lipids)
    ctom = c.atoms.center_of_mass()[2]
    ctom = ctom * 0.1  # change units to nm
    
    # simulation time (ns)
    time = int(system['TRJLENGTH']*0.001 - system['TIMELEFTOUT'])
    #print(time)

    # precision
    # print("'grep "t="' + outputFOLDERS + "/OW.gro  | head -n 2 | awk '{if(NR ==1 ) v1 = $8; if(NR == 2) v2 = $8}END{print v2-v1}' > "  + outputFOLDERS + 'precision.dat'")
    #os.system('grep "t="' + outputFOLDERS + "/OW.gro  | head -n 2 | awk '{if(NR ==1 ) v1 = $8; if(NR == 2) v2 = $8}END{print v2-v1}' > "  + outputFOLDERS + 'precision.dat')
    #with open(outputFOLDERS + 'precision.dat', "r") as presFILE:
    #    precision = presFILE.read()
    #print(precision)

    precision = int(u.trajectory.dt)
    print(precision)

    
    # first water residue
    os.system("awk 'FIELDWIDTHS = '5' {if(NR == 3){ printf $1; exit}}' " +  outputFOLDERS + '/OW.gro > ' + outputFOLDERS + '/Wresidue.dat')
    with open(outputFOLDERS + 'Wresidue.dat', "r") as resFILE:
        Wresidue = resFILE.read()

    os.system('cp ./MD-permeation ' + outputFOLDERS)
    print('Changing to', outputFOLDERS)
    os.chdir(outputFOLDERS)
    print('./MD-permeation ' + str(ctom) + ' ' + str(time) +  ' ' + str(precision) + ' ' + str(Wresidue))
    os.system('./MD-permeation ' + str(ctom) + ' ' + str(time) +  ' ' + str(precision) + ' ' + str(Wresidue))
    print('Changing to', cwd)
    os.chdir(cwd)

    os.system("rm " +  outputFOLDERS + '/OW.gro')
    os.system("rm " +  outputFOLDERS + '/P.gro')
    
