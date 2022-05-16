#Lipid flip-flops
#loops over trajectories and checks if some molecule (lipid or cholesterol) in membrane flips
import os
import sys
import yaml
import json
import matplotlib.pyplot as plt
import numpy as np
import math
from random import randint

from matplotlib import cm
from scipy.stats import norm

import urllib.request
from urllib.error import URLError,HTTPError,ContentTooShortError

import MDAnalysis as mda

sys.path.insert(1, '../../Databank/Scripts/BuildDatabank/')
from databankLibrary import download_link, lipids_dict, databank

def loadMappingFile(path_to_mapping_file):
    # load mapping file into a dictionary
    mapping_dict = {}
    with open('../../Databank/Scripts/BuildDatabank/mapping_files/'+path_to_mapping_file, "r") as yaml_file:
        mapping_dict = yaml.load(yaml_file, Loader=yaml.FullLoader)
    yaml_file.close()
    
    return mapping_dict

def getLipids(readme):
        lipids = []
        for key in lipids_dict.keys():
            try:
                if key in readme['COMPOSITION'].keys():                
       #         if readme['N'+key] != [0,0]: 
                    lipids.append(key)
            except KeyError:
                continue
        return lipids
        
def getAtoms(readme, lipid):
    atoms = ""
    path_to_mapping_file = readme['COMPOSITION'][lipid]['MAPPING']
    mapping_dict = loadMappingFile(path_to_mapping_file)
    for key in mapping_dict:
        atoms = atoms + ' ' + mapping_dict[key]['ATOMNAME']
  
    return atoms

def membraneCentreOfMass(universe, readme):
    lipids = []
    for key_mol in lipids_dict:
        if key_mol in readme['COMPOSITION'].keys():
            selection = ""
            m_file = readme['COMPOSITION'][key_mol]['MAPPING']
            
            with open('../../Databank/Scripts/BuildDatabank/mapping_files/'+m_file,"r") as f:
                for line in f:
                    if len(line.split()) > 2 and "Individual atoms" not in line:
                        selection = selection + "(resname " + line.split()[2] + " and name " + line.split()[1] + ") or "
                    elif "Individual atoms" in line:
                        continue
                    else:
                        selection = "resname " + readme['COMPOSITION'][key_mol]['NAME']
                        #print(selection)
                        break
            selection = selection.rstrip(' or ')
        #    print("selection    " + selection)
            molecules = universe.select_atoms(selection)
          #  print("molecules")
          #  print(molecules)
            if molecules.n_residues > 0:
                lipids.append(universe.select_atoms(selection))
            #print(lipids) 
    # join all the selected the lipids together to make a selection of the entire membrane and calculate the
    # z component of the centre of mass of the membrane
    membrane = universe.select_atoms("")   # causes a warning: UserWarning: Empty string to select atoms, empty group returned, but it's ok
    R_membrane_z = 0
    if lipids != []:
        for i in range(0,len(lipids)):
            membrane = membrane + lipids[i]
            #print("membrane") 
            #print(membrane)  
    R_membrane = membrane.center_of_mass()[2]
  #  print(R_membrane)
    return R_membrane
    
def headgroupAtoms(readme,lipid):
    mapping_file = loadMappingFile(readme['COMPOSITION'][lipid]['MAPPING'])
    #returns atom names labeled as belonging to headgroup in the mapping file
    headgroup = ""
    for key in mapping_file:
        if mapping_file[key]['FRAGMENT'] == 'headgroup':
            headgroup = headgroup + " " + mapping_file[key]['ATOMNAME']
            
    #if names in structure file contain ' characters add escape sign
    if "\'" in headgroup:
        headgroup = headgroup.replace("'","\'")
    
    return headgroup
    

path = '../../Databank/Data/Simulations/006/559/006559139e730fc43b244726992145c2f37a1461/3c99810c45a83b4ba0e54a69fdea8817498a8930/'
db_data = databank(path)
systems = db_data.get_systems()


#compare frames every 500 ps
time_diff = 500 

for system in systems:

    indexingPath = "/".join(system['path'].split("/")[5:9])
    print(indexingPath)
    subdir = '../../Databank/Data/Simulations/' + indexingPath + '/'
    DATAdir = '../Data/Flipflops/' + indexingPath + '/'
    
    if (not os.path.isdir(DATAdir)):
        os.system('mkdir ../Data/Flipflops/' + indexingPath.split('/')[0])
        os.system('mkdir ../Data/Flipflops/' + indexingPath.split('/')[0] + '/' + indexingPath.split('/')[1])
        os.system('mkdir ../Data/Flipflops/' + indexingPath.split('/')[0] + '/' + indexingPath.split('/')[1] + '/' + indexingPath.split('/')[2])
        os.system('mkdir ../Data/Flipflops/' + indexingPath.split('/')[0] + '/' + indexingPath.split('/')[1] + '/' + indexingPath.split('/')[2] + '/' + indexingPath.split('/')[3]) 
        
    flipflop_file = DATAdir + 'flipflop.dat' #change later
                                      
                
    doi = system['DOI']
    trj=system['TRJ'][0][0]
    tpr=system['TPR'][0][0]
    trj_name = subdir +  trj
    tpr_name = subdir + tpr
    trj_url = download_link(doi, trj)
    tpr_url = download_link(doi, tpr)
    EQtime=float(system['TIMELEFTOUT'])*1000

    if (not os.path.isfile(tpr_name)):
        response = urllib.request.urlretrieve(tpr_url, tpr_name)
                        
    if (not os.path.isfile(trj_name)):
        response = urllib.request.urlretrieve(trj_url, trj_name)
                    
    xtcwhole=subdir + '/whole.xtc'
    if (not os.path.isfile(xtcwhole)):
        os.system('echo System | gmx trjconv -f ' + trj_name + ' -s ' + tpr_name + ' -o ' + xtcwhole + ' -pbc mol -b ' + str(EQtime))
                    
    u = mda.Universe(tpr_name,xtcwhole)
                
    lipids = getLipids(system)
   # print(lipids)            
    end = len(u.trajectory)
                
    dt = int(u.trajectory.dt)
    
    start = int(EQtime / dt)
    flipflops = 0
    
    # compare frames i and i - 500 ps
    # if saving frequency is greater than 500 ps then there's no need to skip
    skip = int(time_diff / dt)
    if time_diff < dt:
        skip = 1
    
 #   print("EQtime " + str(EQtime))
 #   print("start " + str(start))
 #   print("dt " + str(dt))
 #   print("end " + str(end))
 #   print("skip " + str(skip))
 #   print("time_diff  " + str(time_diff))
                
    frames_lipids = []
                
    for lipid in lipids:
        #save state of the leaflet in previous frames
        previous_leaflet = []
        lipid_atoms = getAtoms(system, lipid)
        
        headgroup_atoms = headgroupAtoms(system,lipid)

        lipid_resname = system['COMPOSITION'][lipid]['NAME']

        headgroup_atoms_selection = u.select_atoms('resname ' + lipid_resname + ' and name' + headgroup_atoms).split('residue')  
     
        
     #   print(lipid + "  " +lipid_resname)
        
     #   print(lipid_atoms)
     #   print(lipid_atoms_selection)
     
        #every 500th frame for checking flip flops
        
        for ts in u.trajectory[start:end:skip]:
          #  print("frame " + str(ts.frame))
            time = u.trajectory.time
            
            R_m = membraneCentreOfMass(u, system)
            
            l_leaflet_fr = []
            
            for i, hg in enumerate(headgroup_atoms_selection):   #!!!!
                hg_z = hg.center_of_mass()[2]
                hg_resid = i
                if hg_z < R_m:
                    l_leaflet_fr.append([hg_resid, 'l1']) #in leaflet 1
                if hg_z > R_m:
                    l_leaflet_fr.append([hg_resid,'l2']) #in leaflet 2       

                               
           # print("l_leaflet_fr")
           # print(l_leaflet_fr)
            
            #check if lipid has changed leaflet between frames
            if ts.frame != start and ts.frame != start + skip and ts.frame != start + 2*skip: # first two frames
             #   print(previous_leaflet)
                for j in range(0,len(l_leaflet_fr)-1):
                
                    index_1 = int(ts.frame/(skip)-start)-1 
                    index_2 = int(ts.frame/(skip)-start)-2 
                    
                   # print(previous_leaflet[index_1])
                   # print(previous_leaflet[index_2])
                    
                    
                    # is lipid location in current frame different from the lipid location in the two previous frames  
                    if l_leaflet_fr[j][0] != previous_leaflet[index_1][j][0] and l_leaflet_fr[j][1] == previous_leaflet[index_1][j][1] and l_leaflet_fr[j][0] != previous_leaflet[index_2][j][0] and l_leaflet_fr[j][1] == previous_leaflet[index_2][j][1]: # RESID SAME BUT LEAFLET IS DIFFERENT
                        #print("FLIP FLOP")
                        flipflops += 1
                        res_index = l_leaflet_fr[j][1]
                        frames_lipids.append([str(time - time_diff) + ' ' + str(time) + '   ' + lipid + ' ' + str(res_index) + ' ' + str(ts.frame - skip) + ' ' + str(ts.frame)    ])
                                    
                        
                                    
                        
            previous_leaflet.append(l_leaflet_fr.copy())    
                            
       # print(flipflops) 
       
    time_btwn_frames = time_diff
    
    if time_diff < dt:
        time_btwn_frames = dt   
       
    with open(flipflop_file, 'w') as f:
         f.write('# time n - ' + str(time_btwn_frames) + ' (ps)    time n (ps)     molecule\n')
                        
    f.close()               
                       
                            
    if flipflops > 0:
        print('number of flip flops in trajectory is ' + str(flipflops) + " in the simulation in directory " + indexingPath)

       # save frame nr and lipid name
        with open(flipflop_file, 'a') as f:
            for fr in frames_lipids:
                f.write(str(fr).strip("[]'") + '\n')
        f.close()
    else:
        print("no flip flops in trajectory")


                                    
                
                
