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
    m_file = readme['COMPOSITION'][lipid]['MAPPING']
    with open('../../Databank/Scripts/BuildDatabank/mapping_files/'+m_file,"r") as f:
        for line in f:
                atoms = atoms + ' ' + line.split()[1]
  
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
    

path = '../../Databank/Data/Simulations/'
db_data = databank(path)
systems = db_data.get_systems()

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
        
    flipflop_file = DATAdir + 'flipflop.txt' #change later
    with open(flipflop_file, 'w') as f:
         f.write('# time n - ' + str(time_diff) + ' (ps)    time n (ps)     molecule\n')
                        
    f.close()                                  
                
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
    print(lipids)            
    end = len(u.trajectory)
                
    dt = u.trajectory.dt
    
    start = int(EQtime*dt)
    flipflops = 0
    
    skip = int(time_diff / dt)
                
    frames_lipids = []
                
    for lipid in lipids:
        previous_leaflet = []
        lipid_atoms = getAtoms(system, lipid)
        
        #if names in structure file contain ' characters add escape sign 
        if "\'" in lipid_atoms:
            lipid_atoms = lipid_atoms.replace("'","\'")
          #  print(lipid_atoms)
            
        lipid_resname = system['COMPOSITION'][lipid]['NAME']

        lipid_atoms_selection = u.select_atoms('resname ' + lipid_resname + ' and name' + lipid_atoms).split('residue')
               
        
      #  if lipid == 'CHOL':
        print(lipid + "  " +lipid_resname)
     #   print(lipid_atoms)
     #   print(lipid_atoms_selection)
     
        #every 500th frame for checking flip flops
        for ts in u.trajectory[start:end:skip]:
            # print('frame ' + str(ts.frame))
            time = u.trajectory.time
            
            R_m = membraneCentreOfMass(u, system)
            
            l_leaflet_fr = []
                        
            for l in lipid_atoms_selection:  
                l_z = l.center_of_mass()[2]
                if l_z < R_m:
                    l_leaflet_fr.append('l1') #in leaflet 1
                if l_z > R_m:
                    l_leaflet_fr.append('l2') #in leaflet 2
                               
                        
                 #check if lipid has changed leaflet between frames
            if ts.frame != 0:
                for i in range(0,len(l_leaflet_fr)):
                    if l_leaflet_fr[i] != previous_leaflet[i]:
                        print("FLIP FLOP")
                        flipflops += 1
                        frames_lipids.append([str(time - time_diff) + ' ' + str(time) + '   ' + lipid])
                                    
                        
                                    
                        
            previous_leaflet = l_leaflet_fr.copy()    
                            
       # print(flipflops)                
                       
                            
    if flipflops > 0:
        print('number of flip flops in trajectory is ' + str(flipflops) + " in the simulation in directory " + indexingPath)
       # flipflop_file = '../../Data/FlipFlop/flipflop.txt'     #+ indexingPath + '/flipflop.txt'
       # save frame nr and lipid name
        with open(flipflop_file, 'a') as f:
            for fr in frames_lipids:
                f.write(str(fr).strip("[]'") + '\n')
        f.close()
    else:
        print("no flip flops in trajectory")


                                    
                
                
