#Lipid flip-flops
#loops over trajectories and checks if some molecule (lipid or cholesterol) in membrane flips

import os
import sys
import yaml
import json
import math
import pdb
from random import randint
import urllib.request
from urllib.error import URLError,HTTPError,ContentTooShortError

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from scipy.stats import norm
import MDAnalysis as mda
from lipyphilic.lib.assign_leaflets import AssignLeaflets
from lipyphilic.lib.flip_flop import FlipFlop 

sys.path.insert(1, '../../Databank/Scripts/BuildDatabank/')
from databankLibrary import download_link, lipids_dict, databank



def loadMappingFile(path_to_mapping_file):    # load mapping file into a dictionary
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

def headgroupAtom(readme, lipid):
    path_to_mapping_file = readme['COMPOSITION'][lipid]['MAPPING']
    mapping_dict = loadMappingFile(path_to_mapping_file)
    for key in mapping_dict:
        if mapping_dict[key]['FRAGMENT'] == 'headgroup': 
            return mapping_dict[key]['ATOMNAME']

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
    
def centerxtcfile(system, system_path):
    # FIND LAST CARBON OF SN-1 TAIL AND G3 CARBON
    for molecule in system['COMPOSITION']:
        if molecule in lipids_dict:
            mapping_file = f"../../Databank/Scripts/BuildDatabank/mapping_files/{system['COMPOSITION'][molecule]['MAPPING']}"
            with open(mapping_file, "r") as yaml_file:
                mapping = yaml.load(yaml_file,  Loader=yaml.FullLoader)
            try:
                G3atom = mapping['M_G3_M']['ATOMNAME']
            except:
                pass

            for Cindex in range(1,30):
                atom = f'M_G1C{str(Cindex)}_M'
                try:
                    lastAtom = mapping[atom]['ATOMNAME']
                except:
                    continue

    # Center around one lipid tail CH3 to guarantee all lipids in the same box
    if 'gromacs' in system['SOFTWARE']:
        if 'WARNINGS' in system and 'GROMACS_VERSION' in system['WARNINGS'] and system['WARNINGS']['GROMACS_VERSION'] == 'gromacs3':
            trjconvCOMMAND = '/home/osollila/Programs/gromacs/gromacs402/bin/trjconv'
            makendxCOMMAND = '/home/osollila/Programs/gromacs/gromacs402/bin/make_ndx'
        else:
            trjconvCOMMAND = 'gmx trjconv'
            makendxCOMMAND = 'gmx make_ndx'

        os.system('rm foo.ndx')
        os.system(f'echo "a {lastAtom}\nq" | {makendxCOMMAND} -f {tpr_name} -o foo.ndx')
        os.system("tail -n1 foo.ndx | awk '{print $NF}'")
        os.system('echo "[ centralAtom ]" >> foo.ndx')
        os.system("tail -n2 foo.ndx | head -n1 |  awk '{print $NF}' >> foo.ndx")

        xtcwhole = f'{system_path}/whole.xtc'
        xtcfoo = f'{system_path}/foo2.xtc'
        xtccentered = f'{system_path}/centered.xtc'
        if (not os.path.isfile(xtccentered)):
            print("Make molecules whole in the trajectory")
            if (not os.path.isfile(xtcfoo)):
                os.system(f'echo "centralAtom\nSystem" |  {trjconvCOMMAND} -center -pbc mol -n foo.ndx -f {xtcwhole} -s {tpr_name} -o {xtcfoo}')

            os.system('rm foo.ndx')
            os.system(f'rm {xtcwhole}')

            os.system(f'echo "a {G3atom}\nq" | {makendxCOMMAND} -f {tpr_name} -o foo.ndx')
            os.system(f'echo "{G3atom}\nSystem" |  {trjconvCOMMAND} -center -pbc mol -n foo.ndx -f {xtcfoo} -s {tpr_name} -o {xtccentered}')
            os.system(f'rm {xtcfoo}')
    else:
        print('Centering for other than Gromacs may not work if there are jumps over periodic boundary conditions in z-direction.')


print("Loading README files")
path = '../../Databank/Data/Simulations/'
db_data = databank(path)      #searches through every subfolder of a path and finds every trajectory in databank
systems = db_data.get_systems()

for system in systems:   #checks everysystem for flipflop
    #directories
    indexingPath = "/".join(system['path'].split("/")[5:9])
    subdir = f'../../Databank/Data/Simulations/{indexingPath}/'
    DATAdir = f'../Data/Flipflops/{indexingPath}/'

    if 'WARNINGS' in system and 'GROMACS_VERSION' in system['WARNINGS'] and system['WARNINGS']['GROMACS_VERSION'] == 'gromacs3':
        continue
        trjconvCOMMAND = '/home/osollila/Programs/gromacs/gromacs402/bin/trjconv'
    else:
        trjconvCOMMAND = 'gmx trjconv'

    #exception skiping
    if system['SOFTWARE'] == 'openMM':
        print('OpenMM')
        continue
        print('Warning: using openMM which cant be centered by this code')

    if system['SOFTWARE'] == 'gromacs':                 
        if 'SOFTWARE_VERSION' in system:
            if system['SOFTWARE_VERSION'] == '3.x':
                continue                                          #comment this if you have gromacs 3
                pass

    if system['TYPEOFSYSTEM'] == 'miscellaneous':
        continue

    if os.path.isfile(flipflop_file):
        print("skipped")
        continue

#    if indexingPath.split('/')[2] != '304765469045020262110ceea29c4a404b503c63':    #testing trajectory
#        continue

    #building directory and file system
    if (not os.path.isdir(DATAdir)):
        hashes = indexingPath.split('/')
        os.system(f"mkdir ../Data/Flipflops/{hashes[0]}")
        os.system(f"mkdir ../Data/Flipflops/{hashes[0]}/{hashes[1]}")
        os.system(f"mkdir ../Data/Flipflops/{hashes[0]}/{hashes[1]}/{hashes[2]}")
        os.system(f"mkdir ../Data/Flipflops/{hashes[0]}/{hashes[1]}/{hashes[2]}/{hashes[3]}") 

    os.system(f'cp {subdir}README.yaml {DATAdir}/')
   
    #getting data from databank and proprocessing them                
    doi = system['DOI']
    trj = system['TRJ'][0][0]
    tpr = system['TPR'][0][0]
    trj_name = subdir + trj
    tpr_name = subdir + tpr
    trj_url = download_link(doi, trj)
    tpr_url = download_link(doi, tpr)
    EQtime = float(system['TIMELEFTOUT'])*1000
    flipflop_file = f'{DATAdir}flipflop.dat'

    #downloading missing files
    try:
        if (not os.path.isfile(tpr_name)):
            #continue
            print("downloading tpr", doi)
            response = urllib.request.urlretrieve(tpr_url, tpr_name)
        
        if (not os.path.isfile(trj_name)):
            print("downloading trj", doi)
            response = urllib.request.urlretrieve(trj_url, trj_name)
    except:
        os.system(f"rm {subdir}/{trj_name} {subdir}/{tpr_name}")
        continue
                   
    #preparing trajetory
    xtcwhole = f'{subdir}whole.xtc'
    xtccentered = f'{subdir}centered.xtc'
    if not (os.path.isfile(xtcwhole) or os.path.isfile(xtccentered)):
        print("making molecules whole")
        os.system(f'echo System | {trjconvCOMMAND} -f {trj_name} -s {tpr_name} -o {xtcwhole} -pbc mol -b {EQtime}')
        print('sys1')

    if (not os.path.isfile(xtccentered)):
        print("centering")
        centerxtcfile(system,subdir)
    
    #setting MDAnalysis
    try:
        print("Loading system", doi)
        u = mda.Universe(tpr_name,xtccentered)
    except:
        print('Reading of data with MDanalysis failed', indexingPath)
        #pdb.set_trace()
        continue
        
    #reading lipid names and their head-group atoms
    lipids = getLipids(system)
    names = [system['COMPOSITION'][lipid]['NAME'] for lipid in lipids]
    decode = dict(zip(names, lipids))
    headgroupatoms = [headgroupAtom(system, lipid)for lipid in lipids]

    #atom selection in MDAnalysis
    HGnames = f"name {' '.join(headgroupatoms)}"
    lipidnames = f"resname {' '.join(lipids)}"

    lipid_sys = u.select_atoms(lipidnames)

    headatoms_sys = u.select_atoms(HGnames)
    
    #puting trajectory in to a much larger box to remove wraping of molecules
    print("Preprocessing trajectory")
    preparedfile = f'{subdir}smth.xtc'
    with mda.Writer(preparedfile, u.atoms.n_atoms) as w:
        for ts in u.trajectory:
            ts.dimensions = [1000, 1000, 1000, 90, 90, 90]
            lipid_center = headatoms_sys.center_of_geometry(pbc=True)
            dim = ts.triclinic_dimensions
            box_center = np.sum(dim, axis=0) / 2
            u.atoms.translate(box_center - lipid_center)
            w.write(u.atoms)
    u = mda.Universe(tpr_name,preparedfile)

    #analysing trajectory with LiPyphilic
    print("Assigning atoms to leaflets")
    leaflets = AssignLeaflets(
        universe=u,
        lipid_sel=HGnames,          # names of a headgroup atom from each lipid used in simulation
        midplane_sel=HGnames,       # only cholesterol is allowed to flip-flop      ??
        midplane_cutoff=10.0,          # buffer size for assigning molecules to the midplane
        )
    leaflets.run()
    
    print("Counting flipflops")
    flip_flop = FlipFlop(
        universe=u,
        lipid_sel=HGnames,
        leaflets=leaflets.filter_leaflets(HGnames),  # pass only the relevant leaflet data
        frame_cutoff=100,
        )

    flip_flop.run(
        start=None,
        stop=None,
        step=None
        )

    #writing output file
    with open(flipflop_file, 'w') as f:
        f.write(f'# res.; time-begin; time-end; to leaflet; outcome\n')
        f.close()               
                       
    if len(flip_flop.flip_flop_success) > 0:
        output = np.vstack((flip_flop.flip_flops.T,np.transpose(flip_flop.flip_flop_success))).T
        print(f'number of flip flops in trajectory is {len(output)} in the simulation in directory {indexingPath}')

        with open(flipflop_file, 'a') as f:
            for fr in output:
                fr[0] = str(int(fr[0])+1)
                if fr[1] >= fr[2]:
                    continue
                resname = decode[u.select_atoms(f'resid {fr[0]}').residues.resnames[0]]
                print(resname)
                f.write(f'{str(resname)}    {"    ".join(fr)}\n')
            f.close()
    else:
        print("no flip flops in trajectory")
    
    os.system(f'rm {preparedfile}')
