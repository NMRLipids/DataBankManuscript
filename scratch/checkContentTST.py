def getTotalNlipids(system):
    NLIPIDS = 0
    for molecule in system['COMPOSITION']:
        if molecule in lipids_dict:
            NLIPIDS += np.sum(system['COMPOSITION'][molecule]['COUNT'])
    return NLIPIDS

def getTotalNsolvent(system):
    NMOLECULES = 0
    for molecule in system['COMPOSITION']:
        if molecule not in lipids_dict:
            NMOLECULES += np.sum(system['COMPOSITION'][molecule]['COUNT'])
    return NMOLECULES

def checkAvailabilitySIM(system,lipid,counterion):
    status = {}
    TotalNlipids = getTotalNlipids(system)
    TotalNsolvent = getTotalNsolvent(system)
    Nwater = system['COMPOSITION']['SOL']['COUNT']
    if counterion != 'no':
        try:
            Ncounterion = system['COMPOSITION'][counterion]['COUNT']
        except:
            Ncounterion = 0
    
    try:
        Nlipid = np.sum(system['COMPOSITION'][lipid]['COUNT'])
    except:
        Nlipid = 0
    
    path = system['path']
    
    QualityEvaluated = False
    TotalQualityFilePath = path + '/SYSTEM_quality.json'
    if (os.path.isfile(TotalQualityFilePath)):
        with open(TotalQualityFilePath) as json_file:
            Quality = json.load(json_file)
        json_file.close()
        if all(value > 0 for value in Quality.values()):
            #print(Quality)
            QualityEvaluated = True
            
    xrayQualityEvaluated = False
    xrayQualityFilePath = path + '/FormFactorQuality.json'
    if (os.path.isfile(xrayQualityFilePath)):
        with open(xrayQualityFilePath) as json_file:
            xrayQuality = json.load(json_file)
        json_file.close()
        if len(xrayQuality) > 0 and xrayQuality[0] > 0:
            xrayQualityEvaluated = True
    
    SingleComponentSystem = False
    if Nlipid == TotalNlipids:
        if counterion == 'no' and Nwater == TotalNsolvent:
            SingleComponentSystem = True
        if counterion != 'no'and  Nwater == TotalNsolvent-Ncounterion and Nlipid == Ncounterion:    
            SingleComponentSystem = True
        
    if SingleComponentSystem:
        status['Simulation'] = 'yes'
        status['FF'] = system['FF']
        if QualityEvaluated:
            status['Experiment'] = 'yes'
            status['Quality'] = Quality['total']
            
        if not QualityEvaluated:
            status['Experiment'] = 'no'
            status['Quality'] = 0 
            
        if xrayQualityEvaluated:
            status['xrayExperiment'] = 'yes'
            status['xrayQuality'] = xrayQuality[0]
            
        if not xrayQualityEvaluated:
            status['xrayExperiment'] = 'no'
            status['xrayQuality'] = 0 
                
    return status 
        
        
def giveStatus(systems,lipid,counterion):  
    status = {'Simulation': 'no', 'Experiment': 'no', 'Quality': 0, 'xrayExperiment': 'no', 'xrayQuality': 0}
    QualityEvaluatedFound = False
    for system in systems:
        TMPstatus = checkAvailabilitySIM(system,lipid,counterion)
        if TMPstatus and TMPstatus['Quality'] > status['Quality']:
            QualityEvaluatedFound = True
            status = TMPstatus 
        if TMPstatus and not QualityEvaluatedFound and TMPstatus['Simulation'] == 'yes':
            status = TMPstatus
    return status
    
def giveExpStatus(lipids,counterion,status):

    

    return status


    
HGs = {'PC', 'PE', 'PG', 'PS'}
tails = {'PO', 'DO', 'DP'}
table = {}

for tail in tails:
    table[tail] = {}
    for HG in HGs:
        lipid = tail + HG
        
        if lipid == 'POPS' or lipid == 'POPG' or lipid == 'DPPG':
            counterion = 'SOD'
        else:
            counterion = 'no'
        
        status = giveStatus(systems,lipid,counterion)
        statusString = ''
            
        if status['Quality'] > 0:
            if 'ECC-lipids' in status['FF']:
                FF = 'ECClipids'
            else:
                FF = status['FF']
            #print(status['xrayQuality'])
            statusString = FF + '(' + str(round(status['Quality'],2)) + ',' + str(round(status['xrayQuality'])) + ')'
        else:
            if status['Simulation'] == 'no':
                statusString = statusString + 'MD,'
            if status['Experiment'] == 'no':
                statusString = statusString + 'NMR,'
            if status['xrayExperiment'] == 'no':
                statusString = statusString + 'x-ray'
            
        table[tail][HG] = statusString


print(pd.DataFrame(table))
display(pd.DataFrame(table))
