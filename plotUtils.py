import os, sys
sys.path.append('./inputs')
from higgs_xsbr_13TeV import *

def getObsName(obsName):
    doubleDiff = False
    if ' vs ' in obsName: doubleDiff = True

    if doubleDiff:
        obs_name_1st = obsName.split(' vs ')[0]
        obs_name_2nd = obsName.split(' vs ')[1]
        obs_name = obs_name_1st + '_' + obs_name_2nd
    else:
        obs_name = obsName

    return obs_name, doubleDiff

def parseBins(OBSBINS):
    if not 'vs' in OBSBINS: #It is not a double-differential analysis
        obs_bins = {0:(OBSBINS.split("|")[1:(len(OBSBINS.split("|"))-1)]),1:['0','inf']}[OBSBINS=='inclusive']
        obs_bins = [float(i) for i in obs_bins] #Convert a list of str to a list of float
        doubleDiff = False
        print('It is a single-differential measurement, binning', obs_bins)
    else: #It is a double-differential analysis
        doubleDiff = True
        # The structure of obs_bins is:
        # index of the dictionary is the number of the bin
        # [obs_bins_low, obs_bins_high, obs_bins_low_2nd, obs_bins_high_2nd]
        # The first two entries are the lower and upper bound of the first variable
        # The second two entries are the lower and upper bound of the second variable
        if OBSBINS.count('vs')==1 and OBSBINS.count('/')>=1: #Situation like this one '|0|1|2|3|20| vs |0|10|20|45|90|250| / |0|10|20|80|250| / |0|20|90|250| / |0|25|250|'
            obs_bins_tmp = OBSBINS.split(" vs ") #['|0|1|2|3|20|', '|0|10|20|45|90|250| / |0|10|20|80|250| / |0|20|90|250| / |0|25|250|']
            obs_bins_1st = obs_bins_tmp[0].split('|')[1:len(obs_bins_tmp[0].split('|'))-1] #['0', '1', '2', '3', '20']
            obs_bins_1st = [float(i) for i in obs_bins_1st] #Convert a list of str to a list of float
            obs_bins_tmp = obs_bins_tmp[1].split(' / ') #['|0|10|20|45|90|250|', '|0|10|20|80|250|', '|0|20|90|250|', '|0|25|250|']
            obs_bins_2nd = {}
            for i in range(len(obs_bins_tmp)): #At the end of the loop -> obs_bins_2nd {0: ['0', '10', '20', '45', '90', '250'], 1: ['0', '10', '20', '80', '250'], 2: ['0', '20', '90', '250'], 3: ['0', '25', '250']}
                obs_bins_2nd[i] = obs_bins_tmp[i].split('|')[1:len(obs_bins_tmp[i].split('|'))-1]
                obs_bins_2nd[i] = [float(j) for j in obs_bins_2nd[i]] #Convert a list of str to a list of float
            obs_bins = {}
            k = 0 #Bin index
            for i in range(len(obs_bins_1st)-1):
                for j in range(len(obs_bins_2nd[i])-1):
                    obs_bins[k] = []
                    obs_bins[k].append(obs_bins_1st[i])
                    obs_bins[k].append(obs_bins_1st[i+1])
                    obs_bins[k].append(obs_bins_2nd[i][j])
                    obs_bins[k].append(obs_bins_2nd[i][j+1])
                    k +=1
        elif OBSBINS.count('vs')>1 and OBSBINS.count('/')>1: #Situation like this one '|50|80| vs |10|30| / |50|80| vs |30|60| / |80|110| vs |10|25| / |80|110| vs |25|30|'
            obs_bins_tmp = OBSBINS.split(' / ') #['|50|80| vs |10|30|', '|50|80| vs |30|60|', '|80|110| vs |10|25|', '|80|110| vs |25|30|']
            obs_bins_1st={}
            obs_bins_2nd={}
            obs_bins={}
            for i in range(len(obs_bins_tmp)): #At the end of the loop -> obs_bins_1st {0: ['50', '80'], 1: ['50', '80'], 2: ['80', '110'], 3: ['80', '110']} and obs_bins_2nd {0: ['10', '30'], 1: ['30', '60'], 2: ['10', '25'], 3: ['25', '30']}
                obs_bins_tmp_bis = obs_bins_tmp[i].split(' vs ')
                obs_bins_1st[i] = obs_bins_tmp_bis[0].split('|')[1:len(obs_bins_tmp_bis[0].split('|'))-1]
                obs_bins_1st[i] = [float(j) for j in obs_bins_1st[i]] #Convert a list of str to a list of float
                obs_bins_2nd[i] = obs_bins_tmp_bis[1].split('|')[1:len(obs_bins_tmp_bis[1].split('|'))-1]
                obs_bins_2nd[i] = [float(j) for j in obs_bins_2nd[i]] #Convert a list of str to a list of float
                obs_bins[i] = []
                obs_bins[i].append(obs_bins_1st[i][0])
                obs_bins[i].append(obs_bins_1st[i][1])
                obs_bins[i].append(obs_bins_2nd[i][0])
                obs_bins[i].append(obs_bins_2nd[i][1])
        elif OBSBINS.count('vs')==1 and OBSBINS.count('/')==0: #Situation like this one '|0|1|2|3|20| vs |0|10|20|45|90|250|'
            obs_bins_tmp = OBSBINS.split(" vs ") #['|0|1|2|3|20|', '|0|10|20|45|90|250|']
            obs_bins_1st = obs_bins_tmp[0].split('|')[1:len(obs_bins_tmp[0].split('|'))-1] #['0', '1', '2', '3', '20']
            obs_bins_1st = [float(i) for i in obs_bins_1st] #Convert a list of str to a list of float
            obs_bins_2nd = obs_bins_tmp[1].split('|')[1:len(obs_bins_tmp[1].split('|'))-1] #['0', '10', '20', '45', '90', '250']
            obs_bins_2nd = [float(i) for i in obs_bins_2nd] #Convert a list of str to a list of float
            obs_bins = {}
            k = 0 #Bin index
            for i in range(len(obs_bins_1st)-1):
                for j in range(len(obs_bins_2nd)-1):
                    obs_bins[k] = []
                    obs_bins[k].append(obs_bins_1st[i])
                    obs_bins[k].append(obs_bins_1st[i+1])
                    obs_bins[k].append(obs_bins_2nd[j])
                    obs_bins[k].append(obs_bins_2nd[j+1])
                    k +=1
        else:
            print('Problem in the definition of the binning')
            quit()
        print('It is a double-differential measurement, binning for the 1st variable', obs_bins_1st, 'and for the 2nd variable', obs_bins_2nd)
        print(obs_bins)

    return obs_bins

## For 2D measurements
def binBoundaries(obs_bins):
    range_1st = []
    range_2nd = []
    for _bin in obs_bins:
        range_1st.append([obs_bins[_bin][0], obs_bins[_bin][1]])
        range_2nd.append([obs_bins[_bin][2], obs_bins[_bin][3]])
    return range_1st, range_2nd

def getMath(obsName):
    if(obsName == 'mass4l'): label = 'm_{4l}'
    elif(obsName == 'rapidity4l'): label = '|y_{H}|'
    elif(obsName == 'pT4l'): label = 'p_{T}^{H} \\mathrm{ (GeV)}'
    elif(obsName == 'massZ1'): label = 'm_{Z1} \\mathrm{ (GeV)}'
    elif(obsName == 'massZ2'): label = 'm_{Z2} \\mathrm{ (GeV)}'
    elif(obsName == 'njets_pt30_eta2p5'): label = '\\mathrm{nJets} (2.5)' #label = 'nJets, pT>30 GeV, |#eta|<2.5'
    elif(obsName == 'njets_pt30_eta4p7'): label = '\\mathrm{nJets} (4.7)' 
    elif(obsName == 'pTj1'): label = 'p_{T}^{(j1, 2.5)} \\mathrm{ (GeV)}'
    elif(obsName == 'pTj1_eta4p7'): label = 'p_{T}^{(j1, 4.7)} \\mathrm{ (GeV)}'
    elif(obsName == 'pTj2'): label = 'p_{T}^{(Sub. jet)} \\mathrm{ (GeV)}'
    elif(obsName == 'pTHj'): label = 'p_{T}^{(Hj)} \\mathrm{ (GeV)}'
    elif(obsName == 'pTHjj'): label = 'p_{T}^{(Hjj)} \\mathrm{ (GeV)}'
    elif(obsName == 'mHj'): label = 'm(Hj) \\mathrm{ (GeV)}'
    elif(obsName == 'mHjj'): label = 'm(Hjj) \\mathrm{ (GeV)}'
    elif(obsName == 'mjj'): label = 'm(jj) \\mathrm{ (GeV)}'
    elif(obsName == 'mass4l'): label = 'm_{4\\ell} \\mathrm{ (GeV)}'
    elif(obsName == 'costhetastar'): label = 'cos(\\theta^{*})'
    elif(obsName == 'costhetaZ1'): label = 'cos(\\theta_{1})'
    elif(obsName == 'costhetaZ2'): label = 'cos(\\theta_{2})'
    elif(obsName == 'phistar'): label = '\\Phi_{1}'
    elif(obsName == 'phi'): label = '\\Phi'
    elif((obsName == 'TCj') or (obsName == 'TCjmax')): label = '\\mathcal{T}_{Cj} \\mathrm{ (GeV)}'
    elif((obsName == 'TBj') or (obsName == 'TBjmax')): label = '\\mathcal{T}_{Bj} \\mathrm{ (GeV)}'
    else: label = obsName
    return label

def computeXSTH(acc, acc_NNLOPS, theoryMass, obsName, obs_bins, doubleDiff):
    # cross sections
    ggH_powheg = []
    ggH_minloHJ = []
    # XH unc
    XH = []

    nBins=len(obs_bins)
    if not doubleDiff: nBins = nBins -1    
    for obsBin in range(nBins):

        # theory cross sections
        ggH_powheg.append(0.0)
        ggH_minloHJ.append(0.0)
        # XH
        XH.append(0.0)

        for channel in ['4e','4mu','2e2mu']:
            XH_fs = higgs_xs['VBF_'+theoryMass]*higgs4l_br[theoryMass+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['WH_'+theoryMass]*higgs4l_br[theoryMass+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['ZH_'+theoryMass]*higgs4l_br[theoryMass+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['ttH_'+theoryMass]*higgs4l_br[theoryMass+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            
            XH[obsBin]+=XH_fs
            ggH_xsBR = higgs_xs['ggH_'+theoryMass]*higgs4l_br[theoryMass+'_'+channel]
            
            ggH_powheg[obsBin]+=ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            ggH_minloHJ[obsBin]+=ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]

    return ggH_powheg, ggH_minloHJ, XH

## TODO ##
# def computeUNCTH(acc, acc_NNLOPS, theoryMass, obsName, obs_bins):
