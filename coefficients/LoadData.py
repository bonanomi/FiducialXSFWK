from ROOT import *
from array import array
import os


# 80X samples

dirMC_94 = '/eos/user/a/atarabin/MC_samples/2017/merged'
print "samples directory: ", dirMC_94

SamplesMC_94 = [
'ggH125_mergedTree_MC_2017.root',
'ggH125_NNLOPS_mergedTree_MC_2017',
'VBFH125_mergedTree_MC_2017.root',
'ZH125_mergedTree_MC_2017.root',
'ttH125_mergedTree_MC_2017.root',
'WH125_mergedTree_MC_2017.root'
]

SamplesData_94 = [
'DoubleEG_Run2017B-17Nov2017-v1.root','DoubleEG_Run2017C-17Nov2017-v1.root','DoubleEG_Run2017D-17Nov2017-v1.root','DoubleEG_Run2017E-17Nov2017-v1.root','DoubleEG_Run2017F-17Nov2017-v1.root',
'DoubleMuon_Run2017-17Nov2017-v1.root','DoubleMuon_Run2017B-17Nov2017-v1.root','DoubleMuon_Run2017C-17Nov2017-v1.root'
]
###################################################### 
RootFile = {} 
Tree = {} 
nEvents = {} 
sumw = {}


# 80X MC
for i in range(0,len(SamplesMC_94)):

    sample = SamplesMC_94[i].rstrip('.root')

    RootFile[sample] = TFile(dirMC_94+'/'+sample+'.root',"READ")
    Tree[sample] = RootFile[sample].Get("fullTree")
    
    h_nevents = RootFile[sample].Get("Counters")
    h_sumw = RootFile[sample].Get("Counters")

    if (h_nevents): nEvents[sample] = h_nevents.GetBinContent(40) #not really nevents. Will work for now
    else: nEvents[sample] = 0.

    if (h_sumw): sumw[sample] = h_sumw.GetBinContent(40)
    else: sumw[sample] = 0.

    if (not Tree[sample]): print sample+' has no passedEvents tree'
    else:
        print sample,"nevents",nEvents[sample],"sumw",sumw[sample]

