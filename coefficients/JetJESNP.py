import numpy
from tqdm import tqdm
import optparse, os, sys
import json
import numpy      
import ROOT
# from ROOT import *
# from ROOT import TFile, TH1, TH1F, TCanvas, gSystem, TRatioPlot, TPad, TStyle, TChain, gStyle
from binning import binning
from createdf_jes import skim_df

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

sys.path.append('../inputs/')
from observables import observables
_temp = __import__('observables', globals(), locals(), ['observables'], -1)
observables = _temp.observables
sys.path.remove('../inputs/')

def parseOptions():

    global opt, args

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--m4lLower',  dest='LOWER_BOUND',  type='string',default='105.0',   help='Lower bound for m4l')
    parser.add_option('',   '--m4lUpper',  dest='UPPER_BOUND',  type='string',default='140.0',   help='Upper bound for m4l')
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    if (opt.OBSBINS=='' and opt.OBSNAME!='inclusive'):
        parser.error('Bin boundaries not specified for differential measurement. Exiting...')
        sys.exit()


# parse the arguments and options
global opt, args
parseOptions()

def computeJES(obsname, obs_bins, year, fState, m4l_low, m4l_high, doubleDiff, obsname_out):

	sel_m4l = (d_sig[year].ZZMass > m4l_low) and (d_sig[year].ZZMass < m4l_high)

	


    edges = numpy.asarray(obs_bins)
    NBINS = int(len(edges) - 1)

    print(NBINS, edges)

    h = TH1F("h", "", NBINS, edges)
    h_up = TH1F("h_up", "", NBINS, edges)
    h_dn = TH1F("h_dn", "", NBINS, edges)
    
    ch = TChain('candTree')
    for pmode in tqdm(pmodes):
        path = indir + year + '/' + pmode + '/' + pmode + '_reducedTree_MC_' + year + '.root'
        ch.Add(path)



    toCut = "1"; toCut_dn = "1"; toCut_up = "1"
    if 'njets' not in obsname:
        toCut = 'njets_pt30_eta2p5>0'
        toCut_dn = 'njets_pt30_eta2p5_jesdn>0'
        toCut_up = 'njets_pt30_eta2p5_jesup>0'

    ## Otherwise doesn't fill the histos... weird
    h = ROOT.TH1F("h", "", NBINS, edges)
    h_up = ROOT.TH1F("h_up", "", NBINS, edges)
    h_dn = ROOT.TH1F("h_dn", "", NBINS, edges)

    massRange = 'ZZMass>%s && ZZMass<%s' %(m4l_low, m4l_high)

    ch.Draw("%s>>h" %obsname, "(%s && %s)" %(massRange, toCut), "goff")
    ch.Draw("%s_jesup>>h_up" %obsname, "(%s && %s)" %(massRange, toCut_up), "goff")
    ch.Draw("%s_jesdn>>h_dn" %obsname, "(%s && %s)" %(massRange, toCut_dn), "goff")

    ## Normalize to unity, only interested in shapes comparison
    h.Scale(1./h.Integral())
    h_dn.Scale(1./h_dn.Integral())
    h_up.Scale(1./h_up.Integral())

    h_up.SetLineColor(2)
    h_dn.SetLineColor(3)

    h.SetLineWidth(2)
    h_up.SetLineWidth(2)
    h_dn.SetLineWidth(2)

    c = TCanvas()
    h.Draw("HIST")
    h_up.Draw("HIST SAME")
    h_dn.Draw("HIST SAME")
    c.Draw()
    c.SaveAs("JES.pdf")
    
    c1 = TCanvas()
    rp_dn = TRatioPlot(h, h_dn)
    rp_dn.Draw()
    c1.Draw()
    c1.SaveAs("JESDn.pdf")
    c2 = TCanvas()
    rp_up = TRatioPlot(h, h_up)
    rp_up.Draw()
    c2.Draw()
    c1.SaveAs("JESUp.pdf")
    graph_up = rp_up.GetLowerRefGraph()
    graph_dn = rp_dn.GetLowerRefGraph()

    jes_up = numpy.array(graph_up.GetY())
    jes_dn = numpy.array(graph_dn.GetY())
  
    jesNP = {}
    print(jes_up, jes_dn)
    i = 0
    for up, dn in zip(jes_up, jes_dn):
        if up > dn:
            print('%.3f/%.3f' %(up, dn))
        else:
            print('%.3f/%.3f' %(dn, up))
        np = '%.3f/%.3f' %(dn, up)
        jesNP['recobin%i' %i] = np
        i+=1

    with open('../inputs/JESNP_'+year+'_'+obsname_out+'.py', 'w') as f:
        f.write('obsbins = ' + str(obs_bins) + '\n') 
        f.write('JESNP = ' + str(jesNP) + '\n')
    if doubleDiff:
        return jesNP


## ---------------------------- Main ----------------------------
indir = '/eos/user/a/atarabin/MC_samples/'
# pmodes = ['ggH125', 'VBFH125', 'ttH125', 'WminusH125', 'WplusH125', 'ZH125']

obsname = opt.OBSNAME
obsname_out = obsname
doubleDiff = False

if ' vs ' in obsname:
    doubleDiff = True
    obsname_1st = opt.OBSNAME.split(' vs ')[0]
    obsname_2nd = opt.OBSNAME.split(' vs ')[1]
    obsname_out = obsname_1st + '_' + obsname_2nd

if doubleDiff:
    obs_reco_2nd = observables[obsname]['obs_reco_2nd']

obs_reco = observables[obsname]['obs_reco']

year    = opt.YEAR
if (opt.YEAR == '2016'): years = [2016]
if (opt.YEAR == '2017'): years = [2017]
if (opt.YEAR == '2018'): years = [2018]
if (opt.YEAR == 'Full'): years = [2016,2017,2018]

m4l_low = opt.LOWER_BOUND
m4l_high = opt.UPPER_BOUND

obs_bins, doubleDiff = binning(opt.OBSBINS)


# if not doubleDiff: #It is not a double-differential analysis
#     obs_bins = {0:(opt.OBSBINS.split("|")[1:(len(opt.OBSBINS.split("|"))-1)]),1:['0','inf']}[opt.OBSBINS=='inclusive']
#     obs_bins = [float(i) for i in obs_bins] #Convert a list of str to a list of float
#     print 'It is a single-differential measurement, binning', obs_bins
# else: #It is a double-differential analysis
#     # Implementing only the first one as it will be used only for jet-related measurements
#     if opt.OBSBINS.count('vs')==1 and opt.OBSBINS.count('/')>=1:
#         obs_bins_tmp = opt.OBSBINS.split(" vs ")
#         obs_bins_1st = obs_bins_tmp[0].split('|')[1:len(obs_bins_tmp[0].split('|'))-1] #['0', '1', '2', '3', '20']
#         obs_bins_1st = [float(i) for i in obs_bins_1st] #Convert a list of str to a list of float
#         obs_bins_tmp = obs_bins_tmp[1].split(' / ') #['|0|10|20|45|90|250|', '|0|10|20|80|250|', '|0|20|90|250|', '|0|25|250|']
#         obs_bins_2nd = {}
#         for i in range(len(obs_bins_tmp)): #At the end of the loop -> obs_bins_2nd {0: ['0', '10', '20', '45', '90', '250'], 1: ['0', '10', '20', '80', '250'], 2: ['0', '20', '90', '250'], 3: ['0', '25', '250']}
#             obs_bins_2nd[i] = obs_bins_tmp[i].split('|')[1:len(obs_bins_tmp[i].split('|'))-1]
#             obs_bins_2nd[i] = [float(j) for j in obs_bins_2nd[i]] #Convert a list of str to a list of float
#         obs_bins = obs_bins_2nd[1]
      

print obsname, year, m4l_low, m4l_high, obsname_out

# Generate dataframes
d_sig = {}
for year in years:
    if doubleDiff: sig = skim_df(year, doubleDiff, obs_reco, obs_reco_2nd)
    else: sig = skim_df(year, doubleDiff, obs_reco)
    d_sig[year] = sig
	d_sig[year] = pd.concat([d_sig[year]['ggH125'], d_sig[year]['VBFH125'], d_sig[year]['WH125'], d_sig[year]['ZH125'], d_sig[year]['ttH125']])

print(d_sig)


# if not doubleDiff:
#     for year in years:
#     	for fState in ['2e2mu', '4mu' ,'4e']:
# 	        # year = str(year)
# 	        computeJES(obs_reco, obs_bins, year, fState, m4l_low, m4l_high, doubleDiff, obsname_out)
# else:
#     for year in years:
#         year = str(year)
#         jesNP0 = computeJES(obs_reco, obs_bins_1st, year, m4l_low, m4l_high, doubleDiff, obs_reco)
#         jesNP1 = computeJES(obs_reco_2nd, obs_bins, year, m4l_low, m4l_high, doubleDiff, obsname_out)
#         jesNP = {}
#         jesNP['recobin0'] = jesNP0['recobin0']
#         for i in range(len(obs_bins)-1):
#             jesNP['recobin%i' %(i+1)] = jesNP1['recobin%i' %i]
#         with open('../inputs/JESNP_'+year+'_'+obsname_out+'.py', 'w') as f:
#             f.write('obsbins = ' + str(obs_bins) + '\n')
#             f.write('JESNP = ' + str(jesNP) + '\n')
