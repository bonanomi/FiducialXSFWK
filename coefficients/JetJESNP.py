from tqdm import tqdm
import optparse
import json
import ROOT
from ROOT import *
from ROOT import TFile, TH1, TH1F, TCanvas, gSystem, TRatioPlot, TPad, TStyle, TChain, gStyle

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

import numpy 

def parseOptions():

    global opt, args, runAllSteps

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
global opt, args, runAllSteps
parseOptions()

def computeJES(obsname, year, m4l_low, m4l_high):
    obs_bins = {0:(opt.OBSBINS.split("|")[1:(len(opt.OBSBINS.split("|"))-1)]),1:['0','inf']}[opt.OBSBINS=='inclusive']
    obs_bins = [float(i) for i in obs_bins] #Convert a list of str to a list of float
    edges = numpy.asarray(obs_bins)
    NBINS = int(len(edges) - 1)

    print(NBINS, edges)

    h = TH1F("h", "", NBINS, edges)
    h_up = TH1F("h_up", "", NBINS, edges)
    h_dn = TH1F("h_dn", "", NBINS, edges)

    for y in ['2016', '2017', '2018']:
        ch = TChain('candTree')
        for pmode in tqdm(pmodes):
            path = indir + y + '/' + pmode + '/' + pmode + '_reducedTree_MC_' + y + '.root'
            ch.Add(path)
            # if year != 'Full':
            #     path = indir + year + '/' + pmode + '/' + pmode + '_reducedTree_MC_' + year + '.root'
            #     ch.Add(path)
            # else:
            #     for y in ['2016', '2017', '2018']:
            #         path = indir + y + '/' + pmode + '/' + pmode + '_reducedTree_MC_' + y + '.root'
            #         ch.Add(path)

        toCut = "1"; toCut_dn = "1"; toCut_up = "1"
        #if 'njets' not in obsname:
        #    toCut = 'njets_pt30_eta2p5>0'
        #    toCut_dn = 'njets_pt30_eta2p5_jesdn>0'
        #    toCut_up = 'njets_pt30_eta2p5_jesup>0'

        ## Otherwise doesn't fill the histos... weird
        h = TH1F("h", "", NBINS, edges)
        h_up = TH1F("h_up", "", NBINS, edges)
        h_dn = TH1F("h_dn", "", NBINS, edges)

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
                print('%.2f/%.2f' %(up, dn))
            else:
                print('%.2f/%.2f' %(dn, up))
            np = '%.2f/%.2f' %(dn, up)
            jesNP['recobin%i' %i] = np
            i+=1

        with open('../inputs/JESNP_'+y+'_'+obsname+'.py', 'w') as f:
            f.write('obsbins = ' + str(obs_bins) + '\n') 
            f.write('JESNP = ' + str(jesNP) + '\n')

indir = '/eos/user/a/atarabin/MC_samples/'
pmodes = ['ggH125', 'VBFH125', 'ttH125', 'WminusH125', 'WplusH125', 'ZH125']

obsname = opt.OBSNAME
year    = opt.YEAR
m4l_low = opt.LOWER_BOUND
m4l_high = opt.UPPER_BOUND

print(obsname, year, m4l_low, m4l_high)

computeJES(obsname, year, m4l_low, m4l_high)
