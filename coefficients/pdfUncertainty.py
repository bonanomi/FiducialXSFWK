import sys, os, string, re, pwd, commands, ast, optparse, shlex, time
from array import array
from math import *
from decimal import *
from sample_shortnames import *

grootargs = []
def callback_rootargs(option, opt, value, parser):
    grootargs.append(opt)
    
### Define function for parsing options
def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)
    
    # input options
    parser.add_option('-d', '--dir',    dest='SOURCEDIR',  type='string',default='./', help='run from the SOURCEDIR as working area, skip if SOURCEDIR is an empty string')
    parser.add_option('',   '--modelName',dest='MODELNAME',type='string',default='SM', help='Name of the Higgs production or spin-parity model, default is "SM", supported: "SM", "ggH", "VBF", "WH", "ZH", "ttH", "exotic","all"')
    parser.add_option('',   '--obsName',dest='OBSNAME',    type='string',default='',   help='Name of the observalbe, supported: "ZZMass", "pT4l", "massZ2", "rapidity4l", "cosThetaStar", "nets_reco_pt30_eta4p7"')
    parser.add_option('',   '--obsBins',dest='OBSBINS',    type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('-f', '--doFit', action="store_true", dest='DOFIT', default=False, help='doFit, default false')
    parser.add_option('-p', '--doPlots', action="store_true", dest='DOPLOTS', default=False, help='doPlots, default false')
    parser.add_option("-l",action="callback",callback=callback_rootargs)
    parser.add_option("-q",action="callback",callback=callback_rootargs)
    parser.add_option("-b",action="callback",callback=callback_rootargs)
                       
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

# parse the arguments and options
global opt, args, runAllSteps
parseOptions()
sys.argv = grootargs

doFit = opt.DOFIT
doPlots = opt.DOPLOTS

if (not os.path.exists("plots") and doPlots):
    os.system("mkdir plots")

from ROOT import *
from LoadData import *

RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)    

if (opt.DOPLOTS and os.path.isfile('tdrStyle.py')):
    from tdrStyle import setTDRStyle
    setTDRStyle()

Histos = {}
acceptance = {}
qcdUncert = {}
pdfUncert = {}

nnloWeights = {
    0: 'LHEweight_QCDscale_muR1_muF1',
    1: 'LHEweight_QCDscale_muR1_muF2',
    2: 'LHEweight_QCDscale_muR1_muF0p5',
    3: 'LHEweight_QCDscale_muR2_muF1',
    4: 'LHEweight_QCDscale_muR2_muF2',
    5: 'LHEweight_QCDscale_muR2_muF0p5',
    6: 'LHEweight_QCDscale_muR0p5_muF1',
    7: 'LHEweight_QCDscale_muR0p5_muF2',
    8: 'LHEweight_QCDscale_muR0p5_muF0p5'
}

def getunc(channel, List, m4l_bins, m4l_low, m4l_high, obs_reco, obs_gen, obs_bins, genbin):    

    obs_gen_low = obs_bins[genbin]
    obs_gen_high = obs_bins[genbin+1]

    obs_gen_lowest = obs_bins[0]
    obs_gen_highest = obs_bins[len(obs_bins)-1]
    
    if (obs_reco.startswith("ZZMass")):
        m4l_low = float(obs_gen_low)
        m4l_high = float(obs_gen_high)
        m4l_bins = int((m4l_high-m4l_low)/2)

    i_sample = -1

    print(List)

    for Sample in List:
        if (not Sample in Tree): continue
        if (not Tree[Sample]): continue
	print(Tree[Sample].GetName())
        if (obs_reco.startswith("njets")):
            cutobs_gen = "("+obs_gen+">="+str(obs_gen_low)+")"
        else:
            cutobs_gen = "("+obs_gen+">="+str(obs_gen_low)+" && "+obs_gen+"<"+str(obs_gen_high)+")"
        cutm4l_gen     = "(GenHMass>"+str(m4l_low)+" && GenHMass<"+str(m4l_high)+")"
        
        if (channel == "4l"):
            #cutchan_gen      = "((abs(GENlep_id[GENlep_Hindex[0]])==11 || abs(GENlep_id[GENlep_Hindex[0]])==13) && (abs(GENlep_id[GENlep_Hindex[2]])==11 || abs(GENlep_id[GENlep_Hindex[2]])==13) )"
            cutchan_gen      = "(abs(GenLep1Id)!= 15 && abs(GenLep2Id)!= 15 && abs(GenLep3Id)!= 15 && abs(GenLep4Id)!= 15)"
            cutchan_gen_out  = "((abs(GenZ1Flav)==121 || abs(GenZ1Flav)==169) && (abs(GenZ2Flav)==121 || abs(GenZ2Flav)==169))"
            cutm4l_gen       = "(GenHMass>"+str(m4l_low)+" && GenHMass<"+str(m4l_high)+")"
        if (channel == "4e"):
            cutchan_gen      = "(abs(GenLep1Id)==11 && abs(GenLep2Id)==11 && abs(GenLep3Id)==11 && abs(GenLep4Id)==11)"
            cutchan_gen_out  = "(abs(GenZ1Flav)==121 && abs(GenZ2Flav)==121)"
            cutm4l_gen       = "(GenHMass>"+str(m4l_low)+" && GenHMass<"+str(m4l_high)+")"
        if (channel == "4mu"):
            cutchan_gen      = "(abs(GenLep1Id)==13 && abs(GenLep2Id)==13 && abs(GenLep3Id)==13 && abs(GenLep4Id)==13)"
            cutchan_gen_out  = "(abs(GenZ1Flav)==169 && abs(GenZ2Flav)==169)"
            cutm4l_gen       = "(GenHMass>"+str(m4l_low)+" && GenHMass<"+str(m4l_high)+")"
        if (channel == "2e2mu"):
            #cutchan_gen      = "((abs(GenLep1Id)==11 && abs(GenLep2Id)==13 && abs(GenLep3Id)==11 && abs(GenLep4Id)==13) || (abs(GenLep1Id)==13 && abs(GenLep2Id)==11 && abs(GenLep3Id)==13 && abs(GenLep4Id)==11) || (abs(GenLep1Id)==11 && abs(GenLep2Id)==13 && abs(GenLep3Id)==13 && abs(GenLep4Id)==11) || (abs(GenLep1Id)==13 && abs(GenLep2Id)==11 && abs(GenLep3Id)==11 && abs(GenLep4Id)==13))"
            cutchan_gen      = "(abs(GenLep1Id)!=15) && !((abs(GenLep1Id)==13 && abs(GenLep2Id)==13 && abs(GenLep3Id)==13 && abs(GenLep4Id)==13) || (abs(GenLep1Id)==11 && abs(GenLep2Id)==11 && abs(GenLep3Id)==11 && abs(GenLep4Id)==11))"
            cutchan_gen_out  = "((abs(GenZ1Flav)==121 && abs(GenZ2Flav)==169) || (abs(GenZ1Flav)==169 && abs(GenZ2Flav)==121))"   
            cutm4l_gen       = "(GenHMass>"+str(m4l_low)+" && GenHMass<"+str(m4l_high)+")"
        
        cuth4l_gen  = "1==1"
        cutnoth4l_gen  = "(!"+cuth4l_gen+")"
 
        shortname = sample_shortnames[Sample]
        processBin = shortname+'_'+channel+'_'+opt.OBSNAME+'_genbin'+str(genbin)
                                    
        # GEN level        
        Histos[processBin+"fs"] = TH1D(processBin+"fs", processBin+"fs", 100, -1, 10000)
        Histos[processBin+"fs"].Sumw2()
        #TString toCut = ""
        gen_sumWeights = str(sumw[Sample])
	print(gen_sumWeights)
	lumi = '59.7'
        toCut = "((" + nnloWeights[0] + "*1000*" + lumi + "*xsec*genHEPMCweight*PUWeight)/"+gen_sumWeights + ")*(" + cutchan_gen_out + ")"

        Tree[Sample].Draw("GenHMass >> "+processBin+"fs",toCut,"goff")
        #else:
        #    Tree[Sample].Draw("GenHMass >> "+processBin+"fs","(qcdWeights[0])*("+cutchan_gen_out+")","goff")
        
	for i in range(0,9):
            if (i==5 or i==7): continue 
            Histos[processBin+"fs"+str(i)] = TH1D(processBin+"fs"+str(i), processBin+"fs"+str(i), 100, -1, 10000)
            Histos[processBin+"fs"+str(i)].Sumw2()

            toCut = "((" + nnloWeights[i] + "*1000*" + lumi +"*xsec*genHEPMCweight*PUWeight)/"+gen_sumWeights +  ")*("+cutchan_gen_out+")"

            Tree[Sample].Draw("GenHMass >> "+processBin+"fs"+str(i),toCut,"goff")

            Histos[processBin+"fid"+str(i)] = TH1D(processBin+"fid"+str(i), processBin+"fid"+str(i), m4l_bins, m4l_low, m4l_high)  
            Histos[processBin+"fid"+str(i)].Sumw2()

            toCut = "((" + nnloWeights[i] + "*1000*" + lumi + "*xsec*genHEPMCweight*PUWeight)/"+gen_sumWeights + ")*(passedFiducialSelection==1 && "+cutm4l_gen+" && "+cutobs_gen+" && "+cutchan_gen+"  && "+cuth4l_gen+")"

            Tree[Sample].Draw("GenHMass >> "+processBin+"fid"+str(i),toCut,"goff")
            Histos[processBin+"fid"+str(i)].Scale(1.0/Histos[processBin+"fs"].Integral())

            Histos[processBin+"fidraw"+str(i)] = TH1D(processBin+"fidraw"+str(i), processBin+"fidraw"+str(i), m4l_bins, m4l_low, m4l_high)  
            Histos[processBin+"fidraw"+str(i)].Sumw2()
            Tree[Sample].Draw("GenHMass >> "+processBin+"fidraw"+str(i),toCut,"goff")
            Histos[processBin+"fidraw"+str(i)].Scale(1.0/Histos[processBin+"fs"+str(i)].Integral())
            Histos[processBin+"fs"+str(i)].Scale(1.0/Histos[processBin+"fs"+str(i)].Integral())

        Histos[processBin+"fidPDF_up"] = TH1D(processBin+"fidPDF_up", processBin+"fidPDF_up", m4l_bins, m4l_low, m4l_high)
        Histos[processBin+"fidPDF_up"].Sumw2()
        toCut = "((LHEweight_PDFVariation_Up*1000*" + lumi+"*xsec*genHEPMCweight*PUWeight)/"+gen_sumWeights +")*(passedFiducialSelection==1 && "+cutm4l_gen+" && "+cutobs_gen+" && "+cutchan_gen+"  && "+cuth4l_gen+")"
        Tree[Sample].Draw("GenHMass >> "+processBin+"fidPDF_up",toCut,"goff")
        Histos[processBin+"fidPDF_up"].Scale(1.0/Histos[processBin+"fs"].Integral())

        Histos[processBin+"fidPDF_dn"] = TH1D(processBin+"fidPDF_dn", processBin+"fidPDF_dn", m4l_bins, m4l_low, m4l_high)
        Histos[processBin+"fidPDF_dn"].Sumw2()
        toCut = "((LHEweight_PDFVariation_Dn*1000*" + lumi+"*xsec*genHEPMCweight*PUWeight)/"+gen_sumWeights +")*(passedFiducialSelection==1 && "+cutm4l_gen+" && "+cutobs_gen+" && "+cutchan_gen+"  && "+cuth4l_gen+")"
        Tree[Sample].Draw("GenHMass >> "+processBin+"fidPDF_dn",toCut,"goff")
        Histos[processBin+"fidPDF_dn"].Scale(1.0/Histos[processBin+"fs"].Integral())

        fsintegral = Histos[processBin+"fs"].Integral()
        Histos[processBin+"fs"].Scale(1.0/Histos[processBin+"fs"].Integral())

        # GEN level 
        accerrstat=0.0
        if (Histos[processBin+"fs"].Integral()>0):
            # if ("NNLOPS" in processBin):
            print(Histos[processBin+"fs"].Integral(),Histos[processBin+"fid0"].Integral())
            acceptance[processBin] = Histos[processBin+"fid0"].Integral()/Histos[processBin+"fs"].Integral()
            #accerrstat = sqrt(acceptance[processBin]*(1-acceptance[processBin])/fsintegral)
            qcderrup=1.0; qcderrdn=1.0; 
            accerrup=1.0; accerrdn=1.0;
            print(processBin+'fid0', Histos[processBin+"fid0"].Integral())
            #for i in range(9,36):
            for i in range(0,9):
                #if (i==14 or i==16 or i==23 or i==25 or i==32 or i==34): continue                                        
                #if (i==5 or i==7 or i==14 or i==16 or i==23 or i==25): continue                                        
                if (i==5 or i==7): continue 
                ratio = Histos[processBin+"fid"+str(i)].Integral()/Histos[processBin+"fid0"].Integral()
                print(i,'ratio',ratio)
                if (ratio>qcderrup): qcderrup = Histos[processBin+"fid"+str(i)].Integral()/Histos[processBin+"fid0"].Integral()
                if (ratio<qcderrdn): qcderrdn = Histos[processBin+"fid"+str(i)].Integral()/Histos[processBin+"fid0"].Integral()

                acci = Histos[processBin+"fidraw"+str(i)].Integral()/Histos[processBin+"fs"+str(i)].Integral()
                print(i,"acc",acci)
                print(Histos[processBin+"fidraw"+str(i)].Integral(),Histos[processBin+"fs"+str(i)].Integral())
                if (acci/acceptance[processBin]>accerrup): accerrup=acci/acceptance[processBin]
                if (acci/acceptance[processBin]<accerrdn): accerrdn=acci/acceptance[processBin]

            qcdUncert[processBin] = {"uncerDn":abs(qcderrdn-1.0),"uncerUp":abs(qcderrup-1.0)}
            pdferr_up = Histos[processBin+"fidPDF_up"].Integral()/Histos[processBin+"fid0"].Integral()
            pdferr_dn = Histos[processBin+"fidPDF_dn"].Integral()/Histos[processBin+"fid0"].Integral()
            pdfUncert[processBin] = {"uncerDn":abs(pdferr_dn-1.0),"uncerUp":abs(pdferr_up-1.0)}

            print(processBin,acceptance[processBin],accerrstat,qcderrup,qcderrdn,pdferr_up,pdferr_dn)
            print("accerrup",accerrup,"accerrdn",accerrdn)
            
m4l_bins = 35
m4l_low = 105.0
m4l_high = 140.0

# Default to inclusive cross section
obs_reco = 'ZZMass'
obs_gen = 'GenHMass'
obs_reco_low = 105.0
obs_reco_high = 140.0
obs_gen_low = 105.0
obs_gen_high = 140.0

if (opt.OBSNAME == "massZ1"):
    obs_reco = "massZ1"
    obs_gen = "GENmZ1"
if (opt.OBSNAME == "massZ2"):
    obs_reco = "massZ2"
    obs_gen = "GENmZ2"
if (opt.OBSNAME == "pT4l"):
    obs_reco = "ZZPt"
    obs_gen = "GenHPt"
if (opt.OBSNAME == "eta4l"):
    obs_reco = "eta4l"
    obs_gen = "GENeta4l"
if (opt.OBSNAME == "njets_pt30_eta4p7"):
    obs_reco = "njets_pt30_eta4p7"
    obs_gen = "GENnjets_pt30_eta4p7"
if (opt.OBSNAME == "njets_pt30_eta2p5"):
    obs_reco = "njets_pt30_eta2p5"
    obs_gen = "GENnjets_pt30_eta2p5"
if (opt.OBSNAME == "pt_leadingjet_pt30_eta4p7"):
    obs_reco = "pt_leadingjet_pt30_eta4p7"
    obs_gen = "GENpt_leadingjet_pt30_eta4p7"
if (opt.OBSNAME == "pt_leadingjet_pt30_eta2p5"):
    obs_reco = "pt_leadingjet_pt30_eta2p5"
    obs_gen = "GENpt_leadingjet_pt30_eta2p5"
if (opt.OBSNAME == "rapidity4l"):
    obs_reco = "abs(rapidity4l)"
    obs_gen = "abs(GENrapidity4l)"
if (opt.OBSNAME == "cosThetaStar"):
    obs_reco = "abs(cosThetaStar)"
    obs_gen = "abs(GENcosThetaStar)"
if (opt.OBSNAME == "cosTheta1"):
    obs_reco = "abs(cosTheta1)"
    obs_gen = "abs(GENcosTheta1)"
if (opt.OBSNAME == "cosTheta2"):
    obs_reco = "abs(cosTheta2)"
    obs_gen = "abs(GENcosTheta2)"
if (opt.OBSNAME == "Phi"):
    obs_reco = "abs(Phi)"
    obs_gen = "abs(GENPhi)"    
if (opt.OBSNAME == "Phi1"):
    obs_reco = "abs(Phi1)"
    obs_gen = "abs(GENPhi1)"
    
#obs_bins = {0:(opt.OBSBINS.split("|")[1:((len(opt.OBSBINS)-1)/2)]),1:['0','inf']}[opt.OBSNAME=='inclusive'] 
obs_bins = opt.OBSBINS.split("|") 
if (not (obs_bins[0] == '' and obs_bins[len(obs_bins)-1]=='')): 
    print('BINS OPTION MUST START AND END WITH A |') 
obs_bins.pop()
obs_bins.pop(0) 

List = []
for long, short in sample_shortnames.iteritems():
    # if (not "ggH" in short): continue
    print(long, short)
    List.append(long)

if (obs_reco=="ZZMass"):
    chans = ['4e','4mu','2e2mu', '4l']
else:
    #chans = ['4e','4mu']
    chans = ['4e','4mu','2e2mu', '4l']

for chan in chans:
    for genbin in range(len(obs_bins)-1):
        getunc(chan,List, m4l_bins, m4l_low, m4l_high, obs_reco, obs_gen, obs_bins, genbin)  

if (obs_reco.startswith("njets")):
    for chan in chans:
        for genbin in range(len(obs_bins)-2): # last bin is >=3
            for Sample in List:
                shortname = sample_shortnames[Sample]
                processBin = shortname+'_'+chan+'_'+obs_reco+'_genbin'+str(genbin)
                processBinPlus1 = shortname+'_'+chan+'_'+obs_reco+'_genbin'+str(genbin+1)
                acceptance[processBin] = acceptance[processBin]-acceptance[processBinPlus1]
                qcdUncert[processBin]['uncerUp'] = sqrt(qcdUncert[processBin]['uncerUp']*qcdUncert[processBin]['uncerUp']+qcdUncert[processBinPlus1]['uncerUp']*qcdUncert[processBinPlus1]['uncerUp'])
                qcdUncert[processBin]['uncerDn'] = sqrt(qcdUncert[processBin]['uncerDn']*qcdUncert[processBin]['uncerDn']+qcdUncert[processBinPlus1]['uncerDn']*qcdUncert[processBinPlus1]['uncerDn'])

os.system('cp accUnc_'+opt.OBSNAME+'.py cp accUnc_'+opt.OBSNAME+'_ORIG.py')
with open('accUnc_'+opt.OBSNAME+'.py', 'w') as f:
    f.write('acc = '+str(acceptance)+' \n')
    f.write('qcdUncert = '+str(qcdUncert)+' \n')
    f.write('pdfUncert = '+str(pdfUncert)+' \n')
