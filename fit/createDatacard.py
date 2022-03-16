import os,sys



def createDatacard(obsName, channel, nBins, obsBin, observableBins, physicalModel, year, nData, jes, lowerBound, upperBound, yearSetting):
    # Name of the bin (aFINALSTATE_ recobinX)
    if(channel == '4mu'): channelNumber = 1
    if(channel == '4e'): channelNumber = 2
    if(channel == '2e2mu'): channelNumber = 3
    # binName = 'a'+str(channelNumber)+'_recobin'+str(obsBin)

    # _recobin = str(observableBins[obsBin]).replace('.', 'p') + '_' + str(observableBins[obsBin+1]).replace('.', 'p')
    _recobin = str(observableBins[obsBin]).replace('.', 'p')+'_'+str(observableBins[obsBin+1]).replace('.', 'p')
    if int(observableBins[obsBin+1]) > 1000:
        _recobin = 'GT'+str(int(observableBins[obsBin]))

    _obsName = {'pT4l': 'PTH', 'rapidity4l': 'YH', 'pTj1': 'PTJET', 'njets_pt30_eta2p5': 'NJ'}
    binName = 'hzz_' + _obsName[obsName] + '_' + _recobin + '_cat' + channel
    # Root of the name of the process (signal from genBin)
    # processName = 'smH'+channel+'Bin'
    processName = 'smH_' + _obsName[obsName]
    # Background expectations in [105,160]
    sys.path.append('../inputs')
    _temp = __import__('inputs_bkgTemplate_'+obsName, globals(), locals(), ['expected_yield'], -1)
    expected_yield = _temp.expected_yield
    if jes:
        sys.path.append('../coefficients/JES')
        jesNames = ['Abs', 'Abs_year', 'BBEC1', 'BBEC1_year', 'EC2', 'EC2_year', 'FlavQCD', 'HF', 'HF_year', 'RelBal', 'RelSample_year']
        jesNames_datacard = [j.replace('year',year) for j in jesNames] # The name of the nuisance in the datacard should have the correspoding year
        _temp = __import__('JESNP_'+obsName+'_'+str(year), globals(), locals(), ['JESNP'], -1)
        jesnp = _temp.JESNP
        sys.path.remove('../coefficients/JES')
        # print(jesnp)
    sys.path.remove('../inputs')

    #Hard coded values for bkgs
    # bkg_qqzz = {}
    # bkg_qqzz['2016_2e2mu'] = 26.079487245654125
    # bkg_qqzz['2016_4e'] = 9.416905347804663
    # bkg_qqzz['2016_4mu'] = 22.327126146498454
    # bkg_qqzz['2017_2e2mu'] = 29.824908413 #14.998985147 #29.649379102635237
    # bkg_qqzz['2017_4e'] = 10.1053020211 #5.15676840688 #9.912444551671058
    # bkg_qqzz['2017_4mu'] = 26.3103880671 #12.6333277479 #26.57920447504215
    # bkg_qqzz['2018_2e2mu'] = 43.16035102064523
    # bkg_qqzz['2018_4e'] = 13.033793192674157
    # bkg_qqzz['2018_4mu'] = 38.774732781504355
    #
    # bkg_ggzz = {}
    # bkg_ggzz['2016_2e2mu'] = 2.1354071862105046
    # bkg_ggzz['2016_4e'] = 1.2047085561080217
    # bkg_ggzz['2016_4mu'] = 2.4823777973777252
    # bkg_ggzz['2017_2e2mu'] = 1.99032561607 #0.965204958226 #2.334396026294916
    # bkg_ggzz['2017_4e'] = 1.26621968326 #0.565844597635 #1.2379253364867422
    # bkg_ggzz['2017_4mu'] = 2.8809450136  #1.42919245518 #2.8757178358247417
    # bkg_ggzz['2018_2e2mu'] = 3.38606853989503
    # bkg_ggzz['2018_4e'] = 1.9682465482920186
    # bkg_ggzz['2018_4mu'] = 4.483239853533689
    #
    # bkg_zx = {}
    # bkg_zx['2016_2e2mu'] = 12.57494261733606
    # bkg_zx['2016_4e'] = 3.295497544812088
    # bkg_zx['2016_4mu'] = 9.065862569934886
    # bkg_zx['2017_2e2mu'] = 12.07 #12.344443134966895
    # bkg_zx['2017_4e'] = 3.05 #3.1395112978473105
    # bkg_zx['2017_4mu'] = 10.13 #10.874994231420253
    # bkg_zx['2018_2e2mu'] = 18.93941736832197
    # bkg_zx['2018_4e'] = 4.426316103466345
    # bkg_zx['2018_4mu'] = 17.05146905705311

    # lumi
    if yearSetting == 'Full':
        lumi = {}
        lumi['2016'] = '1.01'
        lumi['2017'] = '1.02'
        lumi['2018'] = '1.015'
        lumi_corr_16_17_18 = {}
        lumi_corr_16_17_18['2016'] = '1.006'
        lumi_corr_16_17_18['2017'] = '1.009'
        lumi_corr_16_17_18['2018'] = '1.02'
        lumi_corr_17_18 = {}
        lumi_corr_17_18['2017'] = '1.006'
        lumi_corr_17_18['2018'] = '1.002'
    else:
        lumi = {}
        lumi['2016'] = '1.026'
        lumi['2017'] = '1.025'
        lumi['2018'] = '1.023'

    # Lepton efficiency
    eff_mu = {}
    eff_mu['2016_2e2mu'] = '0.983/1.012'#'1.025'
    eff_mu['2016_4mu'] = '0.977/1.016' #'0.953/1.046'
    eff_mu['2017_2e2mu'] = '0.985/1.008'#'0.985/1.008'
    eff_mu['2017_4mu'] = '0.98/1.011'#'0.98/1.011'
    eff_mu['2018_2e2mu'] = '0.986/1.007'#'0.988/1.01'
    eff_mu['2018_4mu'] = '0.981/1.01'#'0.976/1.018'

    eff_e = {}
    eff_e['2016_2e2mu'] = '0.893/1.105'#'0.96/1.039'
    eff_e['2016_4e'] = '0.835/1.155'#'0.914/1.082'
    eff_e['2017_2e2mu'] = '0.915/1.082'#'0.915/1.082'
    eff_e['2017_4e'] = '0.867/1.121'#'0.867/1.121'
    eff_e['2018_2e2mu'] = '0.923/1.074'#'0.928/1.074'
    eff_e['2018_4e'] = '0.877/1.11'#'0.850/1.161'

    # ZX
    ZX = {}
    ZX['2016_2e2mu'] = '0.65673/1.35484'
    ZX['2016_4e'] = '0.60745/1.42863'
    ZX['2016_4mu'] = '0.69481/1.30542'
    ZX['2017_2e2mu'] = '0.67262/1.33282'
    ZX['2017_4e'] = '0.63816/1.37505'
    ZX['2017_4mu'] = '0.69350/1.30685'
    ZX['2018_2e2mu'] = '0.67618/1.32828'
    ZX['2018_4e'] = '0.64540/1.36539'
    ZX['2018_4mu'] = '0.69559/1.30459'


    # -------------------------------------------------------------------------------------------------

    file = open('../datacard/datacard_'+year+'/hzz4l_'+channel+'S_13TeV_xs_'+obsName+'_bin'+str(obsBin)+'_'+physicalModel+'.txt', 'w+')

    file.write('imax 1 \n')
    file.write('jmax * \n')
    file.write('kmax * \n')

    file.write('------------ \n')

    file.write('shapes * * hzz4l_'+channel+'S_13TeV_xs_SM_125_'+obsName+'_'+physicalModel+'.Databin'+str(obsBin)+'.root w:$PROCESS\n')

    file.write('------------ \n')

    file.write('bin '+binName+'\n')
    file.write('observation '+str(nData)+'\n')

    file.write('------------ \n')
    file.write('## mass window ['+str(lowerBound)+','+str(upperBound)+']\n')
    file.write('bin ')
    # for i in range(nBins+5): # In addition to the observableBins, there are OutsideAcceptance, fakeH, bkg_ggzz, bkg_qqzz, bkg_zjets
    for i in range(nBins+4):
        file.write(binName+' ')
    file.write('\n')
    file.write('process ')
    for i in range(nBins):
        # file.write(processName+str(i)+' ')
        if observableBins[i+1] > 1000:
            file.write(processName+'_GT'+str(int(observableBins[i]))+' ')
        else:
            file.write(processName+'_'+str(observableBins[i]).replace('.', 'p')+'_'+str(observableBins[i+1]).replace('.', 'p')+' ')
    # file.write('OutsideAcceptance nonResH bkg_qqzz bkg_ggzz bkg_zjets')
    file.write('nonResH bkg_qqzz bkg_ggzz bkg_zjets')
    file.write('\n')
    file.write('process ')
    for i in range(nBins):
        file.write('-'+str(i+1)+' ')
    # file.write('1 2 3 4 5')
    file.write('1 2 3 4')
    file.write('\n')
    file.write('rate ')
    # for i in range(nBins+2): # In addition to the observableBins, there are OutsideAcceptance, fakeH
    for i in range(nBins+1):
        file.write('1.0 ')
    # file.write(str(bkg_qqzz[year+'_'+channel])+' '+str(bkg_ggzz[year+'_'+channel])+' '+str(bkg_zx[year+'_'+channel])+'\n') #Old implementation with hard coding bkg expectation values
    file.write(str(expected_yield[int(year),'qqzz',channel])+' '+str(expected_yield[int(year),'ggzz',channel])+' '+str(expected_yield[int(year),'ZX',channel])+'\n')
    file.write('------------ \n')

    # norm_fake
    file.write('norm_nonResH lnU ')
    # for i in range(nBins+1): # Signal + OutsideAcceptance
    for i in range(nBins): # Signal + OutsideAcceptance
        file.write('- ')
    file.write('10.0 - - -    # [/10,*10]\n')

    if yearSetting == 'Full':
        # lumi_uncorrelated
        file.write('lumi_13TeV_'+year+'_uncorrelated lnN ')
        # for i in range(nBins+4): # All except ZX
        for i in range(nBins+3): # All except ZX
            file.write(lumi[year]+' ')
        file.write('-\n') # ZX
        # lumi_correlated_16_17_18
        file.write('lumi_13TeV_correlated_16_17_18 lnN ')
        # for i in range(nBins+4): # All except ZX
        for i in range(nBins+3): # All except ZX
            file.write(lumi_corr_16_17_18[year]+' ')
        file.write('-\n') # ZX
        # lumi_correlated_17_18
        if year == '2017' or year == '2018':
            file.write('lumi_13TeV_correlated_17_18 lnN ')
            # for i in range(nBins+4): # All except ZX
            for i in range(nBins+3): # All except ZX
                file.write(lumi_corr_17_18[year]+' ')
            file.write('-\n') # ZX
    else:
        # lumi
        file.write('lumi_13TeV_'+year+' lnN ')
        # for i in range(nBins+4): # All except ZX
        for i in range(nBins+3): # All except ZX
            file.write(lumi[year]+' ')
        file.write('-\n') # ZX

    # Lepton efficiency
    if channel == '4mu' or channel == '2e2mu':
        file.write('CMS_eff_m lnN ')
        # for i in range(nBins+4): # All except ZX
        for i in range(nBins+3): # All except ZX
            file.write(eff_mu[year+'_'+channel]+' ')
        file.write('-\n') # ZX
    if channel == '4e' or channel == '2e2mu':
        file.write('CMS_eff_e lnN ')
        # for i in range(nBins+4): # All except ZX
        for i in range(nBins+3): # All except ZX
            file.write(eff_e[year+'_'+channel]+' ')
        file.write('-\n') # ZX

    # ZX
    file.write('CMS_hzz'+channel+'_Zjets_'+year+' lnN ')
    # for i in range(nBins+4): # All except ZX
    for i in range(nBins+3): # All except ZX
        file.write('- ')
    file.write(ZX[year+'_'+channel]+'\n')

    # Param
    if(channelNumber != 2):
        file.write('CMS_zz4l_mean_m_sig_'+year+' param 0.0 1.0\n')
        file.write('CMS_zz4l_sigma_m_sig_'+year+' param 0.0 0.2 [-1,1]\n')
    if(channelNumber != 1):
        file.write('CMS_zz4l_mean_e_sig_'+year+' param 0.0 1.0\n')
        file.write('CMS_zz4l_sigma_e_sig_'+year+' param 0.0 0.2 [-1,1]\n')

    file.write('CMS_zz4l_n_sig_'+str(channelNumber)+'_'+year+' param 0.0 0.05\n')

    # Theoretical
    file.write('QCDscale_ggVV lnN ')
    # for i in range(nBins+3): # Signal + out + fake + qqzz
    for i in range(nBins+2): # All except ZX
        file.write('- ')
    file.write('1.039/0.961 -\n')
    file.write('QCDscale_VV lnN ')
    # for i in range(nBins+2): # Signal + out + fake
    for i in range(nBins+1):
        file.write('- ')
    file.write('1.0325/0.958 - -\n')
    file.write('pdf_gg lnN ')
    # for i in range(nBins+3): # Signal + out + fake + qqzz
    for i in range(nBins+2): # All except ZX
        file.write('- ')
    file.write('1.032/0.968 -\n')
    file.write('pdf_qqbar lnN ')
    # for i in range(nBins+2): # Signal + out + fake
    for i in range(nBins+1): # All except ZX
        file.write('- ')
    file.write('1.031/0.966 - -\n')
    file.write('kfactor_ggzz lnN ')
    # for i in range(nBins+3): # Signal + out + fake  + bkg_qqzz
    for i in range(nBins+2): # Signal + out + fake  + bkg_qqzz
        file.write('- ')
    file.write('1.1 -\n')

    # # ZX
    # file.write('CMS_zjets_bkgdcompo_'+str(year)+' lnN ')
    # for i in range(nBins+4): # All except ZX
    #     file.write('- ')
    # file.write('1.32\n')

    # JES
    if jes == True:
        for index,jesName in enumerate(jesNames_datacard):
            file.write('CMS_scale_j_'+jesName+' lnN ')
            for i in range(nBins): # Signals
                file.write(str(jesnp['fiducial_'+jesNames[index]+'_'+channel+'_'+obsName+'_genbin'+str(i)+'_recobin'+str(obsBin)])+' ')
            file.write(str(jesnp['nonFiducial_'+jesNames[index]+'_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])+' ')
            file.write(str(jesnp['nonResonant_'+jesNames[index]+'_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])+' ')
            file.write(str(jesnp['qqzz_'+jesNames[index]+'_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])+' ')
            file.write(str(jesnp['ggzz_'+jesNames[index]+'_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])+' ')
            file.write('-\n') # ZX
        file.write('CMS_scale_j_ZX lnN ')
        for i in range(nBins+3): # All except ZX
            file.write('- ')
        file.write(str(jesnp['ZX_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])+'\n')

            #file.write('JES param 0.0 1.0\n')

    file.close()

    print(os.getcwd())

def createDatacard_ggH(obsName, channel, nBins, obsBin, observableBins, physicalModel, year, nData, jes, lowerBound, upperBound, yearSetting):
    # Name of the bin (aFINALSTATE_ recobinX)
    if(channel == '4mu'): channelNumber = 1
    if(channel == '4e'): channelNumber = 2
    if(channel == '2e2mu'): channelNumber = 3
    # binName = 'a'+str(channelNumber)+'_recobin'+str(obsBin)

    # _recobin = str(observableBins[obsBin]).replace('.', 'p') + '_' + str(observableBins[obsBin+1]).replace('.', 'p')
    _recobin = str(observableBins[obsBin]).replace('.', 'p')+'_'+str(observableBins[obsBin+1]).replace('.', 'p')
    if int(observableBins[obsBin+1]) > 1000:
        _recobin = 'GT'+str(int(observableBins[obsBin]))

    _obsName = {'pT4l': 'PTH', 'rapidity4l': 'YH', 'pTj1': 'PTJET', 'njets_pt30_eta2p5': 'NJ'}
    binName = 'hzz_' + _obsName[obsName] + '_' + _recobin + '_cat' + channel
    # Root of the name of the process (signal from genBin)
    # processName = 'smH'+channel+'Bin'
    processName = 'ggH_' + _obsName[obsName]
    xHName = 'xH_' + _obsName[obsName]
    # Background expectations in [105,160]
    sys.path.append('../inputs')
    _temp = __import__('inputs_bkgTemplate_'+obsName, globals(), locals(), ['expected_yield'], -1)
    expected_yield = _temp.expected_yield
    if jes:
        sys.path.append('../coefficients/JES')
        jesNames = ['Abs', 'Abs_year', 'BBEC1', 'BBEC1_year', 'EC2', 'EC2_year', 'FlavQCD', 'HF', 'HF_year', 'RelBal', 'RelSample_year']
        jesNames_datacard = [j.replace('year',year) for j in jesNames] # The name of the nuisance in the datacard should have the correspoding year
        _temp = __import__('JESNP_'+obsName+'_'+str(year), globals(), locals(), ['JESNP'], -1)
        jesnp = _temp.JESNP
        sys.path.remove('../coefficients/JES')
        # print(jesnp)
    sys.path.remove('../inputs')

    # lumi
    if yearSetting == 'Full':
        lumi = {}
        lumi['2016'] = '1.01'
        lumi['2017'] = '1.02'
        lumi['2018'] = '1.015'
        lumi_corr_16_17_18 = {}
        lumi_corr_16_17_18['2016'] = '1.006'
        lumi_corr_16_17_18['2017'] = '1.009'
        lumi_corr_16_17_18['2018'] = '1.02'
        lumi_corr_17_18 = {}
        lumi_corr_17_18['2017'] = '1.006'
        lumi_corr_17_18['2018'] = '1.002'
    else:
        lumi = {}
        lumi['2016'] = '1.026'
        lumi['2017'] = '1.025'
        lumi['2018'] = '1.023'

    # Lepton efficiency
    eff_mu = {}
    eff_mu['2016_2e2mu'] = '0.983/1.012'#'1.025'
    eff_mu['2016_4mu'] = '0.977/1.016' #'0.953/1.046'
    eff_mu['2017_2e2mu'] = '0.985/1.008'#'0.985/1.008'
    eff_mu['2017_4mu'] = '0.98/1.011'#'0.98/1.011'
    eff_mu['2018_2e2mu'] = '0.986/1.007'#'0.988/1.01'
    eff_mu['2018_4mu'] = '0.981/1.01'#'0.976/1.018'

    eff_e = {}
    eff_e['2016_2e2mu'] = '0.893/1.105'#'0.96/1.039'
    eff_e['2016_4e'] = '0.835/1.155'#'0.914/1.082'
    eff_e['2017_2e2mu'] = '0.915/1.082'#'0.915/1.082'
    eff_e['2017_4e'] = '0.867/1.121'#'0.867/1.121'
    eff_e['2018_2e2mu'] = '0.923/1.074'#'0.928/1.074'
    eff_e['2018_4e'] = '0.877/1.11'#'0.850/1.161'

    # ZX
    ZX = {}
    ZX['2016_2e2mu'] = '0.65673/1.35484'
    ZX['2016_4e'] = '0.60745/1.42863'
    ZX['2016_4mu'] = '0.69481/1.30542'
    ZX['2017_2e2mu'] = '0.67262/1.33282'
    ZX['2017_4e'] = '0.63816/1.37505'
    ZX['2017_4mu'] = '0.69350/1.30685'
    ZX['2018_2e2mu'] = '0.67618/1.32828'
    ZX['2018_4e'] = '0.64540/1.36539'
    ZX['2018_4mu'] = '0.69559/1.30459'


    # -------------------------------------------------------------------------------------------------

    file = open('../datacard/datacard_'+year+'/hzz4l_GGH_'+channel+'S_13TeV_xs_'+obsName+'_bin'+str(obsBin)+'_'+physicalModel+'.txt', 'w+')

    file.write('imax 1 \n')
    file.write('jmax * \n')
    file.write('kmax * \n')

    file.write('------------ \n')

    file.write('shapes * * hzz4l_'+channel+'S_13TeV_xs_SM_125_'+obsName+'_'+physicalModel+'.Databin'+str(obsBin)+'.root w:$PROCESS\n')

    file.write('------------ \n')

    file.write('bin '+binName+'\n')
    file.write('observation '+str(nData)+'\n')

    file.write('------------ \n')
    file.write('## mass window ['+str(lowerBound)+','+str(upperBound)+']\n')
    file.write('bin ')
    # for i in range(nBins+5): # In addition to the observableBins, there are OutsideAcceptance, fakeH, bkg_ggzz, bkg_qqzz, bkg_zjets
    for i in range(2*nBins+4):
        file.write(binName+' ')
    file.write('\n')
    file.write('process ')
    for i in range(nBins):
        # file.write(processName+str(i)+' ')
        if observableBins[i+1] > 1000:
            file.write(processName+'_GT'+str(int(observableBins[i]))+' ')
        else:
            file.write(processName+'_'+str(observableBins[i]).replace('.', 'p')+'_'+str(observableBins[i+1]).replace('.', 'p')+' ')
    for i in range(nBins):
        # file.write(processName+str(i)+' ')
        if observableBins[i+1] > 1000:
            file.write(xHName+'_GT'+str(int(observableBins[i]))+' ')
        else:
            file.write(xHName+'_'+str(observableBins[i]).replace('.', 'p')+'_'+str(observableBins[i+1]).replace('.', 'p')+' ')
    # file.write('OutsideAcceptance nonResH bkg_qqzz bkg_ggzz bkg_zjets')
    file.write('nonResH bkg_qqzz bkg_ggzz bkg_zjets')
    file.write('\n')
    file.write('process ')
    for i in range(nBins):
        file.write('-'+str(i+1)+' ')
    # nonResH bkg_qqzz bkg_ggzz bkg_zjets
    file.write('1 2 3 4 ')
    # xH as bkg processBin
    for i in range(nBins):
        file.write(str(i+5)+' ')
    file.write('\n')
    file.write('rate ')
    # ggH, XH bins and nonRes contribution
    for i in range(2*nBins+1):
        file.write('1.0 ')
    # file.write(str(bkg_qqzz[year+'_'+channel])+' '+str(bkg_ggzz[year+'_'+channel])+' '+str(bkg_zx[year+'_'+channel])+'\n') #Old implementation with hard coding bkg expectation values
    file.write(str(expected_yield[int(year),'qqzz',channel])+' '+str(expected_yield[int(year),'ggzz',channel])+' '+str(expected_yield[int(year),'ZX',channel])+'\n')
    file.write('------------ \n')

    # norm_fake
    file.write('norm_nonResH lnU ')
    # ggH + xH
    for i in range(2*nBins):
        file.write('- ')
    file.write('10.0 - - -    # [/10,*10]\n')

    if yearSetting == 'Full':
        # lumi_uncorrelated
        file.write('lumi_13TeV_'+year+'_uncorrelated lnN ')
        # for i in range(nBins+4): # All except ZX
        for i in range(2*nBins+3): # All except ZX
            file.write(lumi[year]+' ')
        file.write('-\n') # ZX
        # lumi_correlated_16_17_18
        file.write('lumi_13TeV_correlated_16_17_18 lnN ')
        # for i in range(nBins+4): # All except ZX
        for i in range(2*nBins+3): # All except ZX
            file.write(lumi_corr_16_17_18[year]+' ')
        file.write('-\n') # ZX
        # lumi_correlated_17_18
        if year == '2017' or year == '2018':
            file.write('lumi_13TeV_correlated_17_18 lnN ')
            # for i in range(nBins+4): # All except ZX
            for i in range(2*nBins+3): # All except ZX
                file.write(lumi_corr_17_18[year]+' ')
            file.write('-\n') # ZX
    else:
        # lumi
        file.write('lumi_13TeV_'+year+' lnN ')
        # for i in range(nBins+4): # All except ZX
        for i in range(2*nBins+3): # All except ZX
            file.write(lumi[year]+' ')
        file.write('-\n') # ZX

    # Lepton efficiency
    if channel == '4mu' or channel == '2e2mu':
        file.write('CMS_eff_m lnN ')
        # for i in range(nBins+4): # All except ZX
        for i in range(2*nBins+3): # All except ZX
            file.write(eff_mu[year+'_'+channel]+' ')
        file.write('-\n') # ZX
    if channel == '4e' or channel == '2e2mu':
        file.write('CMS_eff_e lnN ')
        # for i in range(nBins+4): # All except ZX
        for i in range(2*nBins+3): # All except ZX
            file.write(eff_e[year+'_'+channel]+' ')
        file.write('-\n') # ZX

    # ZX
    file.write('CMS_hzz'+channel+'_Zjets_'+year+' lnN ')
    # for i in range(nBins+4): # All except ZX
    for i in range(2*nBins+3): # All except ZX
        file.write('- ')
    file.write(ZX[year+'_'+channel]+'\n')

    # Param
    if(channelNumber != 2):
        file.write('CMS_zz4l_mean_m_sig_'+year+' param 0.0 1.0\n')
        file.write('CMS_zz4l_sigma_m_sig_'+year+' param 0.0 0.2 [-1,1]\n')
    if(channelNumber != 1):
        file.write('CMS_zz4l_mean_e_sig_'+year+' param 0.0 1.0\n')
        file.write('CMS_zz4l_sigma_e_sig_'+year+' param 0.0 0.2 [-1,1]\n')

    file.write('CMS_zz4l_n_sig_'+str(channelNumber)+'_'+year+' param 0.0 0.05\n')

    # Theoretical
    file.write('QCDscale_ggVV lnN ')
    # for i in range(nBins+3): # Signal + out + fake + qqzz
    for i in range(2*nBins+2): # All except ZX
        file.write('- ')
    file.write('1.039/0.961 -\n')
    file.write('QCDscale_VV lnN ')
    # for i in range(nBins+2): # Signal + out + fake
    for i in range(2*nBins+1):
        file.write('- ')
    file.write('1.0325/0.958 - -\n')
    file.write('pdf_gg lnN ')
    # for i in range(nBins+3): # Signal + out + fake + qqzz
    for i in range(2*nBins+2): # All except ZX
        file.write('- ')
    file.write('1.032/0.968 -\n')
    file.write('pdf_qqbar lnN ')
    # for i in range(nBins+2): # Signal + out + fake
    for i in range(2*nBins+1): # All except ZX
        file.write('- ')
    file.write('1.031/0.966 - -\n')
    file.write('kfactor_ggzz lnN ')
    # for i in range(nBins+3): # Signal + out + fake  + bkg_qqzz
    for i in range(2*nBins+2): # Signal + out + fake  + bkg_qqzz
        file.write('- ')
    file.write('1.1 -\n')

    # # ZX
    # file.write('CMS_zjets_bkgdcompo_'+str(year)+' lnN ')
    # for i in range(nBins+4): # All except ZX
    #     file.write('- ')
    # file.write('1.32\n')

    # JES
    if jes == True:
        for index,jesName in enumerate(jesNames_datacard):
            file.write('CMS_scale_j_'+jesName+' lnN ')
            for i in range(nBins): # Signals
                file.write(str(jesnp['fiducial_'+jesNames[index]+'_'+channel+'_'+obsName+'_genbin'+str(i)+'_recobin'+str(obsBin)])+' ')
            file.write(str(jesnp['nonFiducial_'+jesNames[index]+'_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])+' ')
            file.write(str(jesnp['nonResonant_'+jesNames[index]+'_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])+' ')
            file.write(str(jesnp['qqzz_'+jesNames[index]+'_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])+' ')
            file.write(str(jesnp['ggzz_'+jesNames[index]+'_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])+' ')
            file.write('-\n') # ZX
            #file.write('JES param 0.0 1.0\n')

    file.close()

    print(os.getcwd())
