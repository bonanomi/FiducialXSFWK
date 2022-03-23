list = {'rapidity4l': '|0.0|0.15|0.3|0.6|0.9|1.2|2.5|',#'|0.0|0.15|0.3|0.45|0.6|0.75|0.9|1.2|1.6|2.0|2.5|',
'pT4l': '|0|10|20|30|45|60|80|120|200|13000|',
'pT4l_kL': '|0.0|45.0|80.0|120.0|200.0|1300.0|',
'massZ1': '|50|64|73|85|106|',
'massZ2': '|12|20|24|28|32|40|55|65|',
'mass4l': '|105|140|',
'njets_pt30_eta2p5': '|0|1|2|3|4|14|',
'pTj1': '|-2|30|55|95|200|1300|',
'njets_pt30_eta2p5 vs pTj1': '|0|1|10| vs |0|1| / |30|60|120|350|',
'njets_pt30_eta4p7 vs pTj1_eta4p7': '|0|1|10| vs |0|1| / |30|60|120|350|',
'njets_pt30_eta2p5 vs pTj2': '|0|2|10| vs |0|1| / |30|60|120|350|',
'costhetastar': '|0.0|0.2|0.4|0.6|0.8|1.0|', #|0.0|0.125|0.25|0.375|0.5|0.625|0.75|0.875|1.0|
'costhetaZ1': '|-1.0|-0.75|-0.50|-0.25|0.0|0.25|0.50|0.75|1.0|', #'|0.0|0.2|0.4|0.6|0.8|1.0|'
'costhetaZ2': '|-1.0|-0.75|-0.50|-0.25|0.0|0.25|0.50|0.75|1.0|', #'|0.0|0.2|0.4|0.6|0.8|1.0|'
'phi': '|-3.14159265359|-2.35619449019|-1.57079632679|-0.785398163397|0.0|0.785398163397|1.57079632679|2.35619449019|3.14159265359|',
'phistar': '|-3.14159265359|-2.35619449019|-1.57079632679|-0.785398163397|0.0|0.785398163397|1.57079632679|2.35619449019|3.14159265359|',
'TCjmax': '|0|15|30|50|70|90|1000|', ## Should treat as 2D? In principle yes, leading jet involved
'TBjmax': '|0|30|60|80|100|120|1000|', ## Should treat as 2D? In principle yes, leading jet involved
'pTHj': '|-2|0|30|70|1000|',
'rapidity4l vs pT4l': '|0.0|0.5| vs |0|40| / |0.0|0.5| vs |40|80| / |0.0|0.5| vs |80|150| / |0.0|0.5| vs |150|1300| / |0.5|1.0| vs |0|45| / |0.5|1.0| vs |45|120| / |0.5|1.0| vs |120|1300| / |1.0|2.5| vs |0|45| / |1.0|2.5| vs |45|120| / |1.0|2.5| vs |120|1300|',
'njets_pt30_eta2p5 vs pT4l': '|0|1|2|3|20| vs |0|15|30|120|350| / |0|15|30|120|350| / |0|120|350| / |0|120|350|',
'massZ1 vs massZ2': '|50|80| vs |10|30| / |50|80| vs |30|60| / |80|110| vs |10|25| / |80|110| vs |25|30|',
'pT4l vs pTj1': '|0|350| vs |0|30| / |0|80| vs |30|60| / |80|350| vs |30|60| / |0|120| vs |60|120| / |120|350| vs |60|120| / |0|120| vs |120|350| / |120|350| vs |120|350|', ## Tricky, it is actually a 3D measurement?
#'pT4l vs pTHj': '|0|350| vs |0|30|', ## Tricky, it is actually a 3D measurement?
# 'mHj vs pTHj': , ## Tricky, it is actually a 3D measurement?
'njets_pt30_eta2p5 vs pTHj': '|0|1|10| vs |0|800| / |0|60|120|350|',
'njets_pt30_eta2p5 vs pTHjj': '|0|2|10| vs |0|800| / |0|60|120|350|',
'njets_pt30_eta2p5 vs mjj': '|0|2|10| vs |0|800| / |0|120|450|3000|', # Add bins of mjj
'njets_pt30_eta2p5 vs mHj': '|0|1|10| vs |0|800| / |120|180|220|300|400|600|2000|',
'njets_pt30_eta2p5 vs mHjj': '|0|2|10| vs |0|800| / |180|320|450|600|1000|2500|', # Add bins of mHjj
'njets_pt30_eta2p5 vs absdetajj': '|0|2|10| vs |0|800| / |0.0|1.0|2.5|9.0|', # Add bins of absdetajj
'pTj1 vs pTj2': '|0|30| vs |0|30| / |30|60| vs |30|60| / |60|120| vs |60|120| / |120|350| vs |120|350|',
'njets_pt30_eta2p5 vs TBjmax': '|0|1|10| vs |-1|1| / |0|30|60|80|100|120|1000|',
'njets_pt30_eta2p5 vs TCjmax': '|0|1|10| vs |-1|1| / |0|15|30|50|70|90|1000|'
# 'njets_pt30_eta2p5 vs absdphijj': '|0|2|10| vs |0|800| / ||' # Decide if we want to go to 2PI or to PI
}


def binning(var):
    obsBins_input = list[var]
    print('input', obsBins_input)
    if not 'vs' in obsBins_input: #It is not a double-differential analysis
        obs_bins = {0:(obsBins_input.split("|")[1:(len(obsBins_input.split("|"))-1)]),1:['0','inf']}[obsBins_input=='inclusive']
        obs_bins = [float(i) for i in obs_bins] #Convert a list of str to a list of float
        doubleDiff = False
        print ('It is a single-differential measurement, binning', obs_bins)
    else: #It is a double-differential analysis
        doubleDiff = True
        # The structure of obs_bins is:
        # index of the dictionary is the number of the bin
        # [obs_bins_low, obs_bins_high, obs_bins_low_2nd, obs_bins_high_2nd]
        # The first two entries are the lower and upper bound of the first variable
        # The second two entries are the lower and upper bound of the second variable
        if obsBins_input.count('vs')==1 and obsBins_input.count('/')>=1: #Situation like this one '|0|1|2|3|20| vs |0|10|20|45|90|250| / |0|10|20|80|250| / |0|20|90|250| / |0|25|250|'
            obs_bins_tmp = obsBins_input.split(" vs ") #['|0|1|2|3|20|', '|0|10|20|45|90|250| / |0|10|20|80|250| / |0|20|90|250| / |0|25|250|']
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
        elif obsBins_input.count('vs')>1 and obsBins_input.count('/')>1: #Situation like this one '|50|80| vs |10|30| / |50|80| vs |30|60| / |80|110| vs |10|25| / |80|110| vs |25|30|'
            obs_bins_tmp = obsBins_input.split(' / ') #['|50|80| vs |10|30|', '|50|80| vs |30|60|', '|80|110| vs |10|25|', '|80|110| vs |25|30|']
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
        elif obsBins_input.count('vs')==1 and obsBins_input.count('/')==0: #Situation like this one '|0|1|2|3|20| vs |0|10|20|45|90|250|'
            obs_bins_tmp = obsBins_input.split(" vs ") #['|0|1|2|3|20|', '|0|10|20|45|90|250|']
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
            print ('Problem in the definition of the binning')
            quit()
        print ('It is a double-differential measurement, binning for the 1st variable', obs_bins_1st, 'and for the 2nd variable', obs_bins_2nd)
        print (obs_bins)
    return obs_bins, doubleDiff
