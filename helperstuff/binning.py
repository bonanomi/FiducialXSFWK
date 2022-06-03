list = {
'mass4l': '|105|160|',
'mass4l_zzfloating': '|105|160|',
'njets_pt30_eta4p7': '|0|1|2|3|4|14|',
'pT4l': '|0|10|20|30|45|60|80|120|200|3000|',
'pT4l_kL': '|0|45|80|120|200|1300|',
'rapidity4l': '|0.0|0.15|0.3|0.45|0.6|0.75|0.9|1.2|1.6|2.5|',
'costhetaZ1': '|-1.0|-0.75|-0.50|-0.25|0.0|0.25|0.50|0.75|1.0|',
'costhetaZ2': '|-1.0|-0.75|-0.50|-0.25|0.0|0.25|0.50|0.75|1.0|',
'phi': '|-3.14159265359|-2.35619449019|-1.57079632679|-0.785398163397|0.0|0.785398163397|1.57079632679|2.35619449019|3.14159265359|',
'phistar': '|-3.14159265359|-2.35619449019|-1.57079632679|-0.785398163397|0.0|0.785398163397|1.57079632679|2.35619449019|3.14159265359|',
'costhetastar': '|-1.0|-0.75|-0.50|-0.25|0.0|0.25|0.50|0.75|1.0|',
'massZ1': '|40|65|73|80|85|90|120|',
'massZ2': '|12|22|25|28|30|32|35|40|50|65|',
'pTj1': '|-2|30|55|95|200|1300|',
'pTHj': '|-2|0|30|60|110|1300|',
'mHj': '|-2|110|180|220|300|400|600|3000|',
'pTj2': '|-2|30|40|65|90|1300|',
'mjj': '|-2|0|120|300|1300|',
'absdetajj': '|-100|0|0.7|1.6|3.0|10.0|',
'dphijj': '|-100|-3.14159265359|-1.5707963267948966|0|1.5707963267948966|3.14159265359|',
'pTHjj': '|-2|0|20|60|1300|',
'TCjmax': '|-2|15|20|30|50|80|1000|',
'TBjmax': '|-2|30|70|130|250|400|1000|',
'D0m': '|0.0|0.4|0.5|0.6|0.7|0.8|0.9|1.0|',
'Dcp': '|-0.75|-0.25|-0.1|0|0.1|0.25|0.75|',
'D0hp': '|0.0|0.35|0.4|0.45|0.55|0.65|0.75|1.0|',
'Dint': '|0.0|0.7|0.8|0.9|0.95|1.0|',
'DL1': '|0.0|0.45|0.5|0.6|0.7|1.0|',
'DL1Zg': '|0.0|0.35|0.45|0.5|0.55|0.65|1.0|',
'rapidity4l vs pT4l': '|0|0.5|1.0|2.5| vs |0|40|80|150|1300| / |0|45|120|1300| / |0|45|120|1300|',
'njets_pt30_eta4p7 vs pT4l': '|0|1|2|20| vs |0|15|30|1300| / |0|60|80|120|1300| / |0|100|170|250|1300|',
'pTj1 vs pTj2': '|-2|30| vs |-2|30| / |30|60| vs |30|60| / |60|350| vs |30|60| / |60|350| vs |60|350|',
'pT4l vs pTHj': '|0|1500| vs |-2|0| / |0|85| vs |0|30| / |85|350| vs |0|45| / |0|85| vs |30|350| / |85|350| vs |45|350|',
'massZ1 vs massZ2': '|40|85| vs |12|35| / |40|70| vs |35|65| / |70|120| vs |35|65| / |85|120| vs |30|35| / |85|120| vs |24|30| / |85|120| vs |12|24|',
'TCjmax vs pT4l': '|-2|15|25|40|1300| vs |0|15|30|45|70|120|1300| / |0|120|1300| / |0|120|1300| / |0|200|1300|',
}

def binning(var):
    obsBins_input = list[var]
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
