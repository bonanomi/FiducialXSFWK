def binning(obsBins_input):
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
