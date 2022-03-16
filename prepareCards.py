import os

bins = [0.0,0.15,0.3,0.45,0.6,0.75,0.9,1.2,1.6,2.5]
obsName = 'rapidity4l'
fitName = 'YH'

card_name = 'card_run2_%s.txt' %fitName

cmd_combCards = 'combineCards.py '

for year in [2016, 2017, 2018]:
  for cat in ['4e', '4mu', '2e2mu']:
    for i in range(len(bins)-1):
      low = str(bins[i]).replace('.','p')
      high = str(bins[i+1]).replace('.','p')
      boundaries = low+'_'+high
      dc_name = 'datacard_%d/hzz4l_%sS_13TeV_xs_%s_bin%d_v3.txt ' %(year,cat,obsName,i)
      cmd_combCards += 'hzz_%s_%s_cat%s_%d=%s' %(fitName,boundaries,cat,year,dc_name)

cmd_combCards += '> %s' %card_name
# print(cmd_combCards)

cmd_t2w = 'text2workspace.py %s -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose ' %card_name
cmd_t2w += "--PO 'higgsMassRange=123,127' "

for i in range(len(bins)-1):
  low = str(bins[i]).replace('.','p')
  high = str(bins[i+1]).replace('.','p')
  boundaries = low+'_'+high
  process = 'smH_%s_%s' %(fitName, boundaries)
  POI = 'r_smH_%s_%d' %(fitName, i)
  cmd_t2w += "--PO 'map=.*/%s:%s[1.0,0.0,3.0]' " %(process, POI)

# print(cmd_t2w)

cmd_fit = 'combine -n _%s_Fit -M MultiDimFit %s ' %(fitName, card_name.replace('txt', 'root'))
cmd_fit += '-m 125.38 --freezeParameters MH --saveWorkspace --algo=singles --cminDefaultMinimizerStrategy 0 -t -1 --setParameters '

for i in range(len(bins)-1):
  POI = 'r_smH_%s_%d' %(fitName, i)
  cmd_fit += '%s=1,' %POI 

print(cmd_fit[:-1])