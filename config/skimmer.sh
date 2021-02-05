#for year in 2016 2017 2018
#do
#    for process in ggH125 VBFH125 WplusH125 WminusH125 ZH125 ttH125
#    do
#	root -q -b skim_MC_tree.C"(\"$process\",\"$year\")"
#    done
#done

#for year in 2016 2017 2018
#do
#    for process in ggH125 VBFH125 WplusH125 WminusH125 ZH125 ttH125 ZZTo4lext ggTo2e2mu_Contin_MCFM701 ggTo2e2tau_Contin_MCFM701 ggTo2mu2tau_Contin_MCFM701 ggTo4e_Contin_MCFM701 ggTo4mu_Contin_MCFM701 ggTo4tau_Contin_MCFM701
#    do
#        root -q -b skim_MC_tree.C"(\"$process\",\"$year\")"
#    done
#done

for year in 2016
do
   for process in ggH0M_M125
   do
       root -q -b skim_MC_tree.C"(\"$process\",\"$year\")"
   done
done

# for year in 2016 2017 2018
# do
#    for process in ggH0M_M125 ttH0M_M125 VBFH0M_M125 WH0M_M125 ZH0M_M125
#    do
#        root -q -b skim_MC_tree.C"(\"$process\",\"$year\")"
#    done
# done
