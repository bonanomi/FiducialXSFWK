for process in ggH125 VBFH125 WplusH125 WminusH125 ZH125 ttH125
do
	root -q -b skim_MC_tree.C"(\"$process\",\"2017\")"
done
