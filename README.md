# Fiducial XS measurements in HZZ4L

Standalone framework for fiducial differential cross section measurements using CJLST TTrees for Run 2 data.

A CMSSW working area with the latest version of [`combine`](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/) is required:

```
# cmssw-el7 (maybe needed)
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.2.0
scramv1 b clean; scramv1 b
cd ../..
```
Once that is done, combine is properly installed. The fiducial framework should be properly set up:
```
git clone https://github.com/AlessandroTarabini/FiducialXSFWK.git
cd FiducialXSFWK
source ./env
```
The `PhysicsModel(s)` used in this analysis (cf. `fit/createDatacard.py`) can be found in the `models` folder and are copied automatically when `source ./env.sh` is called to `$CMSSW_VERSION/src/HiggsAnalysis/CombinedLimit/python`.
If when running `fit/RunFiducialXS.py` the fit crashes because the physics model is not found, you should check that the `HZZ4L_Fiducial*.py` files are present in your `$CMSSW_VERSION/src/HiggsAnalysis/CombinedLimit/python` folder.

## Set Up
For all the imports in the various scripts to work properly, set up the working environment with:

```
cmsenv
export PYTHONPATH="${PYTHONPATH}:<PATH_TO>/FiducialXSFWK/inputs"
export PYTHONPATH="${PYTHONPATH}:<PATH_TO>/FiducialXSFWK/helperstuff"
```

## Workflow
A schematic representation of the framework's workflow is given in the two following sketches:

![Preparation of the reduced trees and datacards](fig/FiducialXS_Workflow.001.png)
![Setting up the datacards and running the fits](fig/FiducialXS_Workflow.002.png)

The input files of the analysis workflow are the HZZ4L ntuples generated with the [CJLST framework](https://github.com/AlessandroTarabini/ZZAnalysis/tree/Run2_CutBased_BTag16). This framework starts from those files and:

1. `config`: Starting from CJLST TTrees, only relevant branches are selected and stored by `skim_MC_tree.C` and `skim_data_tree.C` macros.

Having created these skimmed TTrees, the next steps of the analysis involve the caluclation of the different coefficients needed for the pdf parameterisations and unfolding, as well as the creation of background templates. To do so:

2. `templates`: Templates and normalization coefficients for the backgrounds' pdf are extracted from MC (ggZZ and qqZZ) and data (ZX) using [`RunTemplates.py`](https://github.com/bonanomi/FiducialXSFWK/tree/main/templates)
3. `coefficients`: All the coefficients of the signal parameterization are calculated with [`RunCoefficients.py`](https://github.com/bonanomi/FiducialXSFWK/blob/main/coefficients/RunCoefficients.py) and stored in `inputs` folder.
4. `fit`: The maximum likelihood fit is performed. This step relies on the [`RunFiducialXS.py`](https://github.com/bonanomi/FiducialXSFWK/blob/main/fit/RunFiducialXS.py) script and it can be run either as part of the entire framework, creating the datacards and workspaces from scratch, or using pre-existing datacars as input. Datacards are produced and stored in a `datacard` directory, while fit results (combine `.root` files) are stored in `combine_files` folder.

Additional scripts are provided to plot negative log-likelihood scans and to produce the usual differential xsec plots:

5. [`LHScans`](https://github.com/bonanomi/FiducialXSFWK/tree/main/LHScans): Likelihood scans are plotted, best-fit values and the corresponding uncertainties are calculated using `plotLHScans_compare.py`.
6. `producePlots.py`: Plot of unfolded differential xsec distributions.

### Commands: a gentle introduction
The frameworks picks up the variables to be used in a measurement and the binning to be used in the fit thanks to the dictionaries available in `helperstuff`.
More in detail, to define a new measurement (or to understand what is used in a current measurement):

* Define the variable name and the reco- and gen-level observables modifying the [`observables` `dict` in `helperstuff/observables.py`](https://github.com/bonanomi/FiducialXSFWK/blob/Run3/helperstuff/observables.py). The syntax is:
	```
	observables = {NAME: {"obs_reco": TBranch name, "obs_gen": TBranch name}}
	```
* Define the binning in [`helperstuff/binning.py`](https://github.com/bonanomi/FiducialXSFWK/blob/Run3/helperstuff/binning.py) following the conventions defined [here](https://github.com/bonanomi/FiducialXSFWK/issues/20) for 2D measurements (for 1D measurements the binning definition is intuitive).

1. Skimming of the `TTrees`: TODO. Usually provided centrally

2. **Creation of the templates**: in [`templates`](https://github.com/bonanomi/FiducialXSFWK/tree/Run3/templates) run (the framework detects automatically the binning and the observables to use):
	```
	python RunTemplates.py --obsName NAME (str) --year YEAR (str)
	```

3. **Computation of the coefficients**: in [`coefficients`](https://github.com/bonanomi/FiducialXSFWK/tree/Run3/coefficients) run (the framework detects automatically the binning and the observables to use):
	```
	python RunCoefficients.py --obsName NAME (str) --year YEAR (str)
	```

4. **Run the fits**: in [`fit`](https://github.com/bonanomi/FiducialXSFWK/tree/Run3/fit) run:
	```
	python RunFiducialXS.py --obsName NAME (str) --year YEAR (str) 
	```

5. **Plot the NLL scans**: in [`LHScans`](https://github.com/bonanomi/FiducialXSFWK/tree/Run3/LHScans):
	```
	python plot_LLScan.py --obsName NAME (str) --year YEAR (str)                                
	```