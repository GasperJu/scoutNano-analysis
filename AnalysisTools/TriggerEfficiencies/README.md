# Dijet analysis with scouting data: trigger efficiencies

## Trigger efficiency code

Setup a CMSSW working environment, e.g.:
```
cmsrel CMSSW_14_0_12/
cd CMSSW_14_0_12/src/
cmsenv
```
and copy there the directory provided here (`2024_UtilsDataQuality`), which includes the HLT Jet Energy Corrections json file, included in the analysis script via `correctionlib`, and the Golden json.

A set of scripts to run the trigger efficiency code both for data and MC is provided in the `condor_sub` directory. You can directly test on a single file the code running e.g.:
```
python3 condor_sub/mc_triggEff.py 
python3 condor_sub/data2024_triggEffInclusive.py
```

The code at the moment compute trigger efficiencies inclusively for the JetHT scouting path (`PFScouting_JetHT`), using the scouting single-muon trigger (`PFScouting_SingleMuon`) as a reference, looking at H_T, leading jet pt, and dijet mass as variables.

A script is also provided to study different categories of events:
* VBF-like topology;
* boosted + ISR;
* resolved + ISR;
* rest of pp events;

as a function of:
* m_jj;
* H_T

and it can be tested as it follows:

```                                                                                                                                                            
python3 condor_sub/data2024_triggEff.py
```                                                                                                                                                            

## Condor submission

A setup for condor submission to study trigger efficiency on partial or full dataset is provided in the repository `condor_sub`.
The ingredients (that can be customised according to the needs) are:
* A `data2024_triggEffInclusive.submit` for condor submission. Note: err/out/log files are disabled as they occupy a lot of space and the fill up the disk available in afs quickly.
* A `run_data2024_triggEffInclusive.sh` as executable running the python script to compute trigger efficiencies for every file in the chosen list (see below).
* The python script to compute trigger efficiencies `data2024_triggEffInclusive.py`: note that the only difference with respect to the script provided in the main directory is the parsing of some options to allow to run on every single file separately.
* A list of files to consider. Two examples are provided: a reduced list (`red_list_data2024H.txt`) and the full list for 2024H data (`list_data2024H.txt`). Note that the prefix `root://cms-xrd-global.cern.ch/` (or `root://eoscms.cern.ch//eos/cms` when files are on eos) is needed.

For submission:

* Setup the proxy:

```
voms-proxy-init --voms cms --valid 168:00 
cp /tmp/x509up_u101050 ~/
export X509_USER_PROXY=~/x509up_u101050
```
or:
```
voms-proxy-init --voms cms --valid 168:00 -out $HOME/private/.proxy
export X509_USER_PROXY=$HOME/private/.proxy
```

* Create an output directory to collect the produced histo files, e.g. `inclTrigEff_outputHistos`

* Run:
```
condor_submit data2024_triggEffInclusive.submit
```

Once the jobs are terminated, proceed to add all files together.

Some plotting scripts are provided in the repository:
* For inclusive trigger efficiencies in data: `plotDijet_trigEffIncl.py`
* For inclusive trigger efficiencies in MC (note: proper weighting of the combined samples not included, useful for single-sample plots): `plotDijet_mcTrigEffIncl.py`
* For data trigger efficiencies obtained for each category separately: `plotDijetHTEff.py` 
* For data and MC comparison: `compare_data_mc_eff.py`, see example for running it below. 

```
python3 compare_data_mc_eff.py --data /eos/cms/store/cmst3/user/elfontan/scoutAna/TriggerEff/histos_inclTrgEffData_2024H_wJECs.root --mc histos_TT4Q_JECs.root:762.1 histos_Wto2Q_JECs.root:16100  histos_QCD-HT200to400_JECs.root:1951000.0 histos_QCD-HT400to600_JECs.root:96660.0 histos_QCD-HT600to800_JECs.root:13684.0 histos_QCD-HT800to1000_JECs.root:3047.0 --mc-indir /eos/cms/store/cmst3/user/elfontan/scoutAna/TriggerEff --lumi-pb 1000 --variables ht_inclusive --vline 280 --rebin 2

python3 compare_data_mc_eff.py --data /eos/cms/store/cmst3/user/elfontan/scoutAna/TriggerEff/histos_inclTrgEffData_2024H_wJECs.root --mc histos_TT4Q_JECs.root:762.1 histos_Wto2Q_JECs.root:16100  histos_QCD-HT200to400_JECs.root:1951000.0 histos_QCD-HT400to600_JECs.root:96660.0 histos_QCD-HT600to800_JECs.root:13684.0 histos_QCD-HT800to1000_JECs.root:3047.0 --mc-indir /eos/cms/store/cmst3/user/elfontan/scoutAna/TriggerEff --lumi-pb 1000 --variables pt_leading --vline 180

python3 compare_data_mc_eff.py --data /eos/cms/store/cmst3/user/elfontan/scoutAna/TriggerEff/histos_inclTrgEffData_2024H_wJECs.root --mc histos_TT4Q_JECs.root:762.1 histos_Wto2Q_JECs.root:16100  histos_QCD-HT200to400_JECs.root:1951000.0 histos_QCD-HT400to600_JECs.root:96660.0 histos_QCD-HT600to800_JECs.root:13684.0 histos_QCD-HT800to1000_JECs.root:3047.0 --mc-indir /eos/cms/store/cmst3/user/elfontan/scoutAna/TriggerEff --lumi-pb 1000 --variables dijetMass --vline -1 --rebin 2
```
