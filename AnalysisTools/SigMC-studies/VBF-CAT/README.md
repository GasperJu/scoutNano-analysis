# Run 3 scouting dijet analysis: signal studies 


## Gen-level study: MINIAODSIM

If ScoutingNano are not available yet, a script to perform sanity checks on VBF signals on MINIAODSIM is provided.

How to run:
``` 
python3  studyVBFSignals_miniaod.py --input "/eos/cms/store/cmst3/group/run3Scouting/LowMassDijetSearch/samples/2024/mc/PrivateProduction/SampleFactory/VBFHTo2B_Par-M-1000-W-0p014__chain_RunIII2024Summer24wmLHEGS-RunIII2024Summer24MiniAOD/SampleFactory/260416_160516/0000/RunIII2024Summer24MiniAOD*root" --out VBFHTo2B_M1000_genHiggs.root
```


## Gen-level study: ScoutingNano 

The same type of checks from the script above is implemented to run on ScoutingNano.

Some options can be specified from command line:
* mass;
* flavour can be specified: 2B, 2C, 2Q, 2Glu

How to run:
``` 
python3 studyVBFSignals_scoutNano.py --decay 2B --mass 300
``` 

A plotter to examine the output is provided:
``` 
python3 plot_genVBFHiggs_scoutNano.py --input VBFHTo2B_M300_genHiggs_scoutNano.root --outdir /eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/SignalOptimisationStudies/VBFCat/VBFHTo2B/GENLevel/M300 --label "VBFHTo2B M300 (2024)"
``` 


## Reco-level study (with and without gen-matching): ScoutingNano 

Moving to reco-level studies (with and without gen-matching), a script is provided to study the forward jets from the VBF topology and to optimize the reconstruction of the central dijet pair.

How to run:
``` 
python3 studyVBFReco_scoutNano.py --decay 2B --mass 500  --requireTriggerBaseline --centralPairAlgo leading  --centralEtaMax 2.5  --output VBFHTo2B_M500_recoVBF_scoutNano_leading_massDepSel.root
```

The option `requireTriggerBaseline` applies the trigger selection. There is also the possibility to simply save the tribber bit decision and apply it offline..

The algorithm for pairing the central dijet system can be chosen `--centralPairAlgo`:
* leading
* minDR
* minDPhi
* massClosest

The option `centralPairSelection` allows to choose between a newly tested `massDependent` approach and a more standard angular approach to apply a selection on the central pair:
* angular mode uses `dEta < args.centralPairDEtaMax and dPhi > args.centralPairDPhiMin`;
* massDependent mode uses `dR > dr_min and apt < apt_max`.

By default, a dR and leading/subleading jet compatibility with a dependency on the mass is applied:
* `DeltaR > dR_min(mass)`
* `A_pt < A_pt_max(mass)`, with A_pt = |pT1 - pT2| / (pT1 + pT2)

with reference values of:
* 300 GeV -> (1.5, 0.60)
* 500 GeV -> (2.0, 0.45)
* 1000 GeV -> (2.5, 0.35)

This behaviour can be overwritten:
* `--centralPairDRMin`
* `--centralPairPtImbalanceMax`


A plotter to examine the output is provided:
``` 
python3 plot_recoVBF_scoutNano.py --input VBFHTo2B_M500_recoVBF_scoutNano.root --outdir  /eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/SignalOptimisationStudies/VBFCat/VBFHTo2B//M500/RECOLevel --label "VBFHTo2B M500 (2024)"
``` 

Signal mass histograms (inclusive and matched only) are saved as output of the plotting step, e.g.: `VBF-dijetMass-Histos_ForFIT/vbf-m500-leading.root`, with a 5 GeV bining:
* `h_signal_peak_selected_5GeV`	
* `h_signal_peak_matched_5GeV`


## Standalone signal modeling

From the output of the previous step, signal peaks at different mass hypotheses can be studied together:

```
python3 plot_allSignals_VBF.py --indir VBF-dijetMass-Histos_ForFIT --algo leading  --outdir /eos/user/e/elfontan/www/dijetAnaRun3/SIGModelling/VBF-CAT/
```
and a basic signal modelling based on a double-sided Crystal Ball can be done:
```
python3 fit_signalPeak_VBF.py --indir VBF-dijetMass-Histos_ForFIT --algo leading --outdir /eos/user/e/elfontan/www/dijetAnaRun3/SIGModelling/VBF-CAT/
```