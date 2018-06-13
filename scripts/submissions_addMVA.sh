./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l ttgjets \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_9_4_6/src/flashgg/Jobs/ntuples_HggPresel_bkg_genericTag_v2/           -i TTGJets\*.root   -t TTHGenericTagDumper/trees/ttGJets_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA_genericTag_v2/   -o TTGJets \
    --submit

./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l ttjets \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_9_4_6/src/flashgg/Jobs/ntuples_HggPresel_bkg_genericTag_v2/           -i TTJets\*.root   -t TTHGenericTagDumper/trees/ttJets_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA_genericTag_v2/   -o TTJets \
    -n 2 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l gg \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_9_4_6/src/flashgg/Jobs/ntuples_HggPresel_bkg_genericTag_v2/           -i DiPhotonJetsBox\*.root   -t TTHGenericTagDumper/trees/diphoton_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA_genericTag_v2/   -o DiPhotonJetsBox \
    -n 20 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l gjet \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_9_4_6/src/flashgg/Jobs/ntuples_HggPresel_bkg_genericTag_v2/           -i GJet\*.root   -t TTHGenericTagDumper/trees/gjet_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA_genericTag_v2/   -o GJet \
    -n 20 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l qcd \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_9_4_6/src/flashgg/Jobs/ntuples_HggPresel_bkg_genericTag_v2/           -i QCD\*.root   -t TTHGenericTagDumper/trees/qcd_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA_genericTag_v2/   -o QCD \
    -n 10 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l ttH \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_9_4_6/src/flashgg/Jobs/ntuples_HggPresel_signal_genericTag_v2/           -i ttHJetToGG_M125_13TeV_\*.root   -t TTHGenericTagDumper/trees/tth_125_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA_genericTag_v2/   -o ttHJetToGG_M125_13TeV \
    -n 4 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l ggH \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_9_4_6/src/flashgg/Jobs/ntuples_HggPresel_signal_genericTag_v2/           -i GluGluHToGG_M125_13TeV_\*.root   -t TTHGenericTagDumper/trees/ggh_125_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA_genericTag_v2/   -o GluGluHToGG_M125_13TeV \
    -n 8 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l VBF \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_9_4_6/src/flashgg/Jobs/ntuples_HggPresel_signal_genericTag_v2/           -i VBFHToGG_M125_13TeV\*.root   -t TTHGenericTagDumper/trees/vbf_125_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA_genericTag_v2/   -o VBFHToGG_M125_13TeV \
    -n 4 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l VH \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_9_4_6/src/flashgg/Jobs/ntuples_HggPresel_signal_genericTag_v2/           -i VHToGG_M125_13TeV_\*.root   -t TTHGenericTagDumper/trees/vh_125_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA_genericTag_v2/   -o VHToGG_M125_13TeV \
    -n 2 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l data \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/CMSSW_9_4_6/src/flashgg/Jobs/ntuples_HggPresel_data_genericTag_v2/           -i DoubleEG_Run2017\*.root   -t TTHGenericTagDumper/trees/Data_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA_genericTag_v2/   -o DoubleEG_Run2017 \
    -n 40 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l CS_diLepton \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/           -i controlSample_2017_invBTag.root   -t ControlSample_diLepton \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA_genericTag_v2/   -o controlSample_diLepton_invBTag \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l CS_singleLepton \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/           -i controlSample_2017_invBTag.root   -t ControlSample_singleLepton \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA_genericTag_v2/   -o controlSample_singleLepton_invBTag \
    --submit
