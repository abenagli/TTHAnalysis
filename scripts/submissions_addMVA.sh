./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l ttgjets \
    -I /afs/cern.ch/work/a/abeschi/public/4Davide2017/Ntuples2017/ttH_Bkg_Generic_v2/           -i output_TTGJets\*.root   -t TTHGenericTagDumper/trees/ttGJets_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA/   -o TTGJets \
    --submit

./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l ttjets \
    -I /afs/cern.ch/work/a/abeschi/public/4Davide2017/Ntuples2017/ttH_Bkg_Generic_v2/           -i output_TTJets\*.root   -t TTHGenericTagDumper/trees/ttJets_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA/   -o TTJets \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l gg \
    -I /afs/cern.ch/work/a/abeschi/public/4Davide2017/Ntuples2017/ttH_Bkg_Generic_v2/           -i output_DiPhotonJetsBox\*.root   -t TTHGenericTagDumper/trees/diphoton_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA/   -o DiPhotonJetsBox \
    -n 20 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l gjet \
    -I /afs/cern.ch/work/a/abeschi/public/4Davide2017/Ntuples2017/ttH_Bkg_Generic_v2/           -i output_GJet\*.root   -t TTHGenericTagDumper/trees/gjet_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA/   -o GJet \
    -n 10 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l qcd \
    -I /afs/cern.ch/work/a/abeschi/public/4Davide2017/Ntuples2017/ttH_Bkg_Generic_v2/           -i output_QCD\*.root   -t TTHGenericTagDumper/trees/qcd_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA/   -o QCD \
    -n 10 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l fake-fake \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples/           -i fake-fake.root   -t fake-fake \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA/   -o fake-fake \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l prompt-fake \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples/           -i prompt-fake.root   -t prompt-fake \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA/   -o prompt-fake \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l ttH \
    -I /afs/cern.ch/work/a/abeschi/public/4Davide2017/Ntuples2017/ttH_Signal_Generic_v4/           -i ttHJetToGG_M125_13TeV_\*.root   -t TTHGenericTagDumper/trees/tth_125_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA/   -o ttHJetToGG_M125_13TeV \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l ggH \
    -I /afs/cern.ch/work/a/abeschi/public/4Davide2017/Ntuples2017/ttH_Signal_Generic_v4/           -i GluGluHToGG_M125_13TeV_\*.root   -t TTHGenericTagDumper/trees/ggh_125_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA/   -o GluGluHToGG_M125_13TeV \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l VBF \
    -I /afs/cern.ch/work/a/abeschi/public/4Davide2017/Ntuples2017/ttH_Signal_Generic_v4/           -i VBFHToGG_M125_13TeV\*.root   -t TTHGenericTagDumper/trees/vbf_125_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA/   -o VBFHToGG_M125_13TeV \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l VH \
    -I /afs/cern.ch/work/a/abeschi/public/4Davide2017/Ntuples2017/ttH_Signal_Generic_v4/           -i VHToGG_M125_13TeV_\*.root   -t TTHGenericTagDumper/trees/wzh_125_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA/   -o VHToGG_M125_13TeV \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l data \
    -I /afs/cern.ch/work/a/abeschi/public/4Davide2017/Ntuples2017/ttH_Data_Generic_v2/           -i output_DoubleEG_Run2017\*.root   -t TTHGenericTagDumper/trees/Data_13TeV_all \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA/   -o DoubleEG_Run2017 \
    -n 10 \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l CS_oneCategory \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/           -i controlSample_2017.root   -t ControlSample_oneCategory \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA/   -o controlSample_oneCategory \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l CS_diLepton \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/           -i controlSample_2017.root   -t ControlSample_diLepton \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA/   -o controlSample_diLepton \
    --submit


./sendGenericJob.py   -q 8nh   -b /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/   -e bin/addMVA.exe    -c scripts/addMVA_template.cfg \
    -l CS_singleLepton \
    -I /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/           -i controlSample_2017.root   -t ControlSample_singleLepton \
    -O /afs/cern.ch/work/a/abenagli/HGG/TTH/TTHAnalysis/data/ntuples_withMVA/   -o controlSample_singleLepton \
    --submit


