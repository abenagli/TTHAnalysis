#include "interface/TreeUtils.h"



void InitTreeVars(TChain* chain, TreeVars& treeVars)
{
  chain -> SetBranchAddress("run", &treeVars.run);
  chain -> SetBranchAddress("lumi", &treeVars.lumi);
  chain -> SetBranchAddress("event", &treeVars.event);
  
  chain -> SetBranchAddress("nvtx", &treeVars.nvtx);
  chain -> SetBranchAddress("weight", &treeVars.weight);
  
  chain -> SetBranchAddress("dipho_sumpt", &treeVars.dipho_sumpt);
  chain -> SetBranchAddress("dipho_mass", &treeVars.dipho_mass);
  chain -> SetBranchAddress("dipho_vtxProb", &treeVars.dipho_vtxProb);
  chain -> SetBranchAddress("dipho_sigmaRV", &treeVars.dipho_sigmaRV);
  chain -> SetBranchAddress("dipho_sigmaWV", &treeVars.dipho_sigmaWV);
  chain -> SetBranchAddress("dipho_deltaphi", &treeVars.dipho_deltaphi);
  chain -> SetBranchAddress("dipho_cosDeltaphi", &treeVars.dipho_cosDeltaphi);
  chain -> SetBranchAddress("dipho_leadPt", &treeVars.dipho_leadPt);
  chain -> SetBranchAddress("dipho_leadEta", &treeVars.dipho_leadEta);
  chain -> SetBranchAddress("dipho_leadPhi", &treeVars.dipho_leadPhi);
  chain -> SetBranchAddress("dipho_leadEnergy", &treeVars.dipho_leadEnergy);
  chain -> SetBranchAddress("dipho_leadR9", &treeVars.dipho_leadR9);
  chain -> SetBranchAddress("dipho_lead_ptoM", &treeVars.dipho_lead_ptoM);
  chain -> SetBranchAddress("dipho_lead_sigmaEoE", &treeVars.dipho_lead_sigmaEoE);
  chain -> SetBranchAddress("dipho_leadIDMVA", &treeVars.dipho_leadIDMVA);
  chain -> SetBranchAddress("dipho_lead_PSV", &treeVars.dipho_lead_PSV);
  chain -> SetBranchAddress("dipho_subleadPt", &treeVars.dipho_subleadPt);
  chain -> SetBranchAddress("dipho_subleadEta", &treeVars.dipho_subleadEta);
  chain -> SetBranchAddress("dipho_subleadPhi", &treeVars.dipho_subleadPhi);
  chain -> SetBranchAddress("dipho_subleadEnergy", &treeVars.dipho_subleadEnergy);
  chain -> SetBranchAddress("dipho_subleadR9", &treeVars.dipho_subleadR9);
  chain -> SetBranchAddress("dipho_sublead_ptoM", &treeVars.dipho_sublead_ptoM);
  chain -> SetBranchAddress("dipho_sublead_sigmaEoE", &treeVars.dipho_sublead_sigmaEoE);
  chain -> SetBranchAddress("dipho_subleadIDMVA", &treeVars.dipho_subleadIDMVA);
  chain -> SetBranchAddress("dipho_sublead_PSV", &treeVars.dipho_sublead_PSV);
  chain -> SetBranchAddress("dipho_mva", &treeVars.dipho_mva);
  // chain -> SetBranchAddress("dipho_mva_training2016_best", &treeVars.dipho_mva_training2016_best);
  // chain -> SetBranchAddress("dipho_mva_training2017_v1", &treeVars.dipho_mva_training2017_v1);
  // chain -> SetBranchAddress("dipho_mva_training2017_bug", &treeVars.dipho_mva_training2017_bug);
  
  chain -> SetBranchAddress("nJets", &treeVars.nJets);
  chain -> SetBranchAddress("nJets_bTagLoose", &treeVars.nJets_bTagLoose);
  chain -> SetBranchAddress("nJets_bTagMedium", &treeVars.nJets_bTagMedium);
  chain -> SetBranchAddress("nJets_bTagTight", &treeVars.nJets_bTagTight);
  
  // chain -> SetBranchAddress("MetPt", &treeVars.MetPt);
  // chain -> SetBranchAddress("MetPhi", &treeVars.MetPhi);
  
  for(int i = 1; i <= nJet; i++)
  {
    chain -> SetBranchAddress(("jet_pt"+ std::to_string(i)).c_str(), &treeVars.jet_pt[i-1]);
    chain -> SetBranchAddress(("jet_eta"+ std::to_string(i)).c_str(), &treeVars.jet_eta[i-1]);
    chain -> SetBranchAddress(("jet_phi"+ std::to_string(i)).c_str(), &treeVars.jet_phi[i-1]);
    chain -> SetBranchAddress(("jet_bdiscriminantDeep"+ std::to_string(i)).c_str(), &treeVars.jet_bdiscriminant[i-1]);
  }
  
  for(int i = 1; i <= nLep; i++)
  {
    chain -> SetBranchAddress(("mu_pt"+ std::to_string(i)).c_str(), &treeVars.mu_pt[i-1]);
    chain -> SetBranchAddress(("mu_eta"+ std::to_string(i)).c_str(), &treeVars.mu_eta[i-1]);
    chain -> SetBranchAddress(("mu_phi"+ std::to_string(i)).c_str(), &treeVars.mu_phi[i-1]);
    chain -> SetBranchAddress(("mu_energy"+ std::to_string(i)).c_str(), &treeVars.mu_energy[i-1]);
    chain -> SetBranchAddress(("mu_isLoose"+ std::to_string(i)).c_str(),&treeVars.mu_IDVector[i-1][0]);
    chain -> SetBranchAddress(("mu_isMedium"+ std::to_string(i)).c_str(),&treeVars.mu_IDVector[i-1][1]);
    chain -> SetBranchAddress(("mu_isTight"+ std::to_string(i)).c_str(),&treeVars.mu_IDVector[i-1][2]);
    chain -> SetBranchAddress(("mu_MiniIso"+ std::to_string(i)).c_str(), &treeVars.mu_miniIso[i-1]);
    chain -> SetBranchAddress(("mu_charge"+ std::to_string(i)).c_str(), &treeVars.mu_charge[i-1]);
    chain -> SetBranchAddress(("mu_trackIso"+ std::to_string(i)).c_str(), &treeVars.mu_trackIso[i-1]);
    chain -> SetBranchAddress(("mu_sumChargedHadronPt"+ std::to_string(i)).c_str(), &treeVars.mu_sumChargedHadronPt[i-1]);
    chain -> SetBranchAddress(("mu_sumNeutralHadronEt"+ std::to_string(i)).c_str(), &treeVars.mu_sumNeutralHadronEt[i-1]);
    chain -> SetBranchAddress(("mu_sumPhotonEt"+ std::to_string(i)).c_str(), &treeVars.mu_sumPhotonEt[i-1]);
    chain -> SetBranchAddress(("mu_sumPUPt"+ std::to_string(i)).c_str(), &treeVars.mu_sumPUPt[i-1]);
    
    
    chain -> SetBranchAddress(("ele_pt"+ std::to_string(i)).c_str(), &treeVars.ele_pt[i-1]);
    chain -> SetBranchAddress(("ele_eta"+ std::to_string(i)).c_str(), &treeVars.ele_eta[i-1]);
    chain -> SetBranchAddress(("ele_phi"+ std::to_string(i)).c_str(), &treeVars.ele_phi[i-1]);
    chain -> SetBranchAddress(("ele_energy"+ std::to_string(i)).c_str(), &treeVars.ele_energy[i-1]);
    chain -> SetBranchAddress(("ele_SCeta"+ std::to_string(i)).c_str(), &treeVars.ele_SCeta[i-1]);
    chain -> SetBranchAddress(("ele_SCphi"+ std::to_string(i)).c_str(), &treeVars.ele_SCphi[i-1]);
    chain -> SetBranchAddress(("ele_passVetoId"+ std::to_string(i)).c_str(),&treeVars.ele_IDVector[i-1][0]);
    chain -> SetBranchAddress(("ele_passLooseId"+ std::to_string(i)).c_str(),&treeVars.ele_IDVector[i-1][1]);
    chain -> SetBranchAddress(("ele_passMediumId"+ std::to_string(i)).c_str(),&treeVars.ele_IDVector[i-1][2]);
    chain -> SetBranchAddress(("ele_passTightId"+ std::to_string(i)).c_str(),&treeVars.ele_IDVector[i-1][3]);
    chain -> SetBranchAddress(("ele_MVALooseId"+ std::to_string(i)).c_str(),&treeVars.ele_IDVector[i-1][4]);
    chain -> SetBranchAddress(("ele_MVAMediumId"+ std::to_string(i)).c_str(),&treeVars.ele_IDVector[i-1][5]);
    chain -> SetBranchAddress(("ele_MVATightId"+ std::to_string(i)).c_str(),&treeVars.ele_IDVector[i-1][6]);
    chain -> SetBranchAddress(("ele_MVALooseNoIsoId"+ std::to_string(i)).c_str(),&treeVars.ele_IDVector[i-1][7]);
    chain -> SetBranchAddress(("ele_MVAMediumNoIsoId"+ std::to_string(i)).c_str(),&treeVars.ele_IDVector[i-1][8]);
    chain -> SetBranchAddress(("ele_MVATightNoIsoId"+ std::to_string(i)).c_str(),&treeVars.ele_IDVector[i-1][9]);
    chain -> SetBranchAddress(("ele_MiniIso"+ std::to_string(i)).c_str(), &treeVars.ele_miniIso[i-1]);
    chain -> SetBranchAddress(("ele_ecalEnergy"+ std::to_string(i)).c_str(), &treeVars.ele_ecalEnergy[i-1]);
    chain -> SetBranchAddress(("ele_SCx"+ std::to_string(i)).c_str(), &treeVars.ele_SCx[i-1]);
    chain -> SetBranchAddress(("ele_SCy"+ std::to_string(i)).c_str(), &treeVars.ele_SCy[i-1]);
    chain -> SetBranchAddress(("ele_SCz"+ std::to_string(i)).c_str(), &treeVars.ele_SCz[i-1]);
    chain -> SetBranchAddress(("ele_charge"+ std::to_string(i)).c_str(), &treeVars.ele_charge[i-1]);
    chain -> SetBranchAddress(("ele_dEtaSCTrackAtVtx"+ std::to_string(i)).c_str(), &treeVars.ele_dEtaTrk[i-1]);
    chain -> SetBranchAddress(("ele_dPhiSCTrackAtVtx"+ std::to_string(i)).c_str(), &treeVars.ele_dPhiTrk[i-1]);
  }
}



void InitOutTree(TTree* tree, TreeVars& treeVars)
{
  tree -> Branch("weight",&treeVars.weight);
  
  tree -> Branch("dipho_mass",                  &treeVars.dipho_mass);
  tree -> Branch("dipho_sigmaRV",               &treeVars.dipho_sigmaRV);
  tree -> Branch("dipho_deltaphi",              &treeVars.dipho_deltaphi);
  tree -> Branch("dipho_mva",                   &treeVars.dipho_mva);
  // tree -> Branch("dipho_mva_training2016_best", &treeVars.dipho_mva_training2016_best);
  // tree -> Branch("dipho_mva_training2017_v1",   &treeVars.dipho_mva_training2017_v1);
  // tree -> Branch("dipho_mva_training2017_bug",  &treeVars.dipho_mva_training2017_bug);
  
  tree -> Branch("dipho_leadEta",      &treeVars.dipho_leadEta);
  tree -> Branch("dipho_leadPhi",      &treeVars.dipho_leadPhi);
  tree -> Branch("dipho_lead_ptoM",    &treeVars.dipho_lead_ptoM);
  tree -> Branch("dipho_lead_sigmaEoE",&treeVars.dipho_lead_sigmaEoE);
  tree -> Branch("dipho_leadIDMVA",    &treeVars.dipho_leadIDMVA);
  tree -> Branch("dipho_lead_PSV",     &treeVars.dipho_lead_PSV);
  
  tree -> Branch("dipho_subleadEta",      &treeVars.dipho_subleadEta);
  tree -> Branch("dipho_subleadPhi",      &treeVars.dipho_subleadPhi);
  tree -> Branch("dipho_sublead_ptoM",    &treeVars.dipho_sublead_ptoM);
  tree -> Branch("dipho_sublead_sigmaEoE",&treeVars.dipho_sublead_sigmaEoE);
  tree -> Branch("dipho_subleadIDMVA",    &treeVars.dipho_subleadIDMVA);
  tree -> Branch("dipho_sublead_PSV",     &treeVars.dipho_sublead_PSV);
  
  tree -> Branch("nJets",           &treeVars.nJets);
  tree -> Branch("nJets_bTagLoose", &treeVars.nJets_bTagLoose);
  tree -> Branch("nJets_bTagMedium",&treeVars.nJets_bTagMedium);
  tree -> Branch("nJets_bTagTight", &treeVars.nJets_bTagTight);
  
  tree -> Branch("mu1_pt",    &treeVars.mu1_pt);
  tree -> Branch("mu1_eta",   &treeVars.mu1_eta);
  tree -> Branch("ele1_pt",   &treeVars.ele1_pt);
  tree -> Branch("ele1_eta",  &treeVars.ele1_eta);
  tree -> Branch("jet1_pt",   &treeVars.jet1_pt);
  tree -> Branch("jet1_eta",  &treeVars.jet1_eta);
  tree -> Branch("jet1_bTag", &treeVars.jet1_bTag);
  tree -> Branch("jet2_pt",   &treeVars.jet2_pt);
  tree -> Branch("jet2_eta",  &treeVars.jet2_eta);
  tree -> Branch("jet2_bTag", &treeVars.jet2_bTag);
  tree -> Branch("jet3_pt",   &treeVars.jet3_pt);
  tree -> Branch("jet3_eta",  &treeVars.jet3_eta);
  tree -> Branch("jet3_bTag", &treeVars.jet3_bTag);
  tree -> Branch("jet4_pt",   &treeVars.jet4_pt);
  tree -> Branch("jet4_eta",  &treeVars.jet4_eta);
  tree -> Branch("jet4_bTag", &treeVars.jet4_bTag);
  tree -> Branch("bTag1",     &treeVars.bTag1);
  tree -> Branch("bTag2",     &treeVars.bTag2);
  tree -> Branch("bTag3",     &treeVars.bTag3);
  tree -> Branch("bTag4",     &treeVars.bTag4);
}

void InitOutTreeVars(TreeVars& treeVars)
{
  treeVars.mu1_pt = -1.;
  treeVars.mu1_eta = -5.;
  treeVars.ele1_pt = -1.;
  treeVars.ele1_eta = -5.;
  treeVars.jet1_pt = -1.;
  treeVars.jet1_eta = -5.;
  treeVars.jet1_bTag = -1.;
  treeVars.jet2_pt = -1.;
  treeVars.jet2_eta = -5.;
  treeVars.jet2_bTag = -1.;
  treeVars.jet3_pt = -1.;
  treeVars.jet3_eta = -5.;
  treeVars.jet3_bTag = -1.;
  treeVars.jet4_pt = -1.;
  treeVars.jet4_eta = -5.;
  treeVars.jet4_bTag = -1.;
  treeVars.bTag1 = -1.;
  treeVars.bTag2 = -1.;
  treeVars.bTag3 = -1.;
  treeVars.bTag4 = -1.;
}
