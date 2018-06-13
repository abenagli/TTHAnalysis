#ifndef TREE_UTILS_H
#define TREE_UTILS_H

#include <iostream>
#include <vector>
#include <map>

#include "TChain.h"

#define nLep 8
#define nJet 10



/*** tree variables ***/
struct TreeVars
{
  unsigned int run;
  unsigned int lumi;
  unsigned long long int event;
  
  int nvtx;
  float weight;
  
  float dipho_sumpt;
  float dipho_mass;
  float dipho_vtxProb;
  float dipho_sigmaRV;
  float dipho_sigmaWV;
  float dipho_deltaphi;
  float dipho_cosDeltaphi;
  float dipho_leadPt;
  float dipho_leadEta;
  float dipho_leadPhi;
  float dipho_leadEnergy;
  float dipho_leadR9;
  float dipho_lead_ptoM;
  float dipho_lead_sigmaEoE;
  float dipho_leadIDMVA;
  float dipho_lead_PSV;
  float dipho_leadGenMatch;
  float dipho_subleadPt;
  float dipho_subleadEta;
  float dipho_subleadPhi;
  float dipho_subleadEnergy;
  float dipho_subleadR9;
  float dipho_sublead_ptoM;
  float dipho_sublead_sigmaEoE;
  float dipho_subleadIDMVA;
  float dipho_sublead_PSV;
  float dipho_subleadGenMatch;
  float dipho_mva;
  float dipho_mva_training2016_best;
  float dipho_mva_training2017_v1;
  float dipho_mva_training2017_bug;
  
  float MetPt;
  float MetPhi;
  
  float nJets;
  float nJets_bTagLoose;
  float nJets_bTagMedium;
  float nJets_bTagTight;
  
  float jet_pt[nJet];
  float jet_eta[nJet];
  float jet_phi[nJet];
  float jet_bdiscriminant[nJet];
  
  float mu_pt[nLep];
  float mu_eta[nLep];
  float mu_phi[nLep];
  float mu_energy[nLep];
  float mu_IDVector[nLep][3];
  float mu_miniIso[nLep];
  float mu_trackIso[nLep];
  float mu_charge[nLep];
  float mu_sumChargedHadronPt[nLep];
  float mu_sumNeutralHadronEt[nLep];
  float mu_sumPhotonEt[nLep];
  float mu_sumPUPt[nLep];
  float ele_pt[nLep];
  float ele_eta[nLep];
  float ele_phi[nLep];
  float ele_energy[nLep];
  float ele_IDVector[nLep][10];
  float ele_miniIso[nLep];
  float ele_ecalEnergy[nLep];
  float ele_SCx[nLep];
  float ele_SCy[nLep];
  float ele_SCz[nLep];
  float ele_charge[nLep];
  float ele_SCeta[nLep];
  float ele_SCphi[nLep];
  float ele_dEtaTrk[nLep];
  float ele_dPhiTrk[nLep];
  
  float nElectrons;
  float nMuons;
  float lepton_leadPt;
  float lepton_leadEta;
  
  float mu1_pt;
  float mu1_eta;
  float mu1_phi;
  float mu2_pt;
  float mu2_eta;
  float mu2_phi;
  float ele1_pt;
  float ele1_eta;
  float ele1_phi;
  float ele2_pt;
  float ele2_eta;
  float ele2_phi;
  float jet1_pt;
  float jet1_eta;
  float jet1_phi;
  float jet1_bTag;
  float jet2_pt;
  float jet2_eta;
  float jet2_phi;
  float jet2_bTag;
  float jet3_pt;
  float jet3_eta;
  float jet3_phi;
  float jet3_bTag;
  float jet4_pt;
  float jet4_eta;
  float jet4_phi;
  float jet4_bTag;
  float bTag1;
  float bTag2;
  float bTag3;
  float bTag4;
};
  
void InitTreeVars(TChain* chain, TreeVars& treeVars, const bool& useDeepBDisc = true);
void InitOutTree(TTree* tree, TreeVars& treeVars);
void InitOutTreeVars(TreeVars& treeVars);

#endif
