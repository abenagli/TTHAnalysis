#ifndef TREE_UTILS_H
#define TREE_UTILS_H

#include <iostream>
#include <vector>
#include <map>

#include "TChain.h"



/*** tree variables ***/
struct TreeVars
{
  unsigned int run;
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
  float dipho_subleadPt;
  float dipho_subleadEta;
  float dipho_subleadPhi;
  float dipho_subleadEnergy;
  float dipho_subleadR9;
  float dipho_sublead_ptoM;
  float dipho_sublead_sigmaEoE;
  float dipho_subleadIDMVA;
  float dipho_mva;
  float dipho_mva_training2016_best;
  float dipho_mva_training2017_v1;
  float dipho_mva_training2017_bug;
  
  float MetPt;
  float MetPhi;
  float ttHMVA;
  
  float nJets;
  float nJets_bTagLoose;
  float nJets_bTagMedium;
  float nJets_bTagTight;
  
  float jet_pt[9];
  float jet_eta[9];
  float jet_phi[9];
  float jet_bdiscriminant[9];
  
  float mu_pt[6];
  float mu_eta[6];
  float mu_phi[6];
  float mu_energy[6];
  float mu_IDVector[6][2];
  float mu_miniIso[6];
  float mu_trackIso[6];
  float mu_charge[6];
  float mu_sumChargedHadronPt[6];
  float mu_sumNeutralHadronEt[6];
  float mu_sumPhotonEt[6];
  float mu_sumPUPt[6];
  float ele_pt[6];
  float ele_eta[6];
  float ele_phi[6];
  float ele_energy[6];
  float ele_IDVector[6][10];
  float ele_miniIso[6];
  float ele_ecalEnergy[6];
  float ele_SCx[6];
  float ele_SCy[6];
  float ele_SCz[6];
  float ele_charge[6];
  float ele_SCeta[6];
  float ele_SCphi[6];
  float ele_dEtaTrk[6];
  float ele_dPhiTrk[6];
};
  
void InitTreeVars(TChain* chain, TreeVars& treeVars);
void InitOutTreeVars(TTree* tree, TreeVars& treeVars);

#endif
