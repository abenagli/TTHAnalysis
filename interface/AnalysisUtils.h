#ifndef ANALYSIS_UTILS_H
#define ANALYSIS_UTILS_H

#include "interface/CMS_lumi.h"
#include "interface/TreeUtils.h"

#include <iostream>
#include <iomanip>
#include <map>

#include "TLorentzVector.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

#define MZ 91.187
#define PI 3.14159265359

#define bDiscriminantThresholdLooseCSV 0.5426
#define bDiscriminantThresholdMediumCSV 0.8484
#define bDiscriminantThresholdTightCSV 0.9535

#define bDiscriminantThresholdLooseDeep  0.1522
#define bDiscriminantThresholdMediumDeep 0.4941
#define bDiscriminantThresholdTightDeep  0.8001


#define oneCatMuID   1 // 0=Loose, 1=Medium, 2=tight
#define singleMuID   1 // 0=Loose, 1=Medium, 2=tight
#define diMuID       1 // 0=Loose, 1=Medium, 2=tight
#define mixedMuID    1 // 0=Loose, 1=Medium, 2=tight
#define oneCatEleID  5 // 0=Veto,  1=Loose, 2=Medium, 3= Tight, 4=MVALoose, 5=MVAMedium, 6=MVATight, 7=MVALooseNoIso, 8=MVAMediumNoIso, 9=MVATightNoIso
#define singleEleID  5 // 0=Veto,  1=Loose, 2=Medium, 3= Tight, 4=MVALoose, 5=MVAMedium, 6=MVATight, 7=MVALooseNoIso, 8=MVAMediumNoIso, 9=MVATightNoIso
#define diEleID      4 // 0=Veto,  1=Loose, 2=Medium, 3= Tight, 4=MVALoose, 5=MVAMedium, 6=MVATight, 7=MVALooseNoIso, 8=MVAMediumNoIso, 9=MVATightNoIso
#define mixedEleID   4 // 0=Veto,  1=Loose, 2=Medium, 3= Tight, 4=MVALoose, 5=MVAMedium, 6=MVATight, 7=MVALooseNoIso, 8=MVAMediumNoIso, 9=MVATightNoIso

enum DiLeptonCategories { None, DiMuon, DiElectron, Mixed };

enum ControlSampleType { kNull, kInvertBTag, kInvertLeptons };



float DeltaEta(const float& eta1, const float& eta2);

float DeltaPhi(const float& phi1, const float& phi2);

float DeltaR(const float& eta1, const float& phi1,
             const float& eta2, const float& phi2);


bool DiMuSelections(TLorentzVector mu1, TLorentzVector mu2, float charge1, float charge2, TLorentzVector ph1, TLorentzVector ph2, float Iso1, float Iso2);
bool DiEleSelections(TLorentzVector ele1, TLorentzVector ele2, float charge1, float charge2, TLorentzVector ph1, TLorentzVector ph2, float iso1, float iso2, float dTrk1, float dTrk2);
bool MixedSelections(TLorentzVector mu, TLorentzVector ele, float charge1, float charge2, TLorentzVector ph1, TLorentzVector ph2);
bool SingleMuSelections(TLorentzVector mu1,  TLorentzVector ph1, TLorentzVector ph2, float iso);
bool SingleEleSelections(TLorentzVector ele1, TLorentzVector ph1, TLorentzVector ph2, float iso, float drTrk);

bool OneCategorySelection(const TreeVars& treeVars, const int& type, const bool& bTagSelection, const ControlSampleType& csType = kNull);
bool DiLeptonSelection(const TreeVars& treeVars, const int& type, const bool& bTagSelection, DiLeptonCategories& cat, const ControlSampleType& csType = kNull, const bool& verbosity = false);
bool SingleLeptonSelection(const TreeVars& treeVars, const int& type,
                           const int& nJetMin, const float& leadJetPtMin, const int& nBTagLooseJetMin, const int& nBTagMediumJetMin, const int& nBTagTightJetMin, 
                           const ControlSampleType& csType = kNull, const bool& verbosity = false,
                           std::vector<int>* goodMu = NULL, std::vector<int>* goodEle = NULL, std::vector<int>* goodJet = NULL);

bool CutBasedSelection(const TreeVars& treeVars,
                       const float& min_lead_ptoM, const float& min_sublead_ptoM,
                       const float& min_leadIDMVA, const float& min_subleadIDMVA,
                       const float& max_deltaphi, const float& max_deltaeta);

void MakePlot(TH1F**, TString title);

void MakePlot2(std::map<std::string,TH1F*>& histos, TString title);

#endif
