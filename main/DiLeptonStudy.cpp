#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/CMS_lumi.h"
#include "interface/SetTDRStyle.h"
#include "interface/TreeUtils.h"
#include "interface/AnalysisUtils.h"
#include "interface/RooFitUtils.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TMath.h"

#include "RooMsgService.h"
#include "RooRealVar.h"

using namespace std;
using namespace RooFit;

bool verbosity = false;
bool cutBased = true;

float lumi = 41.7;

ControlSampleType csType = kInvertBTag;



int main(int argc, char* argv[])
{
  if( argc < 2 )
  {
    std::cerr << ">>>>> DiLeptonStudy.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return -1;
  }
  
  
  //------------------
  // graphics settings
  
  writeExtraText = true;
  extraText  = "Preliminary";
  lumi_sqrtS = Form("%.1f fb^{-1} (13 TeV)",lumi);
  
  setTDRStyle();
  gStyle -> SetOptFit(0);
  gStyle -> SetOptStat(0);
  
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  
  
  //----------------------
  // parse the config file
  
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  std::vector<std::string> input = opts.GetOpt<std::vector<std::string> >("Input.input");
  
  std::string outputNtupleFolder = opts.GetOpt<std::string>("Output.outputNtupleFolder");
  std::string outputPlotFolder = opts.GetOpt<std::string>("Output.outputPlotFolder");
  system(Form("mkdir %s",outputNtupleFolder.c_str()));
  
  int doDiphoMVAScan = opts.GetOpt<int>("Cuts.doDiphoMVAScan");
  float diphoMVAMin = opts.GetOpt<float>("Cuts.diphoMVAMin");
  
  
  
  //------------------------
  // define roofit variables
  
  RooRealVar mass_("mass_", "m_{#gamma#gamma}", 100, 180);
  RooRealVar weight_("weight_", "weight", -100, 100);
  RooRealVar mva_("mva_", "mva", -1., 1.);  
  
  RooDataSet dataOneCat("dataOneCat","dataOneCat", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet csOneCat("csOneCat","csOneCat", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet tthOneCat("tthOneCat","tthOneCat", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  
  RooDataSet dataDiLepton("dataDiLepton","dataDiLepton", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet csDiLepton("csDiLepton","csDiLepton", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet tthDiLepton("tthDiLepton","tthDiLepton", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  
  RooDataSet dataSingleLepton("dataSingleLepton","dataSingleLepton", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet csSingleLepton("csSingleLepton","csSingleLepton", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet tthSingleLepton("tthSingleLepton","tthSingleLepton", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  
  RooDataSet dataOneCatNoBTag("dataOneCatNoBTag","dataOneCatNoBTag", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet csOneCatNoBTag("csOneCatNoBTag","csOneCatNoBTag", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet tthOneCatNoBTag("tthOneCatNoBTag","tthOneCatNoBTag", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  
  RooDataSet dataDiLeptonNoBTag("dataDiLeptonNoBTag","dataDiLeptonNoBTag", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet csDiLeptonNoBTag("csDiLeptonNoBTag","csDiLeptonNoBTag", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet tthDiLeptonNoBTag("tthDiLeptonNoBTag","tthDiLeptonNoBTag", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  
  RooDataSet dataSingleLeptonNoBTag("dataSingleLeptonNoBTag","dataSingleLeptonNoBTag", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet csSingleLeptonNoBTag("csSingleLeptonNoBTag","csSingleLeptonNoBTag", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  RooDataSet tthSingleLeptonNoBTag("tthSingleLeptonNoBTag","tthSingleLeptonNoBTag", RooArgSet(mass_, weight_, mva_), WeightVar(weight_));
  
  
  //------------------
  // define histograms
  
  std::map<std::string,TH1F*> mass_oneCat_histo;
  std::map<std::string,TH1F*> mass_diMu_histo;
  std::map<std::string,TH1F*> mass_diEle_histo;
  std::map<std::string,TH1F*> mass_Mixed_histo;
  std::map<std::string,TH1F*> mass_diLepton_histo;
  std::map<std::string,TH1F*> mass_singleLepton_histo;
  
  
  //----------
  // get trees
  
  std::map<std::string,TChain*> trees;
  std::map<std::string,int> types; // -1 = data; -2 = control sample; 1 = MC bkg; 2 = MC signal;

  for(unsigned int n = 0; n < input.size()/4; ++n)
  {
    std::string inFileName = input.at(0+n*4);
    std::string treeName   = input.at(1+n*4);
    std::string label      = input.at(2+n*4);
    int type               = atoi(input.at(3+n*4).c_str());
    
    if( trees[label] == 0 )
    {
      trees[label] = new TChain(Form("tree_%s",label.c_str()),"");
      
      mass_oneCat_histo[label]    = new TH1F(("mass_oneCat_histo"   +label).c_str(), "; m^{#gamma#gamma} (GeV); Counts", 80, 100, 180 );
      mass_diMu_histo[label]      = new TH1F(("mass_diMu_histo"     +label).c_str(), "; m^{#gamma#gamma} (GeV); Counts", 80, 100, 180 );
      mass_diEle_histo[label]     = new TH1F(("mass_diEle_histo"    +label).c_str(), "; m^{#gamma#gamma} (GeV); Counts", 80, 100, 180 );
      mass_Mixed_histo[label]     = new TH1F(("mass_Mixed_histo"    +label).c_str(), "; m^{#gamma#gamma} (GeV); Counts", 80, 100, 180 );
      mass_diLepton_histo[label]  = new TH1F(("mass_diLepton_histo" +label).c_str(), "; m^{#gamma#gamma} (GeV); Counts", 80, 100, 180 );
      mass_singleLepton_histo[label] = new TH1F(("mass_singleLepton_histo"+label).c_str(), "; m^{#gamma#gamma} (GeV); Counts", 80, 100, 180 );
    }
    
    std::cout << ">>> Adding trees " << inFileName+"/"+treeName << " to chain " << "tree_"+label << std::endl;
    trees[label] -> Add((inFileName+"/"+treeName).c_str());
    types[label] = type;
  }
  
  
  //---------------
  // tree variables
  TreeVars treeVars;
  
  
  
  //------------------
  // loop over samples
  std::map<std::string,std::map<float,float> > nEvents_cutBased;
  std::map<std::string,std::map<float,float> > nEvents_mvaCut;
  std::map<std::string,std::map<float,float> > nEvents_mvaCut_new;
  for(std::map<std::string,TChain*>::const_iterator treeIt = trees.begin(); treeIt != trees.end(); ++treeIt)
  {
    std::string label = treeIt -> first;
    TChain* tree = treeIt -> second;
    int type = types[label];
    
    InitTreeVars(tree,treeVars);
    float lumiFactor = type > 0 ? lumi : 1.;
    
    TFile* outFile = TFile::Open(Form("%s/plotTree_%s.root",outputNtupleFolder.c_str(),label.c_str()),"RECREATE");
    outFile -> cd();
    TTree* outTree_1jet = new TTree("plotTree_1jet","plotTree_1jet");
    TTree* outTree_diLepton = new TTree("plotTree_diLepton","plotTree_diLepton");
    TTree* outTree_singleLepton = new TTree("plotTree_singleLepton","plotTree_singleLepton");
    TTree* outTree_oneCategory = new TTree("plotTree_oneCategory","plotTree_oneCategory");
    TTree* outTree_diLepton_noBTag = new TTree("plotTree_diLepton_noBTag","plotTree_diLepton_noBTag");
    TTree* outTree_singleLepton_noBTag = new TTree("plotTree_singleLepton_noBTag","plotTree_singleLepton_noBTag");
    TTree* outTree_oneCategory_noBTag = new TTree("plotTree_oneCategory_noBTag","plotTree_oneCategory_noBTag");
    InitOutTreeVars(outTree_1jet,treeVars);
    InitOutTreeVars(outTree_diLepton,treeVars);
    InitOutTreeVars(outTree_singleLepton,treeVars);
    InitOutTreeVars(outTree_oneCategory,treeVars);
    InitOutTreeVars(outTree_diLepton_noBTag,treeVars);
    InitOutTreeVars(outTree_singleLepton_noBTag,treeVars);
    InitOutTreeVars(outTree_oneCategory_noBTag,treeVars);
    
    
    int nEntries = tree->GetEntries();
    int nEntries_commonCuts = 0;
    int nEntries_diphoMVACut = 0;
    int nEntries_diLepton = 0;
    int nEntries_singleLepton = 0;
    std::cout << std::endl;
    for(int i=0; i<nEntries; i++)
    {
      tree -> GetEntry(i);
      // if( i >= 10 ) break;
      if( i%1000==0 ) std::cout << "Processing tag " << label << ", event " << i << " out of " << nEntries << "\r" << std::flush;
      
      // common cuts
      if( treeVars.dipho_mass < 100 || treeVars.dipho_mass > 180 ) continue;
      if( type == -1 && treeVars.dipho_mass > 115 && treeVars.dipho_mass < 135 ) continue;
      
      int nJets = 0;
      for(int jIndex=0; jIndex<6; jIndex++)
        if(treeVars.jet_pt[jIndex]>25. && abs(treeVars.jet_eta[jIndex])<2.4)
          ++nJets;
      
      if( nJets < 1 ) continue;
      
      ++nEntries_commonCuts;
      
      outTree_1jet -> Fill();
      
      bool passCutBased = CutBasedSelection(treeVars,0.33,0.25,0.0,0.0,3.15,3.);
      
      // fill event counters - cut based
        if( type == 1 )
          nEvents_cutBased["bkg"][0.] += treeVars.weight;
        if( type == 2 && label == "ttH" )
          nEvents_cutBased["sig"][0.] += treeVars.weight;
      if( passCutBased )
      {
        if( type == 1 )
          nEvents_cutBased["bkg"][1.] += treeVars.weight;
        if( type == 2 && label == "ttH" )
          nEvents_cutBased["sig"][1.] += treeVars.weight;
      }
      
      // fill event counters - mva
      for(float mvaCut = -1.; mvaCut < 1.; mvaCut+=0.025)
      {
        if( treeVars.dipho_mva > mvaCut )
        {
          if( type == 1 )
            nEvents_mvaCut["bkg"][mvaCut] += treeVars.weight;
          if( type == 2 && label == "ttH" )
            nEvents_mvaCut["sig"][mvaCut] += treeVars.weight;
        }
        
        if( treeVars.dipho_mva_training2017_v1 > mvaCut )
        {
          if( type == 1 )
            nEvents_mvaCut_new["bkg"][mvaCut] += treeVars.weight;
          if( type == 2 && label == "ttH" )
            nEvents_mvaCut_new["sig"][mvaCut] += treeVars.weight;
        }
      }
      
      if( treeVars.dipho_mva < diphoMVAMin ) continue;
      
      ++nEntries_diphoMVACut;
            
      if( cutBased && !passCutBased ) continue;
      
      
      // if( treeVars.event == 324559 )
      // {
      //   std::cout << "################################" << std::endl;
        
      //   std::cout << ">>> photons: " << std::endl;
        
      //   std::cout << ">>>>>> lead:      pt: " << treeVars.dipho_leadPt << "   eta: " << treeVars.dipho_leadEta << "   phi: " << treeVars.dipho_leadPhi << std::endl;
      //   std::cout << ">>>>>> sublead:   pt: " << treeVars.dipho_subleadPt << "   eta: " << treeVars.dipho_subleadEta << "   phi: " << treeVars.dipho_subleadPhi << std::endl;
        
      //   std::cout << ">>> muons: " << std::endl;
      //   for(int jj = 0; jj < 2; ++jj)
      //   {
      //     std::cout << ">>>>>> jj: " << jj << "   pt: " << treeVars.mu_pt[jj] << "   eta: " << treeVars.mu_eta[jj] << "   phi: " << treeVars.mu_phi[jj];
      //     for(int kk = 0; kk < 2; ++kk)
      //       std::cout << "   IDVector[" << kk << "]: " << treeVars.mu_IDVector[jj][kk];
      //     std::cout << std::endl;
      //   }
        
      //   std::cout << ">>> electrons: " << std::endl;
      //   for(int jj = 0; jj < 2; ++jj)
      //   {
      //     std::cout << ">>>>>> jj: " << jj << "   pt: " << treeVars.ele_pt[jj] << "   eta: " << treeVars.ele_eta[jj] << "   phi: " << treeVars.ele_phi[jj];
      //     for(int kk = 0; kk < 10; ++kk)
      //       std::cout << "   IDVector[" << kk << "]: " << treeVars.ele_IDVector[jj][kk];
      //     std::cout << std::endl;
      //   }
      //   std::cout << "################################" << std::endl;
      // }
      
      //-------------
      // one category
      bool bTagSelection = true;
      if( OneCategorySelection(treeVars,type,bTagSelection,csType) == true )
      {
        mass_oneCat_histo[label] -> Fill(treeVars.dipho_mass, treeVars.weight*lumiFactor);
        
        mass_ = treeVars.dipho_mass;
        mva_ = treeVars.dipho_mva;
        if( label == "ttH"            ) tthOneCat.add(RooArgSet(mass_,mva_),  treeVars.weight*lumiFactor);
        if( label == "CS_oneCategory" ) csOneCat.add(RooArgSet(mass_,mva_),   treeVars.weight*lumiFactor);
        if( label == "data"           ) dataOneCat.add(RooArgSet(mass_,mva_), treeVars.weight*lumiFactor);
        
        outTree_oneCategory -> Fill();
      } // one category
      
      
      //--------------------------
      // two categories - dilepton
      
      DiLeptonCategories cat = None;
      if( DiLeptonSelection(treeVars,type,bTagSelection,cat,csType,false) )
      {
        ++nEntries_diLepton;
        
        mass_diLepton_histo[label] -> Fill(treeVars.dipho_mass, treeVars.weight*lumiFactor);
        
        mass_ = treeVars.dipho_mass;
        mva_ = treeVars.dipho_mva;
        if( label == "ttH"         ) tthDiLepton.add(RooArgSet(mass_,mva_),  treeVars.weight*lumiFactor);
        if( label == "CS_diLepton" ) csDiLepton.add(RooArgSet(mass_,mva_),   treeVars.weight*lumiFactor);
        if( label == "data"        ) dataDiLepton.add(RooArgSet(mass_,mva_), treeVars.weight*lumiFactor);
        
        outTree_diLepton -> Fill();
        
        if( verbosity )
        {
          std::cout << ">>> DiLeptonSelection::pass   run " << treeVars.run << "   event " << treeVars.event << std::endl;
          DiLeptonSelection(treeVars,type,bTagSelection,cat,csType,verbosity);
          std::cout << ">>> DiLeptonSelection::selected category: " << cat << std::endl;
          std::cout << std::endl;
        }
        if( type != -2 ) continue;
      } // two categories - dilepton
      else
      {
        if( verbosity )
        {
          std::cout << ">>> DiLeptonSelection::fail   run " << treeVars.run << "   event " << treeVars.event << std::endl;
          DiLeptonSelection(treeVars,type,bTagSelection,cat,csType,verbosity);
          std::cout << std::endl;
        }
      }
      
      
      //-------------------------------
      // two categories - single lepton
      
      if( SingleLeptonSelection(treeVars,type,bTagSelection,csType) )
      {
        ++nEntries_singleLepton;
        
        mass_singleLepton_histo[label] -> Fill(treeVars.dipho_mass, treeVars.weight*lumiFactor);
        
        mass_ = treeVars.dipho_mass;
        mva_ = treeVars.dipho_mva;
        if( label == "ttH"             ) tthSingleLepton.add(RooArgSet(mass_,mva_),  treeVars.weight*lumiFactor);
        if( label == "CS_singleLepton" ) csSingleLepton.add(RooArgSet(mass_,mva_),   treeVars.weight*lumiFactor);
        if( label == "data"            ) dataSingleLepton.add(RooArgSet(mass_,mva_), treeVars.weight*lumiFactor);
        
        outTree_singleLepton -> Fill();
        
        if( verbosity )
        {
          std::cout << ">>> SingleLeptonSelection::pass   run " << treeVars.run << "   event " << treeVars.event << std::endl;
        }
      } // two categories - single lepton
      // else
      // {
      //   SingleLeptonSelection(treeVars,type,bTagSelection,csType,verbosity);
      // }
      
      
    } // loop over events
    
    outTree_1jet -> AutoSave();
    outTree_diLepton -> AutoSave();
    outTree_singleLepton -> AutoSave();
    outTree_oneCategory -> AutoSave();
    outTree_diLepton_noBTag -> AutoSave();
    outTree_singleLepton_noBTag -> AutoSave();
    outTree_oneCategory_noBTag -> AutoSave();
    outFile -> Close();
    
    std::cout << "\nProcessed tag " << label << ", " << nEntries << " events out of " << nEntries << std::endl;
    std::cout << ">>> nEntries_commonCuts: "   << nEntries_commonCuts   << std::endl;
    std::cout << ">>> nEntries_diphoMVACut: "  << nEntries_diphoMVACut  << std::endl;
    std::cout << ">>> nEntries_diLepton: "     << nEntries_diLepton     << std::endl;
    std::cout << ">>> nEntries_singleLepton: " << nEntries_singleLepton << std::endl;
  }
  
  
  // MakePlot2(mass_oneCat_histo, "oneCategory");
  // MakePlot2(mass_diMu_histo, "diMuon");
  // MakePlot2(mass_diEle_histo, "diElectron");
  // MakePlot2(mass_Mixed_histo, "mixed");
  MakePlot2(mass_diLepton_histo, "diLepton");
  MakePlot2(mass_singleLepton_histo, "singleLepton");
  
  system(Form("mkdir %s",outputPlotFolder.c_str()));
  system(Form("cp /afs/cern.ch/user/a/abenagli/www/index.php %s",outputPlotFolder.c_str()));
  system(Form("mv *.png %s",outputPlotFolder.c_str()));
  system(Form("mv *.pdf %s",outputPlotFolder.c_str()));
  
  
  
  std::cout << "#########################################################################" << std::endl;
  std::cout << "################################## ROC ##################################" << std::endl;
  
  TFile* outFile_global = new TFile("DiLeptonStudy.root","RECREATE");
  outFile_global -> cd();
  
  TGraph* g_ROC_cutBased = new TGraph();
  g_ROC_cutBased -> SetPoint(g_ROC_cutBased->GetN(),nEvents_cutBased["sig"][1.]/nEvents_cutBased["sig"][0.],1.-nEvents_cutBased["bkg"][1.]/nEvents_cutBased["bkg"][0.]);
  
  TGraph* g_ROC_mva = new TGraph();
  TGraph* g_ROC_mva_new = new TGraph();
  for(float mvaCut = -1.; mvaCut < 1.; mvaCut+=0.025)
  {
    g_ROC_mva -> SetPoint(g_ROC_mva->GetN(),nEvents_mvaCut["sig"][mvaCut]/nEvents_mvaCut["sig"][-1.],1.-nEvents_mvaCut["bkg"][mvaCut]/nEvents_mvaCut["bkg"][-1.]);
    g_ROC_mva_new -> SetPoint(g_ROC_mva_new->GetN(),nEvents_mvaCut_new["sig"][mvaCut]/nEvents_mvaCut_new["sig"][-1.],1.-nEvents_mvaCut_new["bkg"][mvaCut]/nEvents_mvaCut_new["bkg"][-1.]);
  }
  
  outFile_global -> cd();
  
  g_ROC_cutBased -> Write("g_ROC_cutBased");
  g_ROC_mva -> Write("g_ROC_mva");
  g_ROC_mva_new -> Write("g_ROC_mva_new");
  
  std::cout << "#########################################################################" << std::endl;
  
  
  
  bool doFit = 1;
  bool doSimultaneous = 0;
  
  if( doFit )
  {
    std::cout << "#################################################################################" << std::endl;
    std::cout << "################################## FIT RESULTS ##################################" << std::endl;
    
    std::vector<float> diLeptonSignificance_cs = makeFits(&tthDiLepton, &csDiLepton, -1., "diLepton_cs", 0, 0);
    std::vector<float> singleLeptonSignificance_cs = makeFits(&tthSingleLepton, &csSingleLepton, -1., "singleLepton_cs", 0, 0);
    std::vector<float> diLeptonSignificance_sb = makeFits(&tthDiLepton, &dataDiLepton, -1., "diLepton_sb", 0, 1, true);
    std::vector<float> singleLeptonSignificance_sb = makeFits(&tthSingleLepton, &dataSingleLepton, -1., "singleLepton_sb", 0, 1, true);
    
    float purityDiLepton = (mass_diLepton_histo["ttH"] -> Integral() /
                            (mass_diLepton_histo["ttH"]->Integral() + mass_diLepton_histo["ggH"]->Integral() + mass_diLepton_histo["VBF"]->Integral() + mass_diLepton_histo["VH"]->Integral() ) ) * 100.;
    float puritySingleLepton = (mass_singleLepton_histo["ttH"] -> Integral() /
                                (mass_singleLepton_histo["ttH"]->Integral() + mass_singleLepton_histo["ggH"]->Integral() + mass_singleLepton_histo["VBF"]->Integral() + mass_singleLepton_histo["VH"]->Integral() ) ) * 100.;
    
    std::cout << "Significance with two leptonic categories (control sample): " << sqrt(diLeptonSignificance_cs[2]*diLeptonSignificance_cs[2] + singleLeptonSignificance_cs[2]*singleLeptonSignificance_cs[2]) << " (" << diLeptonSignificance_cs[2] <<  " from dileptonic category, " << singleLeptonSignificance_cs[2] << " from one lepton category)" << std::endl << std::endl;
    std::cout << "Significance with two leptonic categories (data sideband): " << sqrt(diLeptonSignificance_sb[2]*diLeptonSignificance_sb[2] + singleLeptonSignificance_sb[2]*singleLeptonSignificance_sb[2]) << " (" << diLeptonSignificance_sb[2] <<  " from dileptonic category, " << singleLeptonSignificance_sb[2] << " from one lepton category)" << std::endl << std::endl;
    
    std::cout << "ttH events diLepton:        " <<  mass_diLepton_histo["ttH"] -> Integral() << " events, tag purity: " << purityDiLepton     << std::endl; 
    std::cout << "ttH events one Lepton:      " << mass_singleLepton_histo["ttH"] -> Integral() << " events, tag purity: " << puritySingleLepton << std::endl; 
    std::cout << "#################################################################################" << std::endl;    
    
    
    if( doSimultaneous )
    {
      float significanceCombined = makeFitSimulataneous(&tthDiLepton, &csDiLepton, &tthSingleLepton, &csSingleLepton, "ChiLoSa");
      std::cout << std::endl << "Significance from the simultanoeus fit of the two categories " << significanceCombined << std::endl;
    }
    
    
    if( doDiphoMVAScan )
    {
      std::cout << std::endl;
      std::cout << "####################################################################################" << std::endl;
      std::cout << "################################## DOING MVA SCAN ##################################" << std::endl;
      
      float sigN_diLepton_initial = diLeptonSignificance_cs[0];
      float bkgN_diLepton_initial = diLeptonSignificance_cs[1];
      float sigN_singleLepton_initial = singleLeptonSignificance_cs[0];
      float bkgN_singleLepton_initial = singleLeptonSignificance_cs[1];
      
      TGraph* g_sigN_diLepton = new TGraph();
      TGraph* g_bkgN_diLepton = new TGraph();
      TGraph* g_sigN_singleLepton = new TGraph();
      TGraph* g_bkgN_singleLepton = new TGraph();
      
      TGraph* g_sigEff_diLepton = new TGraph();
      TGraph* g_bkgEff_diLepton = new TGraph();
      TGraph* g_sigEff_singleLepton = new TGraph();
      TGraph* g_bkgEff_singleLepton = new TGraph();
      
      TGraph* g_diLeptonSignificance_cs = new TGraph();
      TGraph* g_singleLeptonSignificance_cs = new TGraph();
      TGraph* g_diLeptonSignificance_sb = new TGraph();
      TGraph* g_singleLeptonSignificance_sb = new TGraph();
      
      TGraph* g_diLeptonNaiveSignificance_cs = new TGraph();
      TGraph* g_singleLeptonNaiveSignificance_cs = new TGraph();
      TGraph* g_diLeptonNaiveSignificance_sb = new TGraph();
      TGraph* g_singleLeptonNaiveSignificance_sb = new TGraph();
      
      for(float mvaMin = -1.; mvaMin < 1.; mvaMin+=0.05)
      {
        diLeptonSignificance_cs = makeFits(&tthDiLepton, &csDiLepton, mvaMin, "diLepton_cs", 0, 0);
        singleLeptonSignificance_cs = makeFits(&tthSingleLepton, &csSingleLepton, mvaMin, "singleLepton_cs", 0, 0);
        diLeptonSignificance_sb = makeFits(&tthDiLepton, &dataDiLepton, mvaMin, "diLepton_sb", 0, 1, true);
        singleLeptonSignificance_sb = makeFits(&tthSingleLepton, &dataSingleLepton, mvaMin, "singleLepton_sb", 0, 1, true);
        
        g_sigN_diLepton     -> SetPoint(g_sigN_diLepton->GetN(),    mvaMin,    diLeptonSignificance_cs[0]);
        g_bkgN_diLepton     -> SetPoint(g_bkgN_diLepton->GetN(),    mvaMin,    diLeptonSignificance_cs[1]/4.);
        g_sigN_singleLepton -> SetPoint(g_sigN_singleLepton->GetN(),mvaMin,singleLeptonSignificance_cs[0]);
        g_bkgN_singleLepton -> SetPoint(g_bkgN_singleLepton->GetN(),mvaMin,singleLeptonSignificance_cs[1]/4.);
        
        g_sigEff_diLepton     -> SetPoint(g_sigEff_diLepton->GetN(),    mvaMin,    diLeptonSignificance_cs[0]/    sigN_diLepton_initial);
        g_bkgEff_diLepton     -> SetPoint(g_bkgEff_diLepton->GetN(),    mvaMin,    diLeptonSignificance_cs[1]/    bkgN_diLepton_initial);
        g_sigEff_singleLepton -> SetPoint(g_sigEff_singleLepton->GetN(),mvaMin,singleLeptonSignificance_cs[0]/sigN_singleLepton_initial);
        g_bkgEff_singleLepton -> SetPoint(g_bkgEff_singleLepton->GetN(),mvaMin,singleLeptonSignificance_cs[1]/bkgN_singleLepton_initial);
        
        g_diLeptonSignificance_cs     -> SetPoint(    g_diLeptonSignificance_cs->GetN(),mvaMin,    diLeptonSignificance_cs[2]);
        g_singleLeptonSignificance_cs -> SetPoint(g_singleLeptonSignificance_cs->GetN(),mvaMin,singleLeptonSignificance_cs[2]);
        g_diLeptonSignificance_sb     -> SetPoint(    g_diLeptonSignificance_sb->GetN(),mvaMin,    diLeptonSignificance_sb[2]);
        g_singleLeptonSignificance_sb -> SetPoint(g_singleLeptonSignificance_sb->GetN(),mvaMin,singleLeptonSignificance_sb[2]);
        
        g_diLeptonNaiveSignificance_cs     -> SetPoint(    g_diLeptonNaiveSignificance_cs->GetN(),mvaMin,    diLeptonSignificance_cs[0]/sqrt(    diLeptonSignificance_cs[1]));
        g_singleLeptonNaiveSignificance_cs -> SetPoint(g_singleLeptonNaiveSignificance_cs->GetN(),mvaMin,singleLeptonSignificance_cs[0]/sqrt(singleLeptonSignificance_cs[1]));
        g_diLeptonNaiveSignificance_sb     -> SetPoint(    g_diLeptonNaiveSignificance_sb->GetN(),mvaMin,    diLeptonSignificance_sb[0]/sqrt(    diLeptonSignificance_sb[1]));
        g_singleLeptonNaiveSignificance_sb -> SetPoint(g_singleLeptonNaiveSignificance_sb->GetN(),mvaMin,singleLeptonSignificance_sb[0]/sqrt(singleLeptonSignificance_sb[1]));
      }
      
      
      TCanvas* c = new TCanvas();
      TH1F* hPad;
  
      c = new TCanvas();
      c -> cd();
      gPad -> SetGridx();
      gPad -> SetGridy();
      
      hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,1.,1.5) );
      hPad -> SetTitle(";diphoton MVA cut;significance");
      
      g_diLeptonSignificance_cs -> SetLineColor(kRed);
      g_diLeptonSignificance_cs -> SetLineStyle(1);
      g_diLeptonSignificance_cs -> SetLineWidth(2);
      g_diLeptonSignificance_cs -> Draw("L,same");
      
      g_diLeptonSignificance_sb -> SetLineColor(kBlack);
      g_diLeptonSignificance_sb -> SetLineStyle(1);
      g_diLeptonSignificance_sb -> SetLineWidth(2);
      g_diLeptonSignificance_sb -> Draw("L,same");
      
      g_diLeptonNaiveSignificance_cs -> SetLineColor(kRed);
      g_diLeptonNaiveSignificance_cs -> SetLineStyle(2);
      g_diLeptonNaiveSignificance_cs -> SetLineWidth(2);
      g_diLeptonNaiveSignificance_cs -> Draw("L,same");
      
      g_diLeptonNaiveSignificance_sb -> SetLineColor(kBlack);
      g_diLeptonNaiveSignificance_sb -> SetLineStyle(2);
      g_diLeptonNaiveSignificance_sb -> SetLineWidth(2);
      g_diLeptonNaiveSignificance_sb -> Draw("L,same");
      
      c -> SaveAs("c_significance_vs_mva_diLepton.png");
      c -> SaveAs("c_significance_vs_mva_diLepton.pdf");
      
      
      c = new TCanvas();
      c -> cd();
      gPad -> SetGridx();
      gPad -> SetGridy();
      
      hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,1.,1.5) );
      hPad -> SetTitle(";diphoton MVA cut;significance");
      
      g_singleLeptonSignificance_cs -> SetLineColor(kRed);
      g_singleLeptonSignificance_cs -> SetLineStyle(1);
      g_singleLeptonSignificance_cs -> SetLineWidth(2);
      g_singleLeptonSignificance_cs -> Draw("L,same");
      
      g_singleLeptonSignificance_sb -> SetLineColor(kBlack);
      g_singleLeptonSignificance_sb -> SetLineStyle(1);
      g_singleLeptonSignificance_sb -> SetLineWidth(2);
      g_singleLeptonSignificance_sb -> Draw("L,same");
      
      g_singleLeptonNaiveSignificance_cs -> SetLineColor(kRed);
      g_singleLeptonNaiveSignificance_cs -> SetLineStyle(2);
      g_singleLeptonNaiveSignificance_cs -> SetLineWidth(2);
      g_singleLeptonNaiveSignificance_cs -> Draw("L,same");
      
      g_singleLeptonNaiveSignificance_sb -> SetLineColor(kBlack);
      g_singleLeptonNaiveSignificance_sb -> SetLineStyle(2);
      g_singleLeptonNaiveSignificance_sb -> SetLineWidth(2);
      g_singleLeptonNaiveSignificance_sb -> Draw("L,same");
      
      c -> SaveAs("c_significance_vs_mva_singleLepton.png");
      c -> SaveAs("c_significance_vs_mva_singleLepton.pdf");
      
      
      c = new TCanvas();
      c -> cd();
      gPad -> SetGridx();
      gPad -> SetGridy();
      
      hPad = (TH1F*)( gPad->DrawFrame(-1.,0.01,1.,5.) );
      hPad -> SetTitle(";diphoton MVA cut;number of events");
      
      g_sigN_diLepton -> SetLineColor(kRed);
      g_sigN_diLepton -> SetLineWidth(2);
      g_sigN_diLepton -> Draw("L,same");
      g_bkgN_diLepton -> SetLineColor(kBlack);
      g_bkgN_diLepton -> SetLineWidth(2);
      g_bkgN_diLepton -> Draw("L,same");
      
      outFile_global -> cd();
      
      c -> SaveAs("c_N_vs_mva_diLepton.png");
      c -> SaveAs("c_N_vs_mva_diLepton.pdf");
      
      
      c = new TCanvas();
      c -> cd();
      gPad -> SetGridx();
      gPad -> SetGridy();
      
      hPad = (TH1F*)( gPad->DrawFrame(-1.,0.01,1.,5.) );
      hPad -> SetTitle(";diphoton MVA cut;number of events");
      
      g_sigN_singleLepton -> SetLineColor(kRed);
      g_sigN_singleLepton -> SetLineWidth(2);
      g_sigN_singleLepton -> Draw("L,same");
      g_bkgN_singleLepton -> SetLineColor(kBlack);
      g_bkgN_singleLepton -> SetLineWidth(2);
      g_bkgN_singleLepton -> Draw("L,same");
      
      outFile_global -> cd();
      
      c -> SaveAs("c_N_vs_mva_singleLepton.png");
      c -> SaveAs("c_N_vs_mva_singleLepton.pdf");
      
      
      c = new TCanvas();
      c -> cd();
      c -> SetLogy();
      gPad -> SetGridx();
      gPad -> SetGridy();
      
      hPad = (TH1F*)( gPad->DrawFrame(-1.,0.01,1.,1.) );
      hPad -> SetTitle(";diphoton MVA cut;efficiency");
      
      g_sigEff_diLepton -> SetLineColor(kRed);
      g_sigEff_diLepton -> SetLineWidth(2);
      g_sigEff_diLepton -> Draw("L,same");
      g_bkgEff_diLepton -> SetLineColor(kBlack);
      g_bkgEff_diLepton -> SetLineWidth(2);
      g_bkgEff_diLepton -> Draw("L,same");
      
      outFile_global -> cd();
      
      c -> SaveAs("c_eff_vs_mva_diLepton.png");
      c -> SaveAs("c_eff_vs_mva_diLepton.pdf");
      
      
      c = new TCanvas();
      c -> cd();
      c -> SetLogy();
      gPad -> SetGridx();
      gPad -> SetGridy();
      
      hPad = (TH1F*)( gPad->DrawFrame(-1.,0.01,1.,1.) );
      hPad -> SetTitle(";diphoton MVA cut;efficiency");
      
      g_sigEff_singleLepton -> SetLineColor(kRed);
      g_sigEff_singleLepton -> SetLineWidth(2);
      g_sigEff_singleLepton -> Draw("L,same");
      g_bkgEff_singleLepton -> SetLineColor(kBlack);
      g_bkgEff_singleLepton -> SetLineWidth(2);
      g_bkgEff_singleLepton -> Draw("L,same");
      
      outFile_global -> cd();
      
      c -> SaveAs("c_eff_vs_mva_singleLepton.png");
      c -> SaveAs("c_eff_vs_mva_singleLepton.pdf");
      
      
      outFile_global -> cd();
      g_diLeptonSignificance_cs -> Write("g_diLeptonSignificance_cs");
      g_singleLeptonSignificance_cs -> Write("g_singleLeptonSignificance_cs");
      g_diLeptonSignificance_sb -> Write("g_diLeptonSignificance_sb");
      g_singleLeptonSignificance_sb -> Write("g_singleLeptonSignificance_sb");
      g_sigEff_diLepton -> Write("g_sigEff_diLepton");
      g_bkgEff_diLepton -> Write("g_bkgEff_diLepton");
      g_sigEff_singleLepton -> Write("g_sigEff_singleLepton");
      g_bkgEff_singleLepton -> Write("g_bkgEff_singleLepton");

    }
    
    system(Form("mkdir %s/Optimization",outputPlotFolder.c_str()));
    system(Form("cp /afs/cern.ch/user/a/abenagli/www/index.php %s/Optimization/",outputPlotFolder.c_str()));
    system(Form("mv *.png %s/Optimization/",outputPlotFolder.c_str()));
    system(Form("mv *.pdf %s/Optimization/",outputPlotFolder.c_str()));
  }
  
  outFile_global -> Close();
}

//  LocalWords:  endl
