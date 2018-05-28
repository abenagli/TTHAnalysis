#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/TreeUtils.h"
#include "interface/AnalysisUtils.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h" 
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TLegend.h"



int main(int argc, char** argv)
{
  if( argc < 2 )
  {
    std::cerr << ">>>>> createControlSample.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return -1;
  }

  
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  //--- open files and get trees
  std::vector<std::string> input = opts.GetOpt<std::vector<std::string> >("Input.input");
  
  TChain* t = new TChain("chain","the chain");
  for(unsigned int ii = 0; ii < input.size()/2; ++ii)
  {
    std::string inFileName = input.at(0+ii*2);
    std::string treeName = input.at(1+ii*2);
    t -> Add((inFileName+"/"+treeName).c_str());
  }
  long int nEntries = t->GetEntries();
  std::cout << "Added " << nEntries << " entries to the chain" << std::endl;
  
  TreeVars treeVars;
  InitTreeVars(t,treeVars);
  
  
  //------------------------------------
  // create histograms for normalization
  std::cout << ">>> creating histograms for normalization" << std::endl;
  
  TH1F* h1_data_oneCategory = new TH1F("h1_data_oneCategory","",80,100.,180.);
  TH1F* h1_cs_oneCategory   = new TH1F("h1_cs_oneCategory","",80,100.,180.);
  
  TH1F* h1_data_diLepton = new TH1F("h1_data_diLepton","",80,100.,180.);
  TH1F* h1_cs_diLepton   = new TH1F("h1_cs_diLepton","",80,100.,180.);
  
  TH1F* h1_data_singleLepton = new TH1F("h1_data_singleLepton","",80,100.,180.);
  TH1F* h1_cs_singleLepton   = new TH1F("h1_cs_singleLepton","",80,100.,180.);
  
  TH1F* h1_data_oneCategoryNoBTag = new TH1F("h1_data_oneCategoryNoBTag","",80,100.,180.);
  TH1F* h1_cs_oneCategoryNoBTag   = new TH1F("h1_cs_oneCategoryNoBTag","",80,100.,180.);
  
  TH1F* h1_data_diLeptonNoBTag = new TH1F("h1_data_diLeptonNoBTag","",80,100.,180.);
  TH1F* h1_cs_diLeptonNoBTag   = new TH1F("h1_cs_diLeptonNoBTag","",80,100.,180.);
  
  TH1F* h1_data_singleLeptonNoBTag = new TH1F("h1_data_singleLeptonNoBTag","",80,100.,180.);
  TH1F* h1_cs_singleLeptonNoBTag   = new TH1F("h1_cs_singleLeptonNoBTag","",80,100.,180.);
  
  for(long int ii = 0; ii < nEntries; ++ii)
  {
    if( ii%1000 == 0 ) std::cout << ">>> Reading entry " << ii << " / " << nEntries << "\r" << std::flush;
    t -> GetEntry(ii);
    
    
    // common cuts                                                                                                                                                                                                                          
    if( treeVars.dipho_mass < 100 || treeVars.dipho_mass > 180 ) continue;
    if( treeVars.dipho_mass > 115 && treeVars.dipho_mass < 135 ) continue;
    
    
    //-------------
    // one category
    
    if( OneCategorySelection(treeVars,-1,true) == true )
    {
      h1_data_oneCategory -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    else if( OneCategorySelection(treeVars,-2,true) == true )
    {
      h1_cs_oneCategory -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    // one category
    
    
    //--------------------------
    // two categories - dilepton
    
    DiLeptonCategories cat;
    if( DiLeptonSelection(treeVars,-1,true,cat) == true )
    {
      h1_data_diLepton -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    else if( DiLeptonSelection(treeVars,-2,true,cat) == true )
    {
      h1_cs_diLepton -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    // two categories - dilepton
    
    
    //-------------------------------
    // two categories - single lepton
    
    if( SingleLeptonSelection(treeVars,-1,true) == true )
    {
      h1_data_singleLepton -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    else if( SingleLeptonSelection(treeVars,-2,true) == true)
    {
      h1_cs_singleLepton -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    // two categories - single lepton
    
    
    //-------------
    // one category
    
    if( OneCategorySelection(treeVars,-1,true) == true )
    {
      h1_data_oneCategory -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    else if( OneCategorySelection(treeVars,-2,true) == true )
    {
      h1_cs_oneCategory -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    // one category
    
    
    
    //------------------------------------
    // two categories - dilepton - no bTag
    
    if( DiLeptonSelection(treeVars,-1,false,cat) == true )
    {
      h1_data_diLeptonNoBTag -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    else if( DiLeptonSelection(treeVars,-2,false,cat) == true )
    {
      h1_cs_diLeptonNoBTag -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    // two categories - dilepton - no bTag
    
    
    //-----------------------------------------
    // two categories - single lepton - no bTag
    
    if( SingleLeptonSelection(treeVars,-1,false) == true )
    {
      h1_data_singleLeptonNoBTag -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    else if( SingleLeptonSelection(treeVars,-2,false) == true)
    {
      h1_cs_singleLeptonNoBTag -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    // two categories - single lepton - no bTag
    
    
    //-----------------------
    // one category - no bTag
    
    if( OneCategorySelection(treeVars,-1,false) == true )
    {
      h1_data_oneCategoryNoBTag -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    else if( OneCategorySelection(treeVars,-2,false) == true )
    {
      h1_cs_oneCategoryNoBTag -> Fill(treeVars.dipho_mass, treeVars.weight);
    }
    // one category - no bTag
  }
  
  float scaleFactor_oneCategory  =  h1_data_oneCategory->Integral() /  h1_cs_oneCategory->Integral();
  float scaleFactor_diLepton     =     h1_data_diLepton->Integral() /     h1_cs_diLepton->Integral();
  float scaleFactor_singleLepton = h1_data_singleLepton->Integral() / h1_cs_singleLepton->Integral();
  
  float scaleFactor_oneCategoryNoBTag  =  h1_data_oneCategoryNoBTag->Integral() /  h1_cs_oneCategoryNoBTag->Integral();
  float scaleFactor_diLeptonNoBTag     =     h1_data_diLeptonNoBTag->Integral() /     h1_cs_diLeptonNoBTag->Integral();
  float scaleFactor_singleLeptonNoBTag = h1_data_singleLeptonNoBTag->Integral() / h1_cs_singleLeptonNoBTag->Integral();
  
  
  //-------------------------------------------------------
  //--- create a new file + a clone of old tree in new file
  std::cout << ">>> creating control sample" << std::endl;
  
  std::string outputFileName = opts.GetOpt<std::string>("Output.outputFileName");
  std::string outputTreeName = opts.GetOpt<std::string>("Output.outputTreeName");
  TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");

  TTree* outputTree_oneCategory = t->CloneTree(0);
  outputTree_oneCategory -> SetName((outputTreeName+"_oneCategory").c_str());
  TTree* outputTree_diLepton = t->CloneTree(0);
  outputTree_diLepton -> SetName((outputTreeName+"_diLepton").c_str());
  TTree* outputTree_singleLepton = t->CloneTree(0);
  outputTree_singleLepton -> SetName((outputTreeName+"_singleLepton").c_str());
  
  TTree* outputTree_oneCategoryNoBTag = t->CloneTree(0);
  outputTree_oneCategoryNoBTag -> SetName((outputTreeName+"_oneCategoryNoBTag").c_str());
  TTree* outputTree_diLeptonNoBTag = t->CloneTree(0);
  outputTree_diLeptonNoBTag -> SetName((outputTreeName+"_diLeptonNoBTag").c_str());
  TTree* outputTree_singleLeptonNoBTag = t->CloneTree(0);
  outputTree_singleLeptonNoBTag -> SetName((outputTreeName+"_singleLeptonNoBTag").c_str());
  
  for(long int ii = 0; ii < nEntries; ++ii)
  {
    if( ii%1000 == 0 ) std::cout << ">>> Reading entry " << ii << " / " << nEntries << "\r" << std::flush;
    t -> GetEntry(ii);
    
    float oldWeight = treeVars.weight;
    
    if( OneCategorySelection(treeVars,-1,true) == false && OneCategorySelection(treeVars,-2,true) == true )
    {
      treeVars.weight = oldWeight * scaleFactor_oneCategory;
      outputTree_oneCategory -> Fill();
    }
    
    DiLeptonCategories cat;
    if( DiLeptonSelection(treeVars,-1,true,cat) == false && DiLeptonSelection(treeVars,-2,true,cat) == true )
    {
      treeVars.weight = oldWeight * scaleFactor_diLepton;
      outputTree_diLepton -> Fill();
    }
    
    if( SingleLeptonSelection(treeVars,-1,true) == false && SingleLeptonSelection(treeVars,-2,true) == true )
    {
      treeVars.weight = oldWeight * scaleFactor_singleLepton;
      outputTree_singleLepton -> Fill();
    }
    
    if( OneCategorySelection(treeVars,-1,false) == false && OneCategorySelection(treeVars,-2,false) == true )
    {
      treeVars.weight = oldWeight * scaleFactor_oneCategoryNoBTag;
      outputTree_oneCategoryNoBTag -> Fill();
    }
    
    if( DiLeptonSelection(treeVars,-1,false,cat) == false && DiLeptonSelection(treeVars,-2,false,cat) == true )
    {
      treeVars.weight = oldWeight * scaleFactor_diLeptonNoBTag;
      outputTree_diLeptonNoBTag -> Fill();
    }
    
    if( SingleLeptonSelection(treeVars,-1,false) == false && SingleLeptonSelection(treeVars,-2,false) == true )
    {
      treeVars.weight = oldWeight * scaleFactor_singleLeptonNoBTag;
      outputTree_singleLeptonNoBTag -> Fill();
    }
  }
  
  outputTree_oneCategory -> AutoSave();
  outputTree_diLepton -> AutoSave();
  outputTree_singleLepton -> AutoSave();
  
  outputTree_oneCategoryNoBTag -> AutoSave();
  outputTree_diLeptonNoBTag -> AutoSave();
  outputTree_singleLeptonNoBTag -> AutoSave();
  
  h1_data_oneCategory  -> Write();
  h1_cs_oneCategory    -> Write();
  h1_data_diLepton     -> Write();
  h1_cs_diLepton       -> Write();
  h1_data_singleLepton -> Write();
  h1_cs_singleLepton   -> Write();
  
  h1_data_oneCategoryNoBTag  -> Write();
  h1_cs_oneCategoryNoBTag    -> Write();
  h1_data_diLeptonNoBTag     -> Write();
  h1_cs_diLeptonNoBTag       -> Write();
  h1_data_singleLeptonNoBTag -> Write();
  h1_cs_singleLeptonNoBTag   -> Write();
     
  outputFile -> Close();
  
  
  return 0;
}
