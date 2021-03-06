#include "interface/AnalysisUtils.h"



float DeltaEta(const float& eta1, const float& eta2)
{
  return fabs( eta1 - eta2 );
}

float DeltaPhi(const float& phi1, const float& phi2)
{
  float dphi = fabs( phi1 - phi2 );
  if( dphi > PI ) dphi = 2*PI - dphi;
  return dphi;
}

float DeltaR(const float& eta1, const float& phi1,
             const float& eta2, const float& phi2)
{
  return sqrt( DeltaEta(eta1,eta2)*DeltaEta(eta1,eta2) + 
               DeltaPhi(phi1,phi2)*DeltaPhi(phi1,phi2) );
}




bool DiMuSelections(TLorentzVector mu1, TLorentzVector mu2, float charge1, float charge2, TLorentzVector ph1, TLorentzVector ph2, float Iso1, float Iso2)
{
  // if(charge1*charge2!=-1) return 0;
  
  if( std::max(mu1.Pt(),mu2.Pt()) < 25. ) return 0;
  if( mu1.Pt() < 10 ) return 0;
  if( mu2.Pt() < 10 ) return 0;
  
  if( fabs(mu1.Eta()) > 2.4 ) return 0;
  if( fabs(mu2.Eta()) > 2.4 ) return 0;
  
  if( fabs((mu1+mu2).M() - MZ) < 5. ) return 0;
  
  if( DeltaR( ph1.Eta(), ph1.Phi(), mu1.Eta(), mu1.Phi() ) < 0.2 ) return 0;
  if( DeltaR( ph1.Eta(), ph1.Phi(), mu2.Eta(), mu2.Phi() ) < 0.2 ) return 0;
  if( DeltaR( ph2.Eta(), ph2.Phi(), mu1.Eta(), mu1.Phi() ) < 0.2 ) return 0;
  if( DeltaR( ph2.Eta(), ph2.Phi(), mu2.Eta(), mu2.Phi() ) < 0.2 ) return 0;
  
  if( DeltaR( mu1.Eta(), mu1.Phi(), mu2.Eta(), mu2.Phi() ) < 0.1 ) return 0;
  
  // if(Iso1>0.25 || Iso2>0.25) return 0;
  // if(Iso1>0.06 || Iso2>0.06) return 0;
  
  return 1;
}


bool DiEleSelections(TLorentzVector ele1, TLorentzVector ele2, float charge1, float charge2, TLorentzVector ph1, TLorentzVector ph2, float iso1, float iso2, float dTrk1, float dTrk2)
{
  // if(charge1*charge2!=-1) return 0;
  
  if( std::max(ele1.Pt(),ele2.Pt()) < 20. ) return 0;
  if( ele1.Pt() < 15. ) return 0;
  if( ele2.Pt() < 15. ) return 0;
  
  if( fabs(ele1.Eta())>2.5 || ( fabs(ele1.Eta())>1.4442 && fabs(ele1.Eta())<1.566) ) return 0;
  if( fabs(ele2.Eta())>2.5 || ( fabs(ele2.Eta())>1.4442 && fabs(ele2.Eta())<1.566) ) return 0;
  
  if( fabs((ele1+ele2).M() - MZ) < 5. ) return 0;
  if( fabs((ph1+ele1).M()  - MZ) < 5. ) return 0;
  if( fabs((ph1+ele2).M()  - MZ) < 5. ) return 0;
  if( fabs((ph2+ele1).M()  - MZ) < 5. ) return 0;
  if( fabs((ph2+ele2).M()  - MZ) < 5. ) return 0;
  
  if( DeltaR( ph1.Eta(), ph1.Phi(), ele1.Eta(), ele1.Phi() ) < 0.2 ) return 0;
  if( DeltaR( ph1.Eta(), ph1.Phi(), ele2.Eta(), ele2.Phi() ) < 0.2 ) return 0;
  if( DeltaR( ph2.Eta(), ph2.Phi(), ele1.Eta(), ele1.Phi() ) < 0.2 ) return 0;
  if( DeltaR( ph2.Eta(), ph2.Phi(), ele2.Eta(), ele2.Phi() ) < 0.2 ) return 0;
  
  if( DeltaR( ele1.Eta(), ele1.Phi(), ele2.Eta(), ele2.Phi() ) < 0.2 ) return 0;
  
  // if((abs(ele1.Eta())<= 1.479 && iso1 > 0.045) || (abs(ele1.Eta())> 1.479 && iso1 > 0.08)) return 0;
  // if((abs(ele2.Eta())<= 1.479 && iso2 > 0.045) || (abs(ele2.Eta())> 1.479 && iso2 > 0.08)) return 0;
  
  // if(dTrk1>0.35) return 0;
  // if(dTrk2>0.35) return 0;
  
  return 1;
}


bool MixedSelections(TLorentzVector mu, TLorentzVector ele, float charge1, float charge2, TLorentzVector ph1, TLorentzVector ph2)
{
  if( charge1*charge2 != -1 ) return 0;
  
  if( ele.Pt() < 10. ) return 0;
  if( mu.Pt()  < 10. ) return 0;
  
  if( fabs(ele.Eta()) > 2.5 || (abs(ele.Eta()) > 1.4442 && fabs(ele.Eta()) < 1.566) ) return 0;
  if( fabs(mu.Eta()) > 2.4 ) return 0;
  
  if( DeltaR( ph1.Eta(), ph1.Phi(),  mu.Eta(),  mu.Phi() ) < 0.2 ) return 0;
  if( DeltaR( ph2.Eta(), ph2.Phi(),  mu.Eta(),  mu.Phi() ) < 0.2 ) return 0;
  if( DeltaR( ph1.Eta(), ph1.Phi(), ele.Eta(), ele.Phi() ) < 0.2 ) return 0;
  if( DeltaR( ph2.Eta(), ph2.Phi(), ele.Eta(), ele.Phi() ) < 0.2 ) return 0;
  
  if( DeltaR( mu.Eta(), mu.Phi(), ele.Eta(), ele.Phi() ) < 0.1 ) return 0;
  
  return 1;
}


bool SingleMuSelections(TLorentzVector mu1,  TLorentzVector ph1, TLorentzVector ph2, float iso)
{
  if( mu1.Pt() < 5. ) return 0;
  
  if( fabs(mu1.Eta()) > 2.4 ) return 0;
  
  if( DeltaR( ph1.Eta(), ph1.Phi(), mu1.Eta(), mu1.Phi() ) < 0.2 ) return 0;
  if( DeltaR( ph2.Eta(), ph2.Phi(), mu1.Eta(), mu1.Phi() ) < 0.2 ) return 0;
  
  if( iso > 0.20 ) return 0;
  
  return 1;
}


bool SingleEleSelections(TLorentzVector ele1, TLorentzVector ph1, TLorentzVector ph2, float iso, float drTrk)
{
  if( ele1.Pt() < 10. ) return 0;
  
  if( fabs(ele1.Eta()) > 2.5 || ( fabs(ele1.Eta()) > 1.4442 && fabs(ele1.Eta()) < 1.566 ) )  return 0;
  
  if( DeltaR( ph1.Eta(), ph1.Phi(), ele1.Eta(), ele1.Phi() ) < 0.2 ) return 0;
  if( DeltaR( ph2.Eta(), ph2.Phi(), ele1.Eta(), ele1.Phi() ) < 0.2 ) return 0;
  
  if( fabs((ph1+ele1).M() - MZ) < 5. ) return 0;
  if( fabs((ph2+ele1).M() - MZ) < 5. ) return 0;
  
  // if( drTrk > 0.35 ) return 0;
  
  return 1;
}



void MakePlot(TH1F** histos, TString title)
{
  histos[0] -> SetLineWidth(3);			//ttH
  histos[0] -> SetLineColor(kRed + 1);
  histos[0] -> SetFillStyle(0);
  
  histos[1] -> SetLineWidth(3);			//ggH
  histos[1] -> SetLineColor(kGreen + 2);
  histos[1] -> SetFillStyle(0);
  
  histos[2] -> SetLineWidth(3);			//VBF
  histos[2] -> SetLineColor(kAzure);
  histos[2] -> SetFillStyle(0);
  
  histos[3] -> SetLineWidth(3);			//VH
  histos[3] -> SetLineColor(kViolet - 2);
  histos[3] -> SetFillStyle(0);
  
  histos[4] -> SetLineWidth(3);			//bbH
  histos[4] -> SetLineColor(kOrange);
  histos[4] -> SetFillStyle(0);
  
  histos[5] -> SetLineWidth(3);			//tHq
  histos[5] -> SetLineColor(kAzure + 8);
  histos[5] -> SetFillStyle(0);
  
  histos[6] -> SetLineWidth(3);			//tHW
  histos[6] -> SetLineColor(kViolet + 2);
  histos[6] -> SetFillStyle(0);
  
  
  
  histos[7] -> SetMarkerStyle(20);		//Data
  histos[7] -> SetMarkerSize(1);
  histos[7] -> SetMarkerColor(kBlack);
  histos[7] -> SetFillStyle(0);
  
  
  
  histos[8] -> SetLineWidth(1);			//Diphoton
  histos[8] -> SetFillColor(kAzure + 1);
  histos[8] -> SetFillStyle(1001);
  
  histos[9] -> SetLineWidth(1);			//Gamma + jets
  histos[9] -> SetFillStyle(1001);
  histos[9] -> SetFillColor(kYellow - 4);
  
  histos[10] -> SetLineWidth(1);			//QCD
  histos[10] -> SetFillColor(kTeal + 9);
  histos[10] -> SetFillStyle(1001);
  
  histos[11] -> SetLineWidth(1);			//ttGG
  histos[11] -> SetFillColor(kMagenta + 1);
  histos[11] -> SetFillStyle(1001);
  
  histos[12] -> SetLineWidth(1);			//ttGJets
  histos[12] -> SetFillColor(kMagenta + 1);
  histos[12] -> SetFillStyle(1001);
  
  histos[13] -> SetLineWidth(1);			//ttJets
  histos[13] -> SetFillColor(kMagenta + 1);
  histos[13] -> SetFillStyle(1001);
  
  
  TLegend* leg = new TLegend(0.65, 0.70, 0.9, 0.85);
  leg -> AddEntry(histos[0], "ttH", "l");
  leg -> AddEntry(histos[1], "ggH", "l");
  leg -> AddEntry(histos[2], "VBF", "l");
  leg -> AddEntry(histos[3], "VH", "l");
  leg -> AddEntry(histos[4], "bbH", "l");
  leg -> AddEntry(histos[5], "tHq", "l");
  leg -> AddEntry(histos[6], "tHW", "l");
  leg -> AddEntry(histos[7], "Data sidebands", "p");
  
  TLegend* leg2 = new TLegend(0.65, 0.70, 0.9, 0.85);
  leg2 -> AddEntry(histos[7], "Data sidebands", "p");
  leg2 -> AddEntry(histos[8], "Diphotons", "f");
  leg2 -> AddEntry(histos[9], "Gamma + jets", "f");
  leg2 -> AddEntry(histos[10], "QCD", "f");
  leg2 -> AddEntry(histos[11], "ttGJets", "f");
  
  
  histos[11] -> Add(histos[12]);
  histos[11] -> Add(histos[13]);
  
  
  TCanvas* c = new TCanvas();
  c -> cd();
  
  float m = std::max(std::max(histos[0]->GetMaximum()/histos[0]->Integral(), histos[1]->GetMaximum()/histos[1]->Integral()), std::max(histos[2]->GetMaximum()/histos[2]->Integral(), histos[3]->GetMaximum()/histos[3]->Integral()));
  float m2 = std::max(std::max(histos[4]->GetMaximum()/histos[4]->Integral(), histos[5]->GetMaximum()/histos[5]->Integral()), std::max(histos[6]->GetMaximum()/histos[6]->Integral(), histos[7]->GetMaximum()/histos[7]->Integral()));
  m = std::max((double)m, (double)m2);
  
  TH1F* axis = new TH1F(*histos[0]);
  axis -> SetMarkerSize(0);
  axis -> SetLineWidth(0);
  axis -> GetYaxis() -> SetTitleOffset(1.5);
  axis -> GetYaxis() -> SetRangeUser(0, 1.1*m);
  
  axis -> Draw("histo");
  histos[1] -> DrawNormalized("histo SAME");
  histos[2] -> DrawNormalized("histo SAME");
  histos[3] -> DrawNormalized("histo SAME");
  histos[4] -> DrawNormalized("histo SAME");
  histos[5] -> DrawNormalized("histo SAME");
  histos[6] -> DrawNormalized("histo SAME");
  histos[0] -> DrawNormalized("histo SAME");
  histos[7] -> DrawNormalized("SAME E1");
  
  leg -> Draw("SAME");
  CMS_lumi(c, 0, 0);
  
  c -> SaveAs("c_" + title + "Signal.png");
  c -> SaveAs("c_" +title + "Signal.pdf");
  
  
  TCanvas* c2 = new TCanvas();
  c2 -> cd();
  
  histos[9] -> Add(histos[8]);
  histos[10] -> Add(histos[9]);
  histos[11] -> Add(histos[10]);
  
  
  double max2 = histos[7]->GetBinContent(histos[7]->GetMaximumBin());
  if(histos[11]->GetBinContent(histos[11]->GetMaximumBin()) > max2 )
    max2 = histos[11]->GetBinContent(histos[11]->GetMaximumBin());
  histos[11]-> GetYaxis() -> SetRangeUser(0, 1.1*max2);
  histos[11]-> GetYaxis() -> SetTitleSize(0.05);
  histos[11]-> GetYaxis() -> SetTitleFont(42);
  histos[11]-> GetYaxis() -> SetLabelSize(0.045);
  histos[11]-> GetYaxis() -> SetLabelFont(42);
  histos[11]-> GetXaxis() -> SetLabelSize(0);
  
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 0.97);
  pad1->SetBottomMargin(0.035);
  pad1->Draw();
  pad1->cd();
  
  histos[11] -> Draw("histo");
  histos[10] -> Draw("histo SAME");
  histos[9] -> Draw("histo SAME");
  histos[8] -> Draw("histo SAME");
  histos[7] -> Draw("SAME E1");
  leg2 -> Draw("SAME");
  
  c2->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.10, 1, 0.35);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0);
  pad2->Draw();
  pad2->cd();
  
  TH1F *h = (TH1F*)histos[7]->Clone("h");
  h->SetLineColor(kBlack);
  h->SetMinimum(0.5);  // Define Y ..
  h->SetMaximum(1.5); // .. range
  h->Sumw2();
  h->SetStats(0);      // No statistics on lower plot
  h->Divide(histos[11]);
  h->SetMarkerStyle(21);
  h -> SetTitle("");
  h-> GetYaxis() -> SetTitle("Data/MC");
  
  
  // Y axis ratio plot settings
  h->GetYaxis()->SetNdivisions(-10);
  h->GetYaxis()->SetTitleSize(0.13);
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(0.5);
  h->GetYaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)
  h->GetYaxis()->SetLabelSize(0.12);
  
  // X axis ratio plot settings
  h->GetXaxis()->SetTitleSize(0.15);
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)
  h->GetXaxis()->SetLabelSize(0.12);
  
  h->Draw("EP");
  
  
  CMS_lumi(c2, 0, 10);
  
  c2 -> SaveAs("c_" +title + "Bkg.png");
  c2 -> SaveAs("c_" +title + "Bkg.pdf");
  
  return;
}



void MakePlot2(std::map<std::string,TH1F*>& histos, TString title)
{
  histos["ttH"] -> SetLineWidth(3);			//ttH
  histos["ttH"] -> SetLineColor(kRed + 1);
  histos["ttH"] -> SetFillStyle(0);
  
  histos["ggH"] -> SetLineWidth(3);			//ggH
  histos["ggH"] -> SetLineColor(kGreen + 2);
  histos["ggH"] -> SetFillStyle(0);
  
  histos["VBF"] -> SetLineWidth(3);			//VBF
  histos["VBF"] -> SetLineColor(kAzure);
  histos["VBF"] -> SetFillStyle(0);
  
  histos["VH"] -> SetLineWidth(3);			//VH
  histos["VH"] -> SetLineColor(kViolet - 2);
  histos["VH"] -> SetFillStyle(0);
  
  // histos["bbH"] -> SetLineWidth(3);			//bbH
  // histos["bbH"] -> SetLineColor(kOrange);
  // histos["bbH"] -> SetFillStyle(0);
  
  // histos["tHq"] -> SetLineWidth(3);			//tHq
  // histos["tHq"] -> SetLineColor(kAzure + 8);
  // histos["tHq"] -> SetFillStyle(0);
  
  // histos["tHW"] -> SetLineWidth(3);			//tHW
  // histos["tHW"] -> SetLineColor(kViolet + 2);
  // histos["tHW"] -> SetFillStyle(0);
  
  histos["data"] -> SetMarkerStyle(20);		//Data
  histos["data"] -> SetMarkerSize(1);
  histos["data"] -> SetMarkerColor(kBlack);
  histos["data"] -> SetFillStyle(0);
  
  /*
  histos[8] -> SetLineWidth(1);			//Diphoton
  histos[8] -> SetFillColor(kAzure + 1);
  histos[8] -> SetFillStyle(1001);
  
  histos[9] -> SetLineWidth(1);			//Gamma + jets
  histos[9] -> SetFillStyle(1001);
  histos[9] -> SetFillColor(kYellow - 4);
  
  histos[10] -> SetLineWidth(1);			//QCD
  histos[10] -> SetFillColor(kTeal + 9);
  histos[10] -> SetFillStyle(1001);
  
  histos[11] -> SetLineWidth(1);			//ttGG
  histos[11] -> SetFillColor(kMagenta + 1);
  histos[11] -> SetFillStyle(1001);
  
  histos[12] -> SetLineWidth(1);			//ttGJets
  histos[12] -> SetFillColor(kMagenta + 1);
  histos[12] -> SetFillStyle(1001);
  
  histos[13] -> SetLineWidth(1);			//ttJets
  histos[13] -> SetFillColor(kMagenta + 1);
  histos[13] -> SetFillStyle(1001);
  */
  
  histos[Form("CS_%s",title.Data())] -> SetLineWidth(3);			//CS
  histos[Form("CS_%s",title.Data())] -> SetLineColor(kAzure);
  histos[Form("CS_%s",title.Data())] -> SetFillStyle(1);
  histos[Form("CS_%s",title.Data())] -> SetFillColor(kAzure-9);
  
  
  TLegend* leg2 = new TLegend(0.65, 0.70, 0.9, 0.85);
  leg2 -> AddEntry(histos["ttH"], "Signal", "l");
  leg2 -> AddEntry(histos["data"], "Data sidebands", "p");
  leg2 -> AddEntry(histos[Form("CS_%s",title.Data())], "Control Sample", "f");
  //	leg2 -> AddEntry(histos[9], "Gamma + jets", "f");
  //	leg2 -> AddEntry(histos[10], "QCD", "f");
  //	leg2 -> AddEntry(histos[11], "ttGJets", "f");
  
  
  // histos[11] -> Add(histos[12]);
  // histos[11] -> Add(histos[13]);
  
  TH1F* signal = new TH1F(*histos["ttH"]);
  signal -> Add(histos["ggH"]);
  signal -> Add(histos["VBF"]);
  signal -> Add(histos["VH"]);
  // signal -> Add(histos["bbH"]);
  // signal -> Add(histos["tHq"]);
  // signal -> Add(histos["tHW"]);
  
  std::cout << "Tag Purity: ttH " << std::setprecision(2) << histos["ttH"]->Integral()/signal->Integral()*100 << "%, (" << std::setprecision(2) << histos["ttH"]->Integral() << ") events" << std::endl;
  std::cout << "Tag Purity: ggH " << std::setprecision(2) << histos["ggH"]->Integral()/signal->Integral()*100 << "%, (" << std::setprecision(4) << histos["ggH"]->Integral() << ") events" << std::endl;
  std::cout << "Tag Purity: VBF " << std::setprecision(2) << histos["VBF"]->Integral()/signal->Integral()*100 << "%, (" << std::setprecision(4) << histos["VBF"]->Integral() << ") events" << std::endl;
  std::cout << "Tag Purity: VH  " << std::setprecision(2) << histos["VH"] ->Integral()/signal->Integral()*100 << "%, (" << std::setprecision(4) << histos["VH"] ->Integral() << ") events" << std::endl;
  // std::cout << "Tag Purity: bbH " << std::setprecision(2) << histos["bbH"]->Integral()/signal->Integral()*100 << "%, (" << std::setprecision(4) << histos["bbH"]->Integral() << ") events" << std::endl;
  // std::cout << "Tag Purity: tHq " << std::setprecision(2) << histos["tHq"]->Integral()/signal->Integral()*100 << "%, (" << std::setprecision(4) << histos["tHq"]->Integral() << ") events" << std::endl;
  // std::cout << "Tag Purity: tHW " << std::setprecision(2) << histos["tHW"]->Integral()/signal->Integral()*100 << "%, (" << std::setprecision(4) << histos["tHW"]->Integral() << ") events" << std::endl;
  
  // histos[9] -> Add(histos[8]);
  // histos[10] -> Add(histos[9]);
  // histos[11] -> Add(histos[10]);
  
  
  TCanvas* c = new TCanvas();
  c -> cd();
  
  TH1F* axis = new TH1F(*histos["ttH"]);
  axis -> SetMarkerSize(0);
  axis -> SetLineWidth(0);
  axis -> GetYaxis() -> SetTitleOffset(1.5);
  axis -> GetYaxis() -> SetRangeUser(0, 1.1*histos["data"]->GetMaximum());
  
  axis -> Draw("histo");
  
  //	histos[11] -> Draw("histo SAME");
  //	histos[10] -> Draw("histo SAME");
  //	histos[9] -> Draw("histo SAME");
  //	histos[8] -> Draw("histo SAME");
  histos[Form("CS_%s",title.Data())] -> Draw("histo,same");
  histos["data"] -> Draw("SAME E1");
  signal -> Draw("histo SAME");
  leg2 -> Draw("SAME");
  
  CMS_lumi(c, 0, 0);
  
  c -> SaveAs("c_" +title + "Signal.png");
  c -> SaveAs("c_" +title + "Signal.pdf");
  
  return;
}



bool OneCategorySelection(const TreeVars& treeVars, const int& type, const bool& bTagSelection, const ControlSampleType& csType)
{
  if( csType == kInvertLeptons && type == -2 )
  {
    int njet = 0;
    int nbjet_loose = 0;
    int nbjet_medium = 0;
    int nbjet_tight = 0;
    float jetPtMax = -999.;
    
    for(int jIndex = 0; jIndex < 6; ++jIndex)
    {
      if( treeVars.jet_pt[jIndex] > 25. && fabs(treeVars.jet_eta[jIndex]) < 2.4 )
      {
        ++njet;
        if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdLooseDeep  ) ++nbjet_loose;
        if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdMediumDeep ) ++nbjet_medium;
        if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdTightDeep  ) ++nbjet_tight;
        
        if( treeVars.jet_pt[jIndex] > jetPtMax ) jetPtMax = treeVars.jet_pt[jIndex];
      }
    }
    
    if( njet >= 3 && jetPtMax > 50. && ( (bTagSelection && nbjet_medium >= 1) || (!bTagSelection && nbjet_medium == 0) ) )
      return true;
    else
      return false;
  }
  
  std::vector<float> goodLeptons;
  
  if(treeVars.mu_pt[0]>10.)
  {
    TLorentzVector mu;
    mu.SetPtEtaPhiE(treeVars.mu_pt[0], treeVars.mu_eta[0], treeVars.mu_phi[0], treeVars.mu_energy[0]);
    
    TLorentzVector ph1, ph2;
    ph1.SetPtEtaPhiE(treeVars.dipho_leadPt, treeVars.dipho_leadEta, treeVars.dipho_leadPhi, treeVars.dipho_leadEnergy);
    ph2.SetPtEtaPhiE(treeVars.dipho_subleadPt, treeVars.dipho_subleadEta, treeVars.dipho_subleadPhi, treeVars.dipho_subleadEnergy);
    
    if(treeVars.mu_IDVector[0][oneCatMuID] && SingleMuSelections(mu, ph1, ph2, treeVars.mu_miniIso[0]))
      goodLeptons.push_back(0);
  }
  
  if(treeVars.mu_pt[1]>10.)
  {
    TLorentzVector mu;
    mu.SetPtEtaPhiE(treeVars.mu_pt[1], treeVars.mu_eta[1], treeVars.mu_phi[1], treeVars.mu_energy[1]);
    
    TLorentzVector ph1, ph2;
    ph1.SetPtEtaPhiE(treeVars.dipho_leadPt, treeVars.dipho_leadEta, treeVars.dipho_leadPhi, treeVars.dipho_leadEnergy);
    ph2.SetPtEtaPhiE(treeVars.dipho_subleadPt, treeVars.dipho_subleadEta, treeVars.dipho_subleadPhi, treeVars.dipho_subleadEnergy);
    
    if(treeVars.mu_IDVector[1][oneCatMuID] && SingleMuSelections(mu, ph1, ph2, treeVars.mu_miniIso[1]))
      goodLeptons.push_back(1);
  }
  
  if(treeVars.ele_pt[0]>10.)
  {
    TLorentzVector ele;
    ele.SetPtEtaPhiE(treeVars.ele_pt[0], treeVars.ele_eta[0], treeVars.ele_phi[0], treeVars.ele_energy[0]);
    
    TLorentzVector ph1, ph2;
    ph1.SetPtEtaPhiE(treeVars.dipho_leadPt, treeVars.dipho_leadEta, treeVars.dipho_leadPhi, treeVars.dipho_leadEnergy);
    ph2.SetPtEtaPhiE(treeVars.dipho_subleadPt, treeVars.dipho_subleadEta, treeVars.dipho_subleadPhi, treeVars.dipho_subleadEnergy);
    
    float dTrk = sqrt(treeVars.ele_dEtaTrk[0]*treeVars.ele_dEtaTrk[0] + treeVars.ele_dPhiTrk[0]*treeVars.ele_dPhiTrk[0]);
    
    bool passSingleEleSelection = SingleEleSelections(ele, ph1, ph2, treeVars.ele_miniIso[0], dTrk);
    if(treeVars.ele_IDVector[0][oneCatEleID] && passSingleEleSelection)
      goodLeptons.push_back(2);
  }
  
  if(treeVars.ele_pt[1]>10.)
  {
    TLorentzVector ele;
    ele.SetPtEtaPhiE(treeVars.ele_pt[1], treeVars.ele_eta[1], treeVars.ele_phi[1], treeVars.ele_energy[1]);
    
    TLorentzVector ph1, ph2;
    ph1.SetPtEtaPhiE(treeVars.dipho_leadPt, treeVars.dipho_leadEta, treeVars.dipho_leadPhi, treeVars.dipho_leadEnergy);
    ph2.SetPtEtaPhiE(treeVars.dipho_subleadPt, treeVars.dipho_subleadEta, treeVars.dipho_subleadPhi, treeVars.dipho_subleadEnergy);
    
    float dTrk = sqrt(treeVars.ele_dEtaTrk[1]*treeVars.ele_dEtaTrk[1] + treeVars.ele_dPhiTrk[1]*treeVars.ele_dPhiTrk[1]);
    
    bool passSingleEleSelection = SingleEleSelections(ele, ph1, ph2, treeVars.ele_miniIso[1], dTrk);
    if(treeVars.ele_IDVector[1][oneCatEleID] && passSingleEleSelection)
      goodLeptons.push_back(3);
  }
  
  
  bool accept = false;
  
  for(unsigned int j = 0; j < goodLeptons.size(); j++)
  { 
    int njet = 0;
    int nbjet_loose = 0;
    int nbjet_medium = 0;
    int nbjet_tight = 0;
    float jetPtMax = -999.;
    
    for(int jIndex = 0; jIndex < 6; ++jIndex)
    {
      if( treeVars.jet_pt[jIndex] > 25. && fabs(treeVars.jet_eta[jIndex]) < 2.4 )
      {
        float DR = -1.;
        if( goodLeptons.at(j)==0 ) DR = DeltaR(treeVars.jet_eta[jIndex], treeVars.jet_phi[jIndex], treeVars.mu_eta[0], treeVars.mu_phi[0]);
        if( goodLeptons.at(j)==1 ) DR = DeltaR(treeVars.jet_eta[jIndex], treeVars.jet_phi[jIndex], treeVars.mu_eta[1], treeVars.mu_phi[1]);
        if( goodLeptons.at(j)==2 ) DR = DeltaR(treeVars.jet_eta[jIndex], treeVars.jet_phi[jIndex], treeVars.ele_eta[0], treeVars.ele_phi[0]);
        if( goodLeptons.at(j)==3 ) DR = DeltaR(treeVars.jet_eta[jIndex], treeVars.jet_phi[jIndex], treeVars.ele_eta[1], treeVars.ele_phi[1]);
        if( DR < 0.4 ) continue;
        
        ++njet;
        if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdLooseDeep  ) ++nbjet_loose;
        if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdMediumDeep ) ++nbjet_medium;
        if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdTightDeep  ) ++nbjet_tight;
        
        if( treeVars.jet_pt[jIndex] > jetPtMax ) jetPtMax = treeVars.jet_pt[jIndex];
      }
    }
    
    if( type != -2 )
    {
      if( njet >= 3 && jetPtMax > 50. && ( (bTagSelection && nbjet_medium >= 1) || (!bTagSelection && nbjet_medium == 0) ) )
      {
        accept = true;
        break;
      }
    }
    else if( csType == kInvertBTag )
    {
      if( njet >= 3 && jetPtMax > 50. && ( (bTagSelection && nbjet_medium == 0) || (!bTagSelection && nbjet_medium == 0) ) )
      {
        accept = true;
        break;
      }      
    }
  }
  
  return accept;
}



bool DiLeptonSelection(const TreeVars& treeVars, const int& type, const bool& bTagSelection, DiLeptonCategories& cat, const ControlSampleType& csType, const bool& verbosity)
{
  if( csType == kInvertLeptons && type == -2 )
  {
    int njet = 0;
    int nbjet_loose = 0;
    int nbjet_medium = 0;
    int nbjet_tight = 0;
    float jetPtMax = -999.;
    
    for(int jIndex = 0; jIndex < nJet; ++jIndex)
    {
      if( treeVars.jet_pt[jIndex] > 25. && fabs(treeVars.jet_eta[jIndex]) < 2.4 )
      {
        ++njet;
        if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdLooseDeep  ) ++nbjet_loose;
        if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdMediumDeep ) ++nbjet_medium;
        if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdTightDeep  ) ++nbjet_tight;
        
        if( treeVars.jet_pt[jIndex] > jetPtMax ) jetPtMax = treeVars.jet_pt[jIndex];
      }
    }
    
    if( njet >= 1 && jetPtMax > 50. && ( (bTagSelection && nbjet_medium >= 1) || (!bTagSelection && nbjet_medium == 0) ) )
      return true;
    else
      return false;
  }
  
  
  std::vector<int> goodLeptons;
  
  for(int ii = 0; ii < nLep; ++ii)
    if( treeVars.mu_pt[0] > 10. && treeVars.mu_IDVector[ii][diMuID] == 1 )
      goodLeptons.push_back(ii);
  
  for(int ii = 0; ii < nLep; ++ii)
    if( treeVars.ele_pt[ii] > 10. && treeVars.ele_IDVector[ii][diEleID] )
      goodLeptons.push_back(nLep+ii);
  
  if( verbosity )
    for(unsigned int ii = 0; ii < goodLeptons.size(); ++ii)
      std::cout << "goodLeptons[" << ii << "] = " << goodLeptons.at(ii) << std::endl;
  
  // the di photon
  TLorentzVector ph1, ph2;
  ph1.SetPtEtaPhiE(treeVars.dipho_leadPt,    treeVars.dipho_leadEta,    treeVars.dipho_leadPhi,    treeVars.dipho_leadEnergy);
  ph2.SetPtEtaPhiE(treeVars.dipho_subleadPt, treeVars.dipho_subleadEta, treeVars.dipho_subleadPhi, treeVars.dipho_subleadEnergy);
  
  // select lepton pairs
  std::vector<std::pair<int,int> > selectedPairs;
  for(unsigned int l1Count = 0; l1Count < goodLeptons.size(); ++l1Count)
  {
    for(unsigned int l2Count = l1Count+1; l2Count < goodLeptons.size(); ++l2Count)
    {
      int id1 = goodLeptons[l1Count]%nLep;
      int id2 = goodLeptons[l2Count]%nLep;
      
      if( goodLeptons[l1Count] < nLep && goodLeptons[l2Count] < nLep ) // di-mu
      {
        TLorentzVector mu1, mu2;
        mu1.SetPtEtaPhiE(treeVars.mu_pt[id1], treeVars.mu_eta[id1], treeVars.mu_phi[id1], treeVars.mu_energy[id1]);
        mu2.SetPtEtaPhiE(treeVars.mu_pt[id2], treeVars.mu_eta[id2], treeVars.mu_phi[id2], treeVars.mu_energy[id2]);
        
        bool passDiMuonSelection = DiMuSelections(mu1, mu2, treeVars.mu_charge[id1], treeVars.mu_charge[id2], ph1, ph2, treeVars.mu_miniIso[id1], treeVars.mu_miniIso[id2]);
        if( !passDiMuonSelection ) continue;
        
        std::pair<int,int> leptons(goodLeptons[l1Count],goodLeptons[l2Count]);
        selectedPairs.push_back(leptons);
      }
      
      else if( goodLeptons[l1Count] >= nLep && goodLeptons[l2Count] >= nLep ) // di-ele
      {
        TLorentzVector ele1, ele2;
        ele1.SetPtEtaPhiE(treeVars.ele_pt[id1], treeVars.ele_eta[id1], treeVars.ele_phi[id1], treeVars.ele_energy[id1]);
        ele2.SetPtEtaPhiE(treeVars.ele_pt[id2], treeVars.ele_eta[id2], treeVars.ele_phi[id2], treeVars.ele_energy[id2]);
        
        float dTrk1 = sqrt(treeVars.ele_dEtaTrk[id1]*treeVars.ele_dEtaTrk[id1] + treeVars.ele_dPhiTrk[id1]*treeVars.ele_dPhiTrk[id1]);
        float dTrk2 = sqrt(treeVars.ele_dEtaTrk[id2]*treeVars.ele_dEtaTrk[id2] + treeVars.ele_dPhiTrk[id2]*treeVars.ele_dPhiTrk[id2]);
        
        bool passDiEleSelection = DiEleSelections(ele1, ele2, treeVars.ele_charge[id1], treeVars.ele_charge[id2], ph1, ph2, treeVars.ele_miniIso[id1], treeVars.ele_miniIso[id2], dTrk1, dTrk2);
        if( !passDiEleSelection ) continue;
        
        std::pair<int,int> leptons(goodLeptons[l1Count],goodLeptons[l2Count]);
        selectedPairs.push_back(leptons);
      }
      
      else //mixed
      {
        int muIndex = id1;
        int eleIndex = id2;
        if( goodLeptons[l1Count] >= nLep && goodLeptons[l2Count] < nLep )
        {
          muIndex = id2;
          eleIndex = id1;
        }
        
        TLorentzVector mu1, ele1;
        mu1.SetPtEtaPhiE(treeVars.mu_pt[muIndex], treeVars.mu_eta[muIndex], treeVars.mu_phi[muIndex], treeVars.mu_energy[muIndex]);
        ele1.SetPtEtaPhiE(treeVars.ele_pt[eleIndex], treeVars.ele_eta[eleIndex], treeVars.ele_phi[eleIndex], treeVars.ele_energy[eleIndex]);
        
        bool passMixedSelection = MixedSelections(mu1, ele1, treeVars.mu_charge[muIndex], treeVars.ele_charge[eleIndex], ph1, ph2);
        if( !passMixedSelection || !treeVars.mu_IDVector[muIndex][mixedMuID] || !treeVars.ele_IDVector[eleIndex][mixedEleID] ) continue;
            
        std::pair<int, int> leptons(goodLeptons[l1Count], goodLeptons[l2Count]);
        selectedPairs.push_back(leptons);
      }
    }
  }
  
  bool accept = false;
  
  if( verbosity )
    std::cout << "selectedPairs: " << selectedPairs.size() << std::endl;
  
  for(unsigned int j = 0; j < selectedPairs.size(); ++j)
  {
    int id1 = (selectedPairs[j].first)%nLep;
    int id2 = (selectedPairs[j].second)%nLep;
    
    int njet = 0;
    int nbjet_loose = 0;
    int nbjet_medium = 0;
    int nbjet_tight = 0;
    float jetPtMax = -999.;
    
    for(int jIndex = 0; jIndex < nJet; ++jIndex)
    {
      if( treeVars.jet_pt[jIndex] > 25. && fabs(treeVars.jet_eta[jIndex]) < 2.4 )
      {
        float DR1 = -999.;
        float DR2 = -999.;
        if( selectedPairs[j].first < nLep && selectedPairs[j].second < nLep )
        {
          DR1 = DeltaR(treeVars.jet_eta[jIndex], treeVars.jet_phi[jIndex], treeVars.mu_eta[id1], treeVars.mu_phi[id1]);
          DR2 = DeltaR(treeVars.jet_eta[jIndex], treeVars.jet_phi[jIndex], treeVars.mu_eta[id2], treeVars.mu_phi[id2]);
          cat = DiMuon;
        }
        else if( selectedPairs[j].first >= nLep && selectedPairs[j].second >= nLep )
        {
          DR1 = DeltaR(treeVars.jet_eta[jIndex], treeVars.jet_phi[jIndex], treeVars.ele_eta[id1], treeVars.ele_phi[id1]);
          DR2 = DeltaR(treeVars.jet_eta[jIndex], treeVars.jet_phi[jIndex], treeVars.ele_eta[id2], treeVars.ele_phi[id2]);
          cat = DiElectron;
        }
        else
        {
          int muIndex = id1;
          int eleIndex = id2;
          if( selectedPairs[j].first >= nLep && selectedPairs[j].second < nLep )
          {
            muIndex = id2;
            eleIndex = id1;
          }
          
          DR1 = DeltaR(treeVars.jet_eta[jIndex], treeVars.jet_phi[jIndex], treeVars.mu_eta[muIndex], treeVars.mu_phi[muIndex]);
          DR2 = DeltaR(treeVars.jet_eta[jIndex], treeVars.jet_phi[jIndex], treeVars.ele_eta[eleIndex], treeVars.ele_phi[eleIndex]);
          cat = Mixed;
        }
        
        if( verbosity )
          std::cout << "jet " << jIndex << ":   pt: " << treeVars.jet_pt[jIndex] << "   DR1 (>0.4): " << DR1 << "   DR2 (>0.4): " << DR2 << std::endl;
        
        if( DR1 < 0.4 || DR2 < 0.4 ) continue;
        
        ++njet;
        
        if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdLooseDeep  ) ++nbjet_loose;
        if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdMediumDeep ) ++nbjet_medium;
        if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdTightDeep  ) ++nbjet_tight;
        
        if( treeVars.jet_pt[jIndex] > jetPtMax ) jetPtMax = treeVars.jet_pt[jIndex];        
      }
    }
    
    if( type != -2 )
    {
      if( verbosity )
        std::cout << "nJet (>1): " << njet << "   jetPtMax (>50): " << jetPtMax << "   nbjet_medium (>1): " <<nbjet_medium << std::endl;
      
      if( njet >= 1 && jetPtMax > 50. && ( (bTagSelection && nbjet_medium >= 1) || (!bTagSelection && nbjet_medium == 0) ) )
      {
        accept = true;
        break;
      }
    }
    else if( csType == kInvertBTag )
    {
      if( njet >= 1 && jetPtMax > 50. && ( (bTagSelection && nbjet_medium < 1) || (!bTagSelection && nbjet_loose == 0) ) )
      {
        accept = true;
        break;
      }    
    }
  }
  
  return accept;
}



bool SingleLeptonSelection(const TreeVars& treeVars, const bool& useDeepBDisc, const int& type,
                           const int& nJetMin, const float& leadJetPtMin, const int& nBTagLooseJetMin, const int& nBTagMediumJetMin, const int& nBTagTightJetMin,
                           const ControlSampleType& csType, const bool& verbosity,
                           std::vector<int>* goodMu, std::vector<int>* goodEle, std::vector<int>* goodJet)

{
  if( csType == kInvertLeptons && type == -2 )
  {
    int njet = 0;
    int nbjet_loose = 0;
    int nbjet_medium = 0;
    int nbjet_tight = 0;
    float jetPtMax = -999.;
    
    for(int jIndex = 0; jIndex < nJet; ++jIndex)
    {
      if( treeVars.jet_pt[jIndex] > 25. && fabs(treeVars.jet_eta[jIndex]) < 2.4 )
      {
        ++njet;
        if( useDeepBDisc )
        {
          if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdLooseDeep  ) ++nbjet_loose;
          if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdMediumDeep ) ++nbjet_medium;
          if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdTightDeep  ) ++nbjet_tight;
        }
        else
        {
          if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdLooseCSV  ) ++nbjet_loose;
          if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdMediumCSV ) ++nbjet_medium;
          if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdTightCSV  ) ++nbjet_tight;        
        }
        
        if( treeVars.jet_pt[jIndex] > jetPtMax ) jetPtMax = treeVars.jet_pt[jIndex];
      }
    }
    
    if( njet >= nJetMin && jetPtMax > leadJetPtMin && nbjet_loose >= nBTagLooseJetMin && nbjet_medium >= nBTagMediumJetMin && nbjet_tight >= nBTagTightJetMin )
      return true;
    else
      return false;
  }
  
  
  // define the diphoton
  TLorentzVector ph1, ph2;
  ph1.SetPtEtaPhiE(treeVars.dipho_leadPt, treeVars.dipho_leadEta, treeVars.dipho_leadPhi, treeVars.dipho_leadEnergy);
  ph2.SetPtEtaPhiE(treeVars.dipho_subleadPt, treeVars.dipho_subleadEta, treeVars.dipho_subleadPhi, treeVars.dipho_subleadEnergy);
  
  
  // select the leptons
  std::vector<int> goodLeptons;
  
  for(int ii = 0; ii < nLep; ++ii)
  {
    if( treeVars.mu_IDVector[ii][singleMuID] != 1 ) continue;
    
    TLorentzVector mu;
    mu.SetPtEtaPhiE(treeVars.mu_pt[ii], treeVars.mu_eta[ii], treeVars.mu_phi[ii], treeVars.mu_energy[ii]);
    
    float iso = ( treeVars.mu_sumChargedHadronPt[ii] + std::max(0.,treeVars.mu_sumNeutralHadronEt[ii]+treeVars.mu_sumPhotonEt[ii]-0.5*treeVars.mu_sumPUPt[ii]) ) / treeVars.mu_pt[ii];
    
    if( SingleMuSelections(mu, ph1, ph2, iso) )
    {
      goodLeptons.push_back(ii);
      if( goodMu ) goodMu -> push_back(ii);
    }
  }
  
  for(int ii = 0; ii < nLep; ++ii)
  {
    if( treeVars.ele_IDVector[ii][singleEleID] != 1 ) continue;
    
    TLorentzVector ele;
    ele.SetPtEtaPhiE(treeVars.ele_pt[ii], treeVars.ele_eta[ii], treeVars.ele_phi[ii], treeVars.ele_energy[ii]);
    
    float dTrk = sqrt(treeVars.ele_dEtaTrk[ii]*treeVars.ele_dEtaTrk[ii] + treeVars.ele_dPhiTrk[ii]*treeVars.ele_dPhiTrk[ii]);
    
    if( SingleEleSelections(ele, ph1, ph2, treeVars.ele_miniIso[ii], dTrk) )
    {
      goodLeptons.push_back(nLep+ii);
      if( goodEle ) goodEle -> push_back(ii);
    }
  }
  
  if( verbosity )
    for(unsigned int ii = 0; ii < goodLeptons.size(); ++ii)
      std::cout << "goodLeptons[" << ii << "] = " << goodLeptons.at(ii) << std::endl;
  
  
  bool accept = false;
  
  for(unsigned int j = 0; j < goodLeptons.size(); ++j)
  { 
    int njet = 0;
    int nbjet_loose = 0;
    int nbjet_medium = 0;
    int nbjet_tight = 0;
    float jetPtMax = -999.;
    
    for(int jIndex = 0; jIndex < nJet; ++jIndex)
    {
      if( treeVars.jet_pt[jIndex] > 25. && abs(treeVars.jet_eta[jIndex]) < 2.4 )
      {
        float DR = -1.;
        if( goodLeptons.at(j) < nLep ) DR = DeltaR(treeVars.jet_eta[jIndex], treeVars.jet_phi[jIndex], treeVars.mu_eta[goodLeptons.at(j)], treeVars.mu_phi[goodLeptons.at(j)]);
        else                           DR = DeltaR(treeVars.jet_eta[jIndex], treeVars.jet_phi[jIndex], treeVars.ele_eta[goodLeptons.at(j)-nLep], treeVars.ele_phi[goodLeptons.at(j)-nLep]);
        
        if( verbosity )
          std::cout << "jet " << jIndex << ":   pt: " << treeVars.jet_pt[jIndex] << "   DR (>0.4): " << DR << std::endl;
        
        if( DR < 0.4 ) continue;
        
        ++njet;
        if( useDeepBDisc )
        {
          if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdLooseDeep  ) ++nbjet_loose;
          if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdMediumDeep ) ++nbjet_medium;
          if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdTightDeep  ) ++nbjet_tight;
        }
        else
        {
          if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdLooseCSV  ) ++nbjet_loose;
          if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdMediumCSV ) ++nbjet_medium;
          if( treeVars.jet_bdiscriminant[jIndex] > bDiscriminantThresholdTightCSV  ) ++nbjet_tight;          
        }
        
        if( treeVars.jet_pt[jIndex] > jetPtMax ) jetPtMax = treeVars.jet_pt[jIndex];
        
        if( goodJet ) goodJet -> push_back(jIndex);
      }
    }
    
    
    if( type != -2 )
    {
      if( verbosity )
        std::cout << "nJet: " << njet << "   jetPtMax: " << jetPtMax << "   nbjet_medium: " <<nbjet_medium << std::endl;
      
      if( njet >= nJetMin && jetPtMax > leadJetPtMin && nbjet_loose >= nBTagLooseJetMin && nbjet_medium >= nBTagMediumJetMin && nbjet_tight >= nBTagTightJetMin)
      {
        accept = true;
        break;
      }
    }
    else if( csType == kInvertBTag )
    {
      if( njet >= nJetMin && jetPtMax > leadJetPtMin && (nbjet_loose < nBTagLooseJetMin || nbjet_medium < nBTagMediumJetMin || nbjet_tight < nBTagTightJetMin) )
      {
        accept = true;
        break;
      }    
    }
  }
  
  return accept;
}



bool CutBasedSelection(const TreeVars& treeVars,
                       const float& min_lead_ptoM, const float& min_sublead_ptoM,
                       const float& min_leadIDMVA, const float& min_subleadIDMVA,
                       const float& max_deltaphi, const float& max_deltaeta)
{
  if( treeVars.dipho_lead_ptoM < min_lead_ptoM || treeVars.dipho_sublead_ptoM < min_sublead_ptoM ) return false;
  if( treeVars.dipho_leadIDMVA < min_leadIDMVA || treeVars.dipho_subleadIDMVA < min_subleadIDMVA ) return false;
  if( treeVars.dipho_deltaphi > max_deltaphi ) return false;
  if( fabs( treeVars.dipho_leadEta - treeVars.dipho_subleadEta ) > max_deltaeta ) return false;
  
  return true;
}
