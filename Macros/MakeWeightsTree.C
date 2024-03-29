#include <TChain.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>

using namespace std;

/**
 * Takes in a file path and file pattern, and "adds" all the TTrees together in a TChain and slims it down for generating event-level weights.
 */
void MakeWeightsTree (const char* path, const char* PbPbFilePattern, const char* ppFilePattern, const char* outFileName = "eventWeightsTree.root") {

  TChain* inTree = new TChain ("bush", "bush");

  {
    auto dir = gSystem->OpenDirectory (path);
    while (auto f = gSystem->GetDirEntry (dir)) {
      if (!strcmp (f, ".") || !strcmp (f, ".."))
        continue;
      if (string (f).find (PbPbFilePattern) != string::npos)
        inTree->Add (TString (path) + "/" + f + "/*.root");
    }
    gSystem->FreeDirectory (dir);
  }

  TFile* outFile = new TFile (Form ("%s/%s", path, outFileName), "recreate");
  outFile->Delete ("PbPbZTrackTree;*");
  outFile->Delete ("ppZTrackTree;*");

  unsigned int event_number = 0, run_number = 0, lumi_block = 0;
  bool passes_toroid = false, isOOTPU = false, BlayerDesyn = false;
  float fcalA_et = 0, fcalC_et = 0;
  float fcalA_et_Cos2 = 0, fcalC_et_Cos2 = 0, fcalA_et_Sin2 = 0, fcalC_et_Sin2 = 0;
  float fcalA_et_Cos3 = 0, fcalC_et_Cos3 = 0, fcalA_et_Sin3 = 0, fcalC_et_Sin3 = 0;
  float fcalA_et_Cos4 = 0, fcalC_et_Cos4 = 0, fcalA_et_Sin4 = 0, fcalC_et_Sin4 = 0;
  vector<float>* mcEventWeights = nullptr;
  int ntrk = 0;
  float trk_pt[10000];
  bool trk_HIloose[10000];
  bool trk_HItight[10000];

  inTree->SetBranchAddress ("event_number", &event_number);
  inTree->SetBranchAddress ("run_number", &run_number);
  inTree->SetBranchAddress ("lumi_block", &lumi_block);
  inTree->SetBranchAddress ("passes_toroid", &passes_toroid);
  inTree->SetBranchAddress ("isOOTPU", &isOOTPU);
  inTree->SetBranchAddress ("BlayerDesyn", &BlayerDesyn);
  inTree->SetBranchAddress ("mcEventWeights", &mcEventWeights);
  inTree->SetBranchAddress ("fcalA_et", &fcalA_et);
  inTree->SetBranchAddress ("fcalC_et", &fcalC_et);
  inTree->SetBranchAddress ("fcalA_et_Cos2", &fcalA_et_Cos2);
  inTree->SetBranchAddress ("fcalC_et_Cos2", &fcalC_et_Cos2);
  inTree->SetBranchAddress ("fcalA_et_Sin2", &fcalA_et_Sin2);
  inTree->SetBranchAddress ("fcalC_et_Sin2", &fcalC_et_Sin2);
  inTree->SetBranchAddress ("fcalA_et_Cos3", &fcalA_et_Cos3);
  inTree->SetBranchAddress ("fcalC_et_Cos3", &fcalC_et_Cos3);
  inTree->SetBranchAddress ("fcalA_et_Sin3", &fcalA_et_Sin3);
  inTree->SetBranchAddress ("fcalC_et_Sin3", &fcalC_et_Sin3);
  inTree->SetBranchAddress ("fcalA_et_Cos4", &fcalA_et_Cos4);
  inTree->SetBranchAddress ("fcalC_et_Cos4", &fcalC_et_Cos4);
  inTree->SetBranchAddress ("fcalA_et_Sin4", &fcalA_et_Sin4);
  inTree->SetBranchAddress ("fcalC_et_Sin4", &fcalC_et_Sin4);
  inTree->SetBranchAddress ("ntrk", &ntrk);
  inTree->SetBranchAddress ("trk_pt", &trk_pt);
  inTree->SetBranchAddress ("trk_HIloose", trk_HIloose);
  inTree->SetBranchAddress ("trk_HItight", trk_HItight);

  outFile->cd ();
  TTree* outTree = new TTree ("PbPbZTrackTree", "PbPbZTrackTree");

  float fcal_et = 0, q2 = 0, psi2 = 0, q3 = 0, psi3 = 0, q4 = 0, psi4 = 0, event_weight = 1;
  int ntrk_cut = 0;

  outTree->Branch ("event_number", &event_number);
  outTree->Branch ("run_number", &run_number);
  outTree->Branch ("lumi_block", &lumi_block);
  outTree->Branch ("fcal_et", &fcal_et);
  outTree->Branch ("q2", &q2);
  outTree->Branch ("psi2", &psi2);
  outTree->Branch ("q3", &q3);
  outTree->Branch ("psi3", &psi3);
  outTree->Branch ("q4", &q4);
  outTree->Branch ("psi4", &psi4);
  outTree->Branch ("event_weight", &event_weight);
  outTree->Branch ("ntrk_cut", &ntrk_cut);
  outTree->Branch ("ntrk_all", &ntrk);

  for (int iEntry = 0; iEntry < inTree->GetEntries (); iEntry++) {
    if (inTree->GetEntries () > 100 && iEntry % (inTree->GetEntries () / 100) == 0)
      cout << iEntry / (inTree->GetEntries () / 100) << "\% done...\r" << flush;
    inTree->GetEntry (iEntry);
    if (BlayerDesyn)
      continue;
    if (isOOTPU)
      continue;

    fcal_et = fcalA_et+fcalC_et;
    q2 = (fcal_et > 0 ? sqrt (pow (fcalA_et_Cos2+fcalC_et_Cos2, 2) + pow (fcalA_et_Sin2+fcalC_et_Sin2, 2)) / (fcal_et) : 0);
    psi2 = atan2 (fcalA_et_Sin2+fcalC_et_Sin2, fcalA_et_Cos2+fcalC_et_Cos2) / 2.;
    q3 = (fcal_et > 0 ? sqrt (pow (fcalA_et_Cos3+fcalC_et_Cos3, 2) + pow (fcalA_et_Sin3+fcalC_et_Sin3, 2)) / (fcal_et) : 0);
    psi3 = atan2 (fcalA_et_Sin3+fcalC_et_Sin3, fcalA_et_Cos3+fcalC_et_Cos3) / 3.;
    q4 = (fcal_et > 0 ? sqrt (pow (fcalA_et_Cos4+fcalC_et_Cos4, 2) + pow (fcalA_et_Sin4+fcalC_et_Sin4, 2)) / (fcal_et) : 0);
    psi4 = atan2 (fcalA_et_Sin4+fcalC_et_Sin4, fcalA_et_Cos4+fcalC_et_Cos4) / 4.;

    if (mcEventWeights)
      event_weight = mcEventWeights->at (0);
    else
      event_weight = 1;

    ntrk_cut = 0;
    for (int iTrk = 0; iTrk < ntrk; iTrk++) {
      if (trk_pt[iTrk] < 1 || !trk_HIloose[iTrk]) continue;
      ntrk_cut++;
    }
    outTree->Fill ();
  }
  cout << "Done loop #1" << endl;

  outTree->SetDirectory (outFile);
  outTree->Write ("", TObject :: kOverwrite);


  inTree->Reset ();
  {
    auto dir = gSystem->OpenDirectory (path);
    while (auto f = gSystem->GetDirEntry (dir)) {
      if (!strcmp (f, ".") || !strcmp (f, ".."))
        continue;
      if (string (f).find (ppFilePattern) != string::npos)
        inTree->Add (TString (path) + "/" + f + "/*.root");
    }
    gSystem->FreeDirectory (dir);
  }

  inTree->SetBranchAddress ("event_number", &event_number);
  inTree->SetBranchAddress ("run_number", &run_number);
  inTree->SetBranchAddress ("lumi_block", &lumi_block);
  inTree->SetBranchAddress ("passes_toroid", &passes_toroid);
  inTree->SetBranchAddress ("isOOTPU", &isOOTPU);
  inTree->SetBranchAddress ("BlayerDesyn", &BlayerDesyn);
  inTree->SetBranchAddress ("mcEventWeights", &mcEventWeights);
  inTree->SetBranchAddress ("fcalA_et", &fcalA_et);
  inTree->SetBranchAddress ("fcalC_et", &fcalC_et);
  inTree->SetBranchAddress ("fcalA_et_Cos2", &fcalA_et_Cos2);
  inTree->SetBranchAddress ("fcalC_et_Cos2", &fcalC_et_Cos2);
  inTree->SetBranchAddress ("fcalA_et_Sin2", &fcalA_et_Sin2);
  inTree->SetBranchAddress ("fcalC_et_Sin2", &fcalC_et_Sin2);
  inTree->SetBranchAddress ("fcalA_et_Cos3", &fcalA_et_Cos3);
  inTree->SetBranchAddress ("fcalC_et_Cos3", &fcalC_et_Cos3);
  inTree->SetBranchAddress ("fcalA_et_Sin3", &fcalA_et_Sin3);
  inTree->SetBranchAddress ("fcalC_et_Sin3", &fcalC_et_Sin3);
  inTree->SetBranchAddress ("fcalA_et_Cos4", &fcalA_et_Cos4);
  inTree->SetBranchAddress ("fcalC_et_Cos4", &fcalC_et_Cos4);
  inTree->SetBranchAddress ("fcalA_et_Sin4", &fcalA_et_Sin4);
  inTree->SetBranchAddress ("fcalC_et_Sin4", &fcalC_et_Sin4);
  inTree->SetBranchAddress ("ntrk", &ntrk);
  inTree->SetBranchAddress ("trk_pt", &trk_pt);
  inTree->SetBranchAddress ("trk_HIloose", trk_HIloose);
  inTree->SetBranchAddress ("trk_HItight", trk_HItight);

  outFile->cd ();
  outTree = new TTree ("ppZTrackTree", "ppZTrackTree");

  outTree->Branch ("event_number", &event_number);
  outTree->Branch ("run_number", &run_number);
  outTree->Branch ("lumi_block", &lumi_block);
  outTree->Branch ("fcal_et", &fcal_et);
  outTree->Branch ("q2", &q2);
  outTree->Branch ("psi2", &psi2);
  outTree->Branch ("q3", &q3);
  outTree->Branch ("psi3", &psi3);
  outTree->Branch ("q4", &q4);
  outTree->Branch ("psi4", &psi4);
  outTree->Branch ("event_weight", &event_weight);
  outTree->Branch ("ntrk_cut", &ntrk_cut);
  outTree->Branch ("ntrk_all", &ntrk);

  for (int iEntry = 0; iEntry < inTree->GetEntries (); iEntry++) {
    if (inTree->GetEntries () > 100 && iEntry % (inTree->GetEntries () / 100) == 0)
      cout << iEntry / (inTree->GetEntries () / 100) << "\% done...\r" << flush;
    inTree->GetEntry (iEntry);
    if (BlayerDesyn)
      continue;
    if (isOOTPU)
      continue;

    fcal_et = fcalA_et+fcalC_et;
    q2 = (fcal_et > 0 ? sqrt (pow (fcalA_et_Cos2+fcalC_et_Cos2, 2) + pow (fcalA_et_Sin2+fcalC_et_Sin2, 2)) / (fcal_et) : 0);
    psi2 = atan2 (fcalA_et_Sin2+fcalC_et_Sin2, fcalA_et_Cos2+fcalC_et_Cos2) / 2.;
    q3 = (fcal_et > 0 ? sqrt (pow (fcalA_et_Cos3+fcalC_et_Cos3, 2) + pow (fcalA_et_Sin3+fcalC_et_Sin3, 2)) / (fcal_et) : 0);
    psi3 = atan2 (fcalA_et_Sin3+fcalC_et_Sin3, fcalA_et_Cos3+fcalC_et_Cos3) / 3.;
    q4 = (fcal_et > 0 ? sqrt (pow (fcalA_et_Cos4+fcalC_et_Cos4, 2) + pow (fcalA_et_Sin4+fcalC_et_Sin4, 2)) / (fcal_et) : 0);
    psi4 = atan2 (fcalA_et_Sin4+fcalC_et_Sin4, fcalA_et_Cos4+fcalC_et_Cos4) / 4.;

    if (mcEventWeights)
      event_weight = mcEventWeights->at (0);
    else
      event_weight = 1;

    ntrk_cut = 0;
    for (int iTrk = 0; iTrk < ntrk; iTrk++) {
      if (trk_pt[iTrk] < 1 || !trk_HIloose[iTrk]) continue;
      ntrk_cut++;
    }
    outTree->Fill ();
  }
  cout << "Done loop #2" << endl;

  outTree->SetDirectory (outFile);
  outTree->Write ("", TObject :: kOverwrite);

  outFile->Close ();
}
  


void MakeppMCWeightsTree () {
  MakeWeightsTree ("/atlasgpfs01/usatlas/data/jeff/ZHadronAnalysis/data/mc_100", ".PbPb.", ".pp.", "ppEventWeightsTree.root");
}

void MakePbPbMCWeightsTree () {
  MakeWeightsTree ("/atlasgpfs01/usatlas/data/jeff/ZHadronAnalysis/data/mc_100", ".PbPb.", ".pp.", "PbPbEventWeightsTree.root");
}

void MakeHijingWeightsTree () {
  MakeWeightsTree ("/atlasgpfs01/usatlas/data/jeff/ZHadronAnalysis/data/hijing_200", ".PbPb.", ".pp.", "PbPbEventWeightsTree.root");
}

//void Make2015HijingWeightsTree () {
//  MakeWeightsTree ("/atlasgpfs01/usatlas/data/jeff/ZHadronAnalysis/data/hijing_110", "user.jeouelle.21.2.82.minbias.110.mc16_5TeV.PbPb.420000.Hijing.Flow_JJFV6.recon.AOD.e4962_a884_s3136_r11318_myOutput.root", "*.pp.*", "hijing2015_eventWeightsTree.root");
//}

//void Make2018HijingWeightsTree () {
//  MakeWeightsTree ("/atlasgpfs01/usatlas/data/jeff/ZHadronAnalysis/data/hijing_110", "user.jeouelle.21.2.82.minbias.110.mc16_5TeV.PbPb.420000.Hijing.Flow_JJFV6.recon.AOD.e4962_a882_s3136_r11157_myOutput.root", "*.pp.*", "hijing2018_eventWeightsTree.root");
//}
