#ifndef __MinbiasAnalysis_h__
#define __MinbiasAnalysis_h__

#include "Params.h"
#include "PhysicsAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;
using namespace atlashi;

class MinbiasAnalysis : public PhysicsAnalysis {

  public:

  bool takeNonTruthTracks = false;
  bool doQ2Mixing = false;
  bool doPsi2Mixing = false;
  bool doPsi3Mixing = false;

  int numQ2MixBins = 1;
  double* q2MixBins = nullptr;
  int numPsi2MixBins = 1;
  double* psi2MixBins = nullptr;
  int numPsi3MixBins = 1;
  double* psi3MixBins = nullptr;
  

  short GetQ2MixBin (const float q2) {
    if (!q2MixBins)
      return -1;
    short i = 0;
    while (i < numQ2MixBins) {
      if (q2 < q2MixBins[i+1])
        break;
      i++;
    }
    return i;
  }

  short GetPsi2MixBin (const float psi2) {
    if (!psi2MixBins)
      return -1;
    short i = 0;
    while (i < numPsi2MixBins) {
      if (psi2 < psi2MixBins[i+1])
        break;
      i++;
    }
    return i;
  }

  short GetPsi3MixBin (const float psi3) {
    if (!psi3MixBins)
      return -1;
    short i = 0;
    while (i < numPsi3MixBins) {
      if (psi3 < psi3MixBins[i+1])
        break;
      i++;
    }
    return i;
  }

  MinbiasAnalysis (const char* _name = "bkg") : PhysicsAnalysis () {
    name = _name;
    plotFill = true;
    plotSignal = false;
    useAltMarker = false;
    hasBkg = false;
    backgroundSubtracted = true;
    histsUnfolded = true;
    iaaCalculated = true;
    //icpCalculated = true;
  }

  void LoadEventWeights () override;
  void GenerateWeights (const char* inFileName = "outFile.root", const char* outFileName = "eventWeightsFile.root");

  virtual void LoadHists (const char* histFileName = "savedHists.root", const bool _finishHists = true) override;
  void Execute (const char* inFileName, const char* outFileName) override {
    cout << "Error: In MinbiasAnalysis :: Execute: Called invalid function! 4 arguments are required to specify the collision system and mixing file name. Exiting." << endl;
  }
  void Execute (const bool isPbPb, const char* inFileName, const char* mbInFileName, const char* outFileName);

  void PlotCentralityDists ();
};




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Load histograms into memory, then combine channels.
//////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: LoadHists (const char* histFileName, const bool _finishHists) {
  PhysicsAnalysis :: LoadHists (histFileName, _finishHists);

  PhysicsAnalysis :: CombineHists ();
}




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Load the event weights into memory
//////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: LoadEventWeights () {
  if (eventWeightsLoaded)
    return;

  SetupDirectories ("", "ZTrackAnalysis/");
  TDirectory* _gDirectory = gDirectory;

  eventWeightsFile = new TFile (Form ("%s/MinbiasAnalysis/Nominal/eventWeightsFile.root", rootPath.Data ()), "read");

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_PbPbQ2_weights[iSpc][iFineCent][iPtZ] = (TH1D*) eventWeightsFile->Get (Form ("h_z_q2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()));
        h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ] = (TH1D*) eventWeightsFile->Get (Form ("h_z_psi2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()));
      }
    }
  }

  eventWeightsLoaded = true;

  _gDirectory-> cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over minbias trees and fills histograms appropriately. (NEW VERSION)
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: Execute (const bool isPbPb, const char* inFileName, const char* mbInFileName, const char* outFileName) {

  //LoadEventWeights ();

  SetupDirectories ("", "ZTrackAnalysis/");

  cout << "Arguments provided: " << endl;
  cout << "isPbPb = " << isPbPb << endl;
  cout << "inFileName = " << inFileName << endl;
  cout << "mbInFileName = " << mbInFileName << endl;
  cout << "outFileName = " << outFileName << endl;

  const bool doSameFileMixing = (string (inFileName) == string (mbInFileName));
  if (doSameFileMixing) {
    cout << "Attempting to mix events within the same file!" << endl;
  }

  if (doQ2Mixing) cout << "Attempting to mix events in q2 magnitude with " << numQ2MixBins << " bins" << endl;
  if (doPsi2Mixing) cout << "Attempting to mix events in psi2 angle with " << numPsi2MixBins << " bins" << endl;
  if (doPsi3Mixing) cout << "Attempting to mix events in psi3 angle with " << numPsi3MixBins << " bins" << endl;
  q2MixBins = linspace (0, 0.2, numQ2MixBins);
  psi2MixBins = linspace (-pi/2, pi/2, numPsi2MixBins);
  psi3MixBins = linspace (-pi/3, pi/3, numPsi3MixBins);


  CreateHists ();

  //eventPlaneCalibrator = EventPlaneCalibrator (Form ("%s/FCalCalibration/Nominal/data18hi.root", rootPath.Data ()));


  // Get TTree with mixed event pool
  TFile* mbInFile = new TFile (Form ("%s/%s", rootPath.Data (), mbInFileName), "read");
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), mbInFileName) << endl;
  TTree* mbTree = (TTree*) mbInFile->Get (doSameFileMixing ? (isPbPb ? "PbPbZTrackTree" : "ppZTrackTree") : (isPbPb ? "PbPbTrackTree" : "ppZTrackTree"));
  const int nMBEvts = mbTree->GetEntries ();
  cout << "N events in mixing pool = " << nMBEvts << endl;

  int iMBEvt = 0;
  mbTree->GetEntry (iMBEvt);

  std::vector <int> mbEventOrder = {};
  std::vector <bool> mbEventsUsed = {};
  for (int i = 0; i < nMBEvts; i++) {
    mbEventOrder.push_back (i);
    mbEventsUsed.push_back (false);
  }
  if (!doSameFileMixing) {
    std::srand (1);
    std::random_shuffle (mbEventOrder.begin (), mbEventOrder.end ());

    cout << "Shuffled minimum bias events randomly; first ten events are ";
    for (int i = 0; i < 10; i++) {
      cout << mbEventOrder[i];
      if (i < 9)
        cout << ", ";
    }
    cout << endl;
  }


  // Get TTree with Z bosons
  cout << Form ("Reading Z bosons from %s/%s", rootPath.Data (), inFileName) << endl;
  TFile* inFile = nullptr;
  TTree* zTree = nullptr;
  if (doSameFileMixing)
    zTree = mbTree;
  else {
    inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
    zTree = (TTree*) inFile->Get (isPbPb ? "PbPbZTrackTree" : "ppZTrackTree");
  }
  if (!zTree)
    cout << "Got a null Z boson tree!" << endl;

  int nZEvts = zTree->GetEntries ();

  if (nZEvts == 0)
    cout << "Warning! No Z's to mix with in this run!" << endl;
  cout << "N Z-tagged events = " << nZEvts << endl;
  cout << "For this Z tree, maximum mixing fraction = " << nMBEvts / nZEvts << endl;

  bool doShuffle = false;
  if (mixingFraction * nZEvts > nMBEvts && !doSameFileMixing) {
    cout << "Warning! Mixing fraction too high, will use " << (float)(nMBEvts / mixingFraction) / (float)(nZEvts) * 100. << "% of Z events" << endl;
    nZEvts = nMBEvts / mixingFraction;
    doShuffle = true;
  }

  std::vector <int> zEventOrder = {};
  std::vector <int> zEventsUsed = {};
  for (int i = 0; i < nZEvts; i++) {
    zEventOrder.push_back (i);
    zEventsUsed.push_back (0);
  }

  if (doShuffle) {
    std::srand (1);
    std::random_shuffle (zEventOrder.begin (), zEventOrder.end ());

    cout << "Shuffled Z boson events randomly; first ten events are ";
    for (int i = 0; i < 10; i++) {
      cout << mbEventOrder[i];
      if (i < 9)
        cout << ", ";
    }
    cout << endl;
  }


  // output TTree to store mixed event
  TFile* mixedEventsFile = new TFile (Form ("%s/%s_mixed.root", rootPath.Data (), inFileName), "recreate");
  TTree* mixedEventsTree = new TTree (isPbPb ? "PbPbMixedTree" : "ppMixedTree", isPbPb ? "PbPbMixedTree" : "ppMixedTree");


  // variables for branches
  int event_number = 0, z_event_number = 0;
  int run_number = 0, z_run_number = 0;
  bool isEE = false;//, passes_toroid = false;
  bool _isEE = false;
  float event_weight = 1, z_event_weight = 1;// q2_weight = 1, psi2_weight = 1; // vz_weight = 1, nch_weight = 1;
  float fcal_et = 0, vz = 0, zdcEnergy = 0, z_fcal_et = 0, z_vz = 0, z_zdcEnergy = 0;
  float phi_transmin = 0, phi_transmax = 0;
  //float q2x_a = 0, q2y_a = 0, q2x_c = 0, q2y_c = 0, z_q2x_a = 0, z_q2y_a = 0, z_q2x_c = 0, z_q2y_c = 0;
  float q2 = 0., q3 = 0., q4 = 0., z_q2 = 0., z_q3 = 0., z_q4 = 0.;
  float psi2 = 0., psi3 = 0., psi4 = 0., z_psi2 = 0., z_psi3 = 0., z_psi4 = 0.;
  float _z_pt = 0, _z_y = 0, _z_phi = 0, _z_m = 0;
  float z_pt = 0, z_y = 0, z_phi = 0, z_m = 0;
  float l1_pt = 0, l1_eta = 0, l1_phi = 0, l1_trk_pt = 0, l1_trk_eta = 0, l1_trk_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0, l2_trk_pt = 0, l2_trk_eta = 0, l2_trk_phi = 0;
  float _l1_pt = 0, _l1_eta = 0, _l1_phi = 0, _l1_trk_pt = 0, _l1_trk_eta = 0, _l1_trk_phi = 0, _l2_pt = 0, _l2_eta = 0, _l2_phi = 0, _l2_trk_pt = 0, _l2_trk_eta = 0, _l2_trk_phi = 0;
  int l1_charge = 0, l2_charge = 0;
  int _l1_charge = 0, _l2_charge = 0;
  int ntrk = 0, z_ntrk = 0;
  float trk_pt[10000], trk_eta[10000], trk_phi[10000];
  bool trk_truth_matched[10000];
  float z_trk_pt[10000], z_trk_eta[10000], z_trk_phi[10000];
  bool HLT_mb_sptrk_L1ZDC_A_C_VTE50 = false;
  bool HLT_noalg_pc_L1TE50_VTE600_0ETA49 = false;
  bool HLT_noalg_cc_L1TE600_0ETA49 = false;

  int**   trks_counts   = Get2DArray <int> (2, 6);
  float** trks_weights1 = Get2DArray <float> (2, 6);
  float** trks_weights2 = Get2DArray <float> (2, 6);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Do this if TTree is PbPb
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (isPbPb) {
    mbTree->LoadBaskets (4000000000); //2000000000 = 2GB
    mbTree->SetBranchAddress ("run_number",   &run_number);
    mbTree->SetBranchAddress ("event_number", &event_number);
    mbTree->SetBranchAddress ("event_weight", &event_weight);
    mbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    mbTree->SetBranchAddress ("zdcEnergy",    &zdcEnergy);
    //mbTree->SetBranchAddress ("q2x_a",         &q2x_a);
    //mbTree->SetBranchAddress ("q2y_a",         &q2y_a);
    //mbTree->SetBranchAddress ("q2x_c",         &q2x_c);
    //mbTree->SetBranchAddress ("q2y_c",         &q2y_c);
    mbTree->SetBranchAddress ("q2",           &q2);
    mbTree->SetBranchAddress ("psi2",         &psi2);
    mbTree->SetBranchAddress ("q3",           &q3);
    mbTree->SetBranchAddress ("psi3",         &psi3);
    mbTree->SetBranchAddress ("q4",           &q4);
    mbTree->SetBranchAddress ("psi4",         &psi4);
    mbTree->SetBranchAddress ("vz",           &vz);
    mbTree->SetBranchAddress ("ntrk",         &ntrk);
    mbTree->SetBranchAddress ("trk_pt",       &trk_pt);
    mbTree->SetBranchAddress ("trk_eta",      &trk_eta);
    mbTree->SetBranchAddress ("trk_phi",      &trk_phi);
    mbTree->SetBranchAddress ("HLT_mb_sptrk_L1ZDC_A_C_VTE50",       &HLT_mb_sptrk_L1ZDC_A_C_VTE50);
    mbTree->SetBranchAddress ("HLT_noalg_pc_L1TE50_VTE600.0ETA49",  &HLT_noalg_pc_L1TE50_VTE600_0ETA49);
    mbTree->SetBranchAddress ("HLT_noalg_cc_L1TE600_0ETA49",        &HLT_noalg_cc_L1TE600_0ETA49);

    if (takeNonTruthTracks)
      mbTree->SetBranchAddress ("trk_truth_matched", &trk_truth_matched);


    zTree->LoadBaskets (2000000000);
    if (!doSameFileMixing) {
      zTree->SetBranchAddress ("run_number",    &z_run_number);
      zTree->SetBranchAddress ("event_number",  &z_event_number);
      zTree->SetBranchAddress ("event_weight",  &z_event_weight);
      zTree->SetBranchAddress ("fcal_et",       &z_fcal_et);
      zTree->SetBranchAddress ("zdcEnergy",     &z_zdcEnergy);
      //zTree->SetBranchAddress ("q2x_a",          &z_q2x_a);
      //zTree->SetBranchAddress ("q2y_a",          &z_q2y_a);
      //zTree->SetBranchAddress ("q2x_c",          &z_q2x_c);
      //zTree->SetBranchAddress ("q2y_c",          &z_q2y_c);
      zTree->SetBranchAddress ("q2",            &z_q2);
      zTree->SetBranchAddress ("psi2",          &z_psi2);
      zTree->SetBranchAddress ("q3",            &z_q3);
      zTree->SetBranchAddress ("psi3",          &z_psi3);
      zTree->SetBranchAddress ("q4",            &z_q4);
      zTree->SetBranchAddress ("psi4",          &z_psi4);
      zTree->SetBranchAddress ("vz",            &z_vz);
      zTree->SetBranchAddress ("ntrk",          &z_ntrk);
      zTree->SetBranchAddress ("trk_pt",        &z_trk_pt);
      zTree->SetBranchAddress ("trk_eta",       &z_trk_eta);
      zTree->SetBranchAddress ("trk_phi",       &z_trk_phi);
    }
    zTree->SetBranchAddress ("isEE",          &_isEE);
    zTree->SetBranchAddress ("z_pt",          &_z_pt);
    zTree->SetBranchAddress ("z_y",           &_z_y);
    zTree->SetBranchAddress ("z_phi",         &_z_phi);
    zTree->SetBranchAddress ("z_m",           &_z_m);
    zTree->SetBranchAddress ("l1_pt",         &_l1_pt);
    zTree->SetBranchAddress ("l1_eta",        &_l1_eta);
    zTree->SetBranchAddress ("l1_phi",        &_l1_phi);
    zTree->SetBranchAddress ("l1_trk_pt",     &_l1_trk_pt);
    zTree->SetBranchAddress ("l1_trk_eta",    &_l1_trk_eta);
    zTree->SetBranchAddress ("l1_trk_phi",    &_l1_trk_phi);
    zTree->SetBranchAddress ("l1_charge",     &_l1_charge);
    zTree->SetBranchAddress ("l2_pt",         &_l2_pt);
    zTree->SetBranchAddress ("l2_eta",        &_l2_eta);
    zTree->SetBranchAddress ("l2_phi",        &_l2_phi);
    zTree->SetBranchAddress ("l2_trk_pt",     &_l2_trk_pt);
    zTree->SetBranchAddress ("l2_trk_eta",    &_l2_trk_eta);
    zTree->SetBranchAddress ("l2_trk_phi",    &_l2_trk_phi);
    zTree->SetBranchAddress ("l2_charge",     &_l2_charge);


    mixedEventsTree->Branch ("run_number",      &run_number,      "run_number/I");
    mixedEventsTree->Branch ("event_number",    &event_number,    "event_number/I");
    mixedEventsTree->Branch ("fcal_et",         &fcal_et,         "fcal_et/F");
    mixedEventsTree->Branch ("zdcEnergy",       &zdcEnergy,       "zdcEnergy/F");
    mixedEventsTree->Branch ("q2",              &q2,              "q2/F");
    mixedEventsTree->Branch ("psi2",            &psi2,            "psi2/F");
    mixedEventsTree->Branch ("vz",              &vz,              "vz/F");
    mixedEventsTree->Branch ("ntrk",            &ntrk,            "ntrk/I");
    mixedEventsTree->Branch ("trk_pt",          &trk_pt,          "trk_pt[ntrk]/F");
    mixedEventsTree->Branch ("trk_eta",         &trk_eta,         "trk_eta[ntrk]/F");
    mixedEventsTree->Branch ("trk_phi",         &trk_phi,         "trk_phi[ntrk]/F");
    mixedEventsTree->Branch ("isEE",            &isEE,            "isEE/O");
    mixedEventsTree->Branch ("z_run_number",    &z_run_number,    "z_run_number/I");
    mixedEventsTree->Branch ("z_event_number",  &z_event_number,  "z_event_number/I");
    mixedEventsTree->Branch ("z_event_weight",  &z_event_weight,  "z_event_weight/F");
    mixedEventsTree->Branch ("z_fcal_et",       &z_fcal_et,       "z_fcal_et/F");
    mixedEventsTree->Branch ("z_zdcEnergy",     &z_zdcEnergy,     "z_zdcEnergy/F");
    mixedEventsTree->Branch ("z_q2",            &z_q2,            "z_q2/F");
    mixedEventsTree->Branch ("z_psi2",          &z_psi2,          "z_psi2/F");
    mixedEventsTree->Branch ("z_vz",            &z_vz,            "z_vz/F");
    mixedEventsTree->Branch ("z_pt",            &z_pt,            "z_pt/F");
    mixedEventsTree->Branch ("z_y",             &z_y,             "z_y/F");
    mixedEventsTree->Branch ("z_phi",           &z_phi,           "z_phi/F");
    mixedEventsTree->Branch ("z_m",             &z_m,             "z_m/F");
    mixedEventsTree->Branch ("l1_pt",           &l1_pt,           "l1_pt/F");
    mixedEventsTree->Branch ("l1_eta",          &l1_eta,          "l1_eta/F");
    mixedEventsTree->Branch ("l1_phi",          &l1_phi,          "l1_phi/F");
    mixedEventsTree->Branch ("l1_trk_pt",       &l1_trk_pt,       "l1_trk_pt/F");
    mixedEventsTree->Branch ("l1_trk_eta",      &l1_trk_eta,      "l1_trk_eta/F");
    mixedEventsTree->Branch ("l1_trk_phi",      &l1_trk_phi,      "l1_trk_phi/F");
    mixedEventsTree->Branch ("l1_charge",       &l1_charge,       "l1_charge/I");
    mixedEventsTree->Branch ("l2_pt",           &l2_pt,           "l2_pt/F");
    mixedEventsTree->Branch ("l2_eta",          &l2_eta,          "l2_eta/F");
    mixedEventsTree->Branch ("l2_phi",          &l2_phi,          "l2_phi/F");
    mixedEventsTree->Branch ("l2_trk_pt",       &l2_trk_pt,       "l2_trk_pt/F");
    mixedEventsTree->Branch ("l2_trk_eta",      &l2_trk_eta,      "l2_trk_eta/F");
    mixedEventsTree->Branch ("l2_trk_phi",      &l2_trk_phi,      "l2_trk_phi/F");
    mixedEventsTree->Branch ("l2_charge",       &l2_charge,       "l2_charge/I");
    mixedEventsTree->Branch ("z_ntrk",          &z_ntrk,          "z_ntrk/I");
    mixedEventsTree->Branch ("z_trk_pt",        &z_trk_pt,        "z_trk_pt[z_ntrk]/F");
    mixedEventsTree->Branch ("z_trk_eta",       &z_trk_eta,       "z_trk_eta[z_ntrk]/F");
    mixedEventsTree->Branch ("z_trk_phi",       &z_trk_phi,       "z_trk_phi[z_ntrk]/F");


    for (int iZEvt = 0; iZEvt < mixingFraction*nZEvts; iZEvt++) {
      if (mixingFraction*nZEvts > 100 && iZEvt % (mixingFraction*nZEvts / 100) == 0)
        cout << iZEvt / (mixingFraction*nZEvts / 100) << "\% done...\r" << flush;

      zTree->GetEntry (zEventOrder[iZEvt % nZEvts]);
      zEventsUsed[iZEvt % nZEvts]++;

      isEE = _isEE;
      z_pt = _z_pt;
      z_y = _z_y;
      z_phi = _z_phi;
      z_m = _z_m;
      l1_pt = _l1_pt;
      l1_eta = _l1_eta;
      l1_phi = _l1_phi;
      l1_trk_pt = _l1_trk_pt;
      l1_trk_eta = _l1_trk_eta;
      l1_trk_phi = _l1_trk_phi;
      l1_charge = _l1_charge;
      l2_pt = _l2_pt;
      l2_eta = _l2_eta;
      l2_phi = _l2_phi;
      l2_trk_pt = _l2_trk_pt;
      l2_trk_eta = _l2_trk_eta;
      l2_trk_phi = _l2_trk_phi;
      l2_charge = _l2_charge;
      if (doSameFileMixing) {
        z_run_number = run_number;
        z_event_number = event_number;
        z_event_weight = event_weight;
        z_fcal_et = fcal_et;
        z_zdcEnergy = zdcEnergy;
        //z_q2x_a = q2x_a;
        //z_q2x_c = q2x_c;
        //z_q2y_a = q2y_a;
        //z_q2y_c = q2y_c;
        z_q2 = q2;
        z_psi2 = psi2;
        z_q3 = q3;
        z_psi3 = psi3;
        z_q4 = q4;
        z_psi4 = psi4;
        z_vz = vz;
        z_ntrk = 0; // placeholder
      }
     

      if (fabs (z_vz) > 150)
        continue;

      if (z_event_weight == 0)
        continue;
      if (isnan (z_event_weight))
        cout << "Warning: z_event_weight = NAN" << endl;

      const short iPtZ = GetPtZBin (z_pt);
      if (iPtZ < 0 || iPtZ > nPtZBins-1)
        continue;

      if (iPtZ < 2)
        continue; // skip unneeded events

      //{
      //  CorrectQ2Vector (z_q2x_a, z_q2y_a, z_q2x_c, z_q2y_c);
      //  const float z_q2x = z_q2x_a + z_q2x_c;
      //  const float z_q2y = z_q2y_a + z_q2y_c;
      //  z_q2 = sqrt (z_q2x*z_q2x + z_q2y*z_q2y) / fcal_et;
      //  z_psi2 = 0.5 * atan2 (z_q2y, z_q2x);
      //}


      // Find the next unused minimum bias event
      {
        const short iFCalEt = GetSuperFineCentBin (z_fcal_et);
        if (iFCalEt < 1 || iFCalEt > numSuperFineCentBins-1)
          continue;
        const short iQ2 = GetQ2MixBin (z_q2);
        if (doQ2Mixing && (iQ2 < 0 || iQ2 > numQ2MixBins-1)) {
          cout << "Out-of-bounds q2, skipping this Z!" << endl;
          continue;
        }
        const short iPsi2 = GetPsi2MixBin (z_psi2);
        if (doPsi2Mixing && (iPsi2 < 0 || iPsi2 > numPsi2MixBins-1)) {
          cout << "Out-of-bounds psi2, skipping this Z!" << endl;
          continue;
        }
        const short iPsi3 = GetPsi3MixBin (z_psi3);
        if (doPsi3Mixing && (iPsi3 < 0 || iPsi3 > numPsi3MixBins-1)) {
          cout << "Out-of-bounds psi3, skipping this Z!" << endl;
          continue;
        }

        // disable big branches during event mixing
        mbTree->SetBranchStatus ("trk_pt", 0);
        mbTree->SetBranchStatus ("trk_eta", 0);
        mbTree->SetBranchStatus ("trk_phi", 0);
        
        bool goodMixEvent = false;
        const int _iMBEvt = iMBEvt;
        do {
          iMBEvt = (iMBEvt+1) % nMBEvts;
          mbTree->GetEntry (mbEventOrder[iMBEvt]);
          //{
          //  CorrectQ2Vector (q2x_a, q2y_a, q2x_c, q2y_c);
          //  const float q2x = q2x_a + q2x_c;
          //  const float q2y = q2y_a + q2y_c;
          //  q2 = sqrt (q2x*q2x + q2y*q2y) / fcal_et;
          //  psi2 = 0.5 * atan2 (q2y, q2x);
          //}
          goodMixEvent = (fabs (vz) < 150 && event_weight != 0); // always require these conditions

          if (!doSameFileMixing || mixingFraction == 1)
            goodMixEvent = (goodMixEvent && !mbEventsUsed[iMBEvt]); // checks for uniqueness (if applicable)

          if (doSameFileMixing)
            goodMixEvent = (goodMixEvent && iMBEvt != iZEvt); // don't mix with the exact same event

          goodMixEvent = (goodMixEvent && iFCalEt == GetSuperFineCentBin (fcal_et)); // always check for centrality matchiing
          if (doQ2Mixing)
            goodMixEvent = (goodMixEvent && iQ2 == GetQ2MixBin (q2)); // potentially match also in q2
          if (doPsi2Mixing)
            goodMixEvent = (goodMixEvent && DeltaPhi (psi2, z_psi2) < (pi / numPsi2MixBins));
            //goodMixEvent = (goodMixEvent && iPsi2 == GetPsi2MixBin (psi2)); // potentially match also in psi2
          if (doPsi3Mixing)
            goodMixEvent = (goodMixEvent && DeltaPhi (psi3, z_psi3) < (2.*pi/3. / numPsi3MixBins));
            //goodMixEvent = (goodMixEvent && iPsi3 == GetPsi3MixBin (psi3)); // potentially match also in psi3

        } while (!goodMixEvent && iMBEvt != _iMBEvt); // only check each event once
        if (_iMBEvt == iMBEvt) {
          cout << "No minbias event to mix with!!! Wrapped around on the same Z!!! Sum Et = " << z_fcal_et << ", q2 = " << z_q2 << ", psi2 = " << z_psi2 << endl;
          continue;
        }
        mbEventsUsed[iMBEvt] = true;

        // reenable big branches after finding a suitable event
        mbTree->SetBranchStatus ("trk_pt", 1);
        mbTree->SetBranchStatus ("trk_eta", 1);
        mbTree->SetBranchStatus ("trk_phi", 1);
        mbTree->GetEntry (mbEventOrder[iMBEvt]);
      }

      // at this point we have a Z boson and a new (unique & random) event to mix with
      mixedEventsTree->Fill ();

      event_weight = z_event_weight;

      const short iSpc = (isEE ? 0 : 1); // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = GetCentBin (z_fcal_et);
      if (iCent < 1 || iCent > numCentBins-1)
        continue;
      const short iFineCent = GetFineCentBin (z_fcal_et);
      if (iFineCent < 1 || iFineCent > numFineCentBins-1)
        continue;

      //{
      //  float dphi = DeltaPhi (z_phi, psi2, false);
      //  if (dphi > pi/2)
      //    dphi = pi - dphi;
      //  //q2_weight = h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->FindBin (q2));
      //  psi2_weight = h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->FindBin (dphi));

      //  //z_event_weight = fcal_weight * q2_weight * psi2_weight * vz_weight;
      //  z_event_weight *= psi2_weight;
      //}
      //z_event_weight = 1;

      {
        const int centrality = GetPercentileCentrality (fcal_et);
        h_fcal_et->Fill (fcal_et);
        h_fcal_et_reweighted->Fill (fcal_et, event_weight);
        if (HLT_mb_sptrk_L1ZDC_A_C_VTE50) {
          h_centrality[0]->Fill (centrality);
          h_centrality_reweighted[0]->Fill (centrality, event_weight);
        }
        if (HLT_noalg_pc_L1TE50_VTE600_0ETA49) {
          h_centrality[1]->Fill (centrality);
          h_centrality_reweighted[1]->Fill (centrality, event_weight);

        }
        if (HLT_noalg_cc_L1TE600_0ETA49) {
          h_centrality[2]->Fill (centrality);
          h_centrality_reweighted[2]->Fill (centrality, event_weight);
        }
        h_q2[iFineCent]->Fill (q2);
        h_q2_reweighted[iFineCent]->Fill (q2, event_weight);
        h_psi2[iFineCent]->Fill (psi2);
        h_psi2_reweighted[iFineCent]->Fill (psi2, event_weight);
        h_PbPb_vz->Fill (vz);
        h_PbPb_vz_reweighted->Fill (vz, event_weight);
      }

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (2.5, pow (event_weight, 2));

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];
        const float xhz = trkpt / z_pt;

        if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
          continue;

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], true);
        const double trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], true);
        if (trkPur == 0 || trkEff == 0)
          continue;
        const float trkWeight = trkPur / trkEff;

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        float dphi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, event_weight * trkWeight);
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt);
            h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, event_weight * trkWeight);
            h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt / z_pt, event_weight * trkWeight);
          }
        }

        if (3*pi/4 <= dphi) {
          if (ptTrkBins[iPtZ][0] <= trkpt) {
            short iPt = 0;
            while (iPt < nPtTrkBins[iPtZ] && ptTrkBins[iPtZ][iPt+1] < trkpt) iPt++;
            if (iPt < 6) {
              trks_counts[0][iPt]   += 1;
              trks_weights1[0][iPt] += trkWeight;
              trks_weights2[0][iPt] += pow (trkWeight, 2);
            }
          }
          if (xHZBins[iPtZ][0] <= xhz) {
            short iX = 0;
            while (iX < nXHZBins[iPtZ] && xHZBins[iPtZ][iX+1] < xhz) iX++;
            if (iX < 6) {
              trks_counts[1][iX]   += 1;
              trks_weights1[1][iX] += trkWeight;
              trks_weights2[1][iX] += pow (trkWeight, 2);
            }
          }
        }
      } // end loop over tracks

      // fill yield histograms and covariance matrices
      for (int i1 = 0; i1 < nPtTrkBins[iPtZ]; i1++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinContent (i1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinContent (i1+1) + event_weight*trks_weights1[0][i1]);
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinError (i1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinError (i1+1) + event_weight*trks_weights1[0][i1]);
        //h_trk_pt_ptz_wgts[iSpc][iPtZ][iCent]->SetBinContent (i1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinContent (i1+1) + event_weight*sqrt(trks_weights2[0][i1]));
        for (int i2 = 0 ; i2 < nPtTrkBins[iPtZ]; i2++)
          h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (i1+1, i2+1, h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (i1+1, i2+1) + event_weight * (trks_weights1[0][i1]) * (trks_weights1[0][i2]));
          //h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (i1+1, i2+1, h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (i1+1, i2+1) + event_weight*(trks_weights1[0][i1])*(trks_weights1[0][i2]));
      } // end loop over i1
      for (int i1 = 0; i1 < nXHZBins[iPtZ]; i1++) {
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinContent (i1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinContent (i1+1) + event_weight*trks_weights1[1][i1]);
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinError (i1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinError (i1+1) + event_weight*trks_weights1[1][i1]);
        //h_trk_xhz_ptz_wgts[iSpc][iPtZ][iCent]->SetBinContent (i1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinContent (i1+1) + event_weight*sqrt (trks_weights2[1][i1]));
        for (int i2 = 0 ; i2 < nXHZBins[iPtZ]; i2++)
          h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (i1+1, i2+1, h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (i1+1, i2+1) + event_weight * (trks_weights1[1][i1]) * (trks_weights1[1][i2]));
          //h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (i1+1, i2+1, h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (i1+1, i2+1) + event_weight*(trks_weights1[1][i1])*(trks_weights1[1][i2]));
      } // end loop over i1

      // reset trk count measurements for next event
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 6; j++) {
          trks_counts[i][j] = 0;
          trks_weights1[i][j] = 0.;
          trks_weights2[i][j] = 0.;
        } // end loop over j
      } // end loop over i

    } // end loop over Pb+Pb tree

    cout << endl << "Done minbias Pb+Pb loop." << endl;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Do this if TTree is pp
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (!isPbPb) {
    assert (doSameFileMixing == doPPMixingVar);

    mbTree->LoadBaskets (3000000000);
    mbTree->SetBranchAddress ("run_number",   &run_number);
    mbTree->SetBranchAddress ("event_number", &event_number);
    mbTree->SetBranchAddress ("event_weight", &event_weight);
    mbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    mbTree->SetBranchAddress ("vz",           &vz);
    mbTree->SetBranchAddress ("ntrk",         &ntrk);
    mbTree->SetBranchAddress ("trk_pt",       &trk_pt);
    mbTree->SetBranchAddress ("trk_eta",      &trk_eta);
    mbTree->SetBranchAddress ("trk_phi",      &trk_phi);

    zTree->LoadBaskets (3000000000);
    if (!doSameFileMixing) {
      zTree->SetBranchAddress ("run_number",    &z_run_number);
      zTree->SetBranchAddress ("event_number",  &z_event_number);
      zTree->SetBranchAddress ("event_weight",  &z_event_weight);
      zTree->SetBranchAddress ("fcal_et",       &z_fcal_et);
      //zTree->SetBranchAddress ("zdcEnergy",     &z_zdcEnergy);
      //zTree->SetBranchAddress ("q2",            &z_q2);
      //zTree->SetBranchAddress ("psi2",          &z_psi2);
      zTree->SetBranchAddress ("vz",            &z_vz);
      zTree->SetBranchAddress ("ntrk",          &z_ntrk);
      zTree->SetBranchAddress ("trk_pt",        &z_trk_pt);
      zTree->SetBranchAddress ("trk_eta",       &z_trk_eta);
      zTree->SetBranchAddress ("trk_phi",       &z_trk_phi);
    }
    zTree->SetBranchAddress ("isEE",          &_isEE);
    zTree->SetBranchAddress ("z_pt",          &_z_pt);
    zTree->SetBranchAddress ("z_y",           &_z_y);
    zTree->SetBranchAddress ("z_phi",         &_z_phi);
    zTree->SetBranchAddress ("z_m",           &_z_m);
    zTree->SetBranchAddress ("phi_transmin",  &phi_transmin);
    zTree->SetBranchAddress ("phi_transmax",  &phi_transmax);
    zTree->SetBranchAddress ("l1_pt",         &_l1_pt);
    zTree->SetBranchAddress ("l1_eta",        &_l1_eta);
    zTree->SetBranchAddress ("l1_phi",        &_l1_phi);
    zTree->SetBranchAddress ("l1_trk_pt",     &_l1_trk_pt);
    zTree->SetBranchAddress ("l1_trk_eta",    &_l1_trk_eta);
    zTree->SetBranchAddress ("l1_trk_phi",    &_l1_trk_phi);
    zTree->SetBranchAddress ("l1_charge",     &_l1_charge);
    zTree->SetBranchAddress ("l2_pt",         &_l2_pt);
    zTree->SetBranchAddress ("l2_eta",        &_l2_eta);
    zTree->SetBranchAddress ("l2_phi",        &_l2_phi);
    zTree->SetBranchAddress ("l2_trk_pt",     &_l2_trk_pt);
    zTree->SetBranchAddress ("l2_trk_eta",    &_l2_trk_eta);
    zTree->SetBranchAddress ("l2_trk_phi",    &_l2_trk_phi);
    zTree->SetBranchAddress ("l2_charge",     &_l2_charge);

    mixedEventsTree->Branch ("run_number",      &run_number,      "run_number/i");
    mixedEventsTree->Branch ("event_number",    &event_number,    "event_number/i");
    mixedEventsTree->Branch ("isEE",            &isEE,            "isEE/O");
    mixedEventsTree->Branch ("fcal_et",         &fcal_et,         "fcal_et/F");
    mixedEventsTree->Branch ("vz",              &vz,              "vz/F");
    mixedEventsTree->Branch ("ntrk",            &ntrk,            "ntrk/I");
    mixedEventsTree->Branch ("trk_pt",          &trk_pt,          "trk_pt[ntrk]/F");
    mixedEventsTree->Branch ("trk_eta",         &trk_eta,         "trk_eta[ntrk]/F");
    mixedEventsTree->Branch ("trk_phi",         &trk_phi,         "trk_phi[ntrk]/F");
    mixedEventsTree->Branch ("z_run_number",    &z_run_number,    "z_run_number/i");
    mixedEventsTree->Branch ("z_event_number",  &z_event_number,  "z_event_number/i");
    mixedEventsTree->Branch ("z_event_weight",  &z_event_weight,  "z_event_weight/F");
    mixedEventsTree->Branch ("z_vz",            &z_vz,            "z_vz/F");
    mixedEventsTree->Branch ("z_pt",            &z_pt,            "z_pt/F");
    mixedEventsTree->Branch ("z_y",             &z_y,             "z_y/F");
    mixedEventsTree->Branch ("z_phi",           &z_phi,           "z_phi/F");
    mixedEventsTree->Branch ("z_m",             &z_m,             "z_m/F");
    mixedEventsTree->Branch ("l1_pt",           &l1_pt,           "l1_pt/F");
    mixedEventsTree->Branch ("l1_eta",          &l1_eta,          "l1_eta/F");
    mixedEventsTree->Branch ("l1_phi",          &l1_phi,          "l1_phi/F");
    mixedEventsTree->Branch ("l1_trk_pt",       &l1_trk_pt,       "l1_trk_pt/F");
    mixedEventsTree->Branch ("l1_trk_eta",      &l1_trk_eta,      "l1_trk_eta/F");
    mixedEventsTree->Branch ("l1_trk_phi",      &l1_trk_phi,      "l1_trk_phi/F");
    mixedEventsTree->Branch ("l1_charge",       &l1_charge,       "l1_charge/I");
    mixedEventsTree->Branch ("l2_pt",           &l2_pt,           "l2_pt/F");
    mixedEventsTree->Branch ("l2_eta",          &l2_eta,          "l2_eta/F");
    mixedEventsTree->Branch ("l2_phi",          &l2_phi,          "l2_phi/F");
    mixedEventsTree->Branch ("l2_trk_pt",       &l2_trk_pt,       "l2_trk_pt/F");
    mixedEventsTree->Branch ("l2_trk_eta",      &l2_trk_eta,      "l2_trk_eta/F");
    mixedEventsTree->Branch ("l2_trk_phi",      &l2_trk_phi,      "l2_trk_phi/F");
    mixedEventsTree->Branch ("l2_charge",       &l2_charge,       "l2_charge/I");
    mixedEventsTree->Branch ("z_ntrk",          &z_ntrk,          "z_ntrk/I");
    mixedEventsTree->Branch ("z_trk_pt",        &z_trk_pt,        "z_trk_pt[z_ntrk]/F");
    mixedEventsTree->Branch ("z_trk_eta",       &z_trk_eta,       "z_trk_eta[z_ntrk]/F");
    mixedEventsTree->Branch ("z_trk_phi",       &z_trk_phi,       "z_trk_phi[z_ntrk]/F");


    for (int iZEvt = 0; iZEvt < mixingFraction*nZEvts; iZEvt++) {
      if (mixingFraction*nZEvts > 100 && iZEvt % (mixingFraction*nZEvts / 100) == 0)
        cout << iZEvt / (mixingFraction*nZEvts / 100) << "\% done...\r" << flush;

      if (doPPMixingVar && iZEvt >= nZEvts) // temporary max mixing fraction of 10 in central-most bin
        continue;

      zTree->GetEntry (zEventOrder[iZEvt % nZEvts]);
      zEventsUsed[iZEvt % nZEvts]++;

      isEE = _isEE;
      z_pt = _z_pt;
      z_y = _z_y;
      z_phi = _z_phi;
      z_m = _z_m;
      l1_pt = _l1_pt;
      l1_eta = _l1_eta;
      l1_phi = _l1_phi;
      l1_trk_pt = _l1_trk_pt;
      l1_trk_eta = _l1_trk_eta;
      l1_trk_phi = _l1_trk_phi;
      l1_charge = _l1_charge;
      l2_pt = _l2_pt;
      l2_eta = _l2_eta;
      l2_phi = _l2_phi;
      l2_trk_pt = _l2_trk_pt;
      l2_trk_eta = _l2_trk_eta;
      l2_trk_phi = _l2_trk_phi;
      l2_charge = _l2_charge;
      if (doSameFileMixing) {
        z_run_number = run_number;
        z_event_number = event_number;
        z_event_weight = event_weight;
        z_fcal_et = fcal_et;
        //z_zdcEnergy = zdcEnergy;
        //z_q2x_a = q2x_a;
        //z_q2x_c = q2x_c;
        //z_q2y_a = q2y_a;
        //z_q2y_c = q2y_c;
        z_vz = vz;
        z_ntrk = 0; // placeholder
      }

      if (fabs (z_vz) > 150)
        continue;

      if (z_event_weight == 0)
        continue;
      if (isnan (z_event_weight))
        cout << "Warning: z_event_weight = NAN" << endl;

      const short iPtZ = GetPtZBin (z_pt);
      if (iPtZ < 0 || iPtZ > nPtZBins-1)
        continue;

      if (iPtZ < 2)
        continue; // skip unneeded events

      // Find the next unused minimum bias event
      {
        //const short iFCalEt = GetPPCentBin (z_fcal_et);
        //if (iFCalEt < 0 || iFCalEt > numppCentBins)
        //  continue;

        // disable big branches during event mixing
        mbTree->SetBranchStatus ("trk_pt", 0);
        mbTree->SetBranchStatus ("trk_eta", 0);
        mbTree->SetBranchStatus ("trk_phi", 0);

        bool goodMixEvent = false;
        const int _iMBEvt = iMBEvt;
        do {
          iMBEvt = (iMBEvt+1) % nMBEvts;
          mbTree->GetEntry (mbEventOrder[iMBEvt]);
          goodMixEvent = (fabs (vz) < 150 && event_weight != 0); // always require these conditions

          if (!doSameFileMixing || mixingFraction == 1)
            goodMixEvent = (goodMixEvent && !mbEventsUsed[iMBEvt]); // checks for uniqueness (if applicable)

          if (doSameFileMixing)
            goodMixEvent = (goodMixEvent && iMBEvt != iZEvt); // don't mix with the exact same event

          if (doPPMixingVar) {
            //goodMixEvent = (goodMixEvent && 1 < _z_pt && _z_pt < 8 && DeltaPhi (phi_transmin, z_phi) < pi/6.);
            goodMixEvent = (goodMixEvent && 1 < _z_pt && _z_pt < 8 && DeltaPhi (phi_transmax, z_phi) < pi/6.);
            //goodMixEvent = (goodMixEvent && 1 < _z_pt && _z_pt < 8 && pi/4. < dphi && dphi < 3.*pi/4.); // variation on mixing: only mix with perpendicular, low-pT Z's
          }
        } while (!goodMixEvent && iMBEvt != _iMBEvt); // only check each event once
        if (_iMBEvt == iMBEvt) {
          cout << "No minbias event to mix with!!! Wrapped around on the same Z!!!" << endl;
          continue;
        }
        mbEventsUsed[iMBEvt] = true;

        // reenable big branches after finding a suitable event
        mbTree->SetBranchStatus ("trk_pt", 1);
        mbTree->SetBranchStatus ("trk_eta", 1);
        mbTree->SetBranchStatus ("trk_phi", 1);
        mbTree->GetEntry (mbEventOrder[iMBEvt]);
      }

      // at this point we have a Z boson and a new (unique & random) event to mix with
      mixedEventsTree->Fill ();

      event_weight *= z_event_weight;

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp

      h_pp_vz->Fill (vz);
      h_pp_vz_reweighted->Fill (vz, event_weight);

      h_pp_nch->Fill (ntrk);
      h_pp_nch_reweighted->Fill (ntrk, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (2.5, pow (event_weight, 2));

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];
        const float xhz = trkpt / z_pt;

        if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
          continue;

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], false);
        const double trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], false);
        if (trkEff == 0 || trkPur == 0)
          continue;
        const float trkWeight = trkPur / trkEff;

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        float dphi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, event_weight * trkWeight);
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (doPPMixingVar && idPhi == numPhiBins && 7.*pi/8. < dphi && dphi <= pi) {
            h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt);
            h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, event_weight * trkWeight);
            h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt / z_pt, event_weight * trkWeight);
          }
          else if (!doPPMixingVar && phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt);
            h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, event_weight * trkWeight);
            h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt / z_pt, event_weight * trkWeight);
          }
        }

        if ((doPPMixingVar && 7.*pi/8. < dphi && dphi <= pi) || (!doPPMixingVar && 3*pi/4 <= dphi)) {
          if (ptTrkBins[iPtZ][0] <= trkpt) {
            short iPt = 0;
            while (iPt < nPtTrkBins[iPtZ] && ptTrkBins[iPtZ][iPt+1] < trkpt) iPt++;
            if (iPt < 6) {
              trks_counts[0][iPt]   += 1;
              trks_weights1[0][iPt] += trkWeight;
              trks_weights2[0][iPt] += pow (trkWeight, 2);
            }
          }
          if (xHZBins[iPtZ][0] <= xhz) {
            short iX = 0;
            while (iX < nXHZBins[iPtZ] && xHZBins[iPtZ][iX+1] < xhz) iX++;
            if (iX < 6) {
              trks_counts[1][iX]   += 1;
              trks_weights1[1][iX] += trkWeight;
              trks_weights2[1][iX] += pow (trkWeight, 2);
            }
          }
        }
      } // end loop over tracks

      // fill yield histograms and covariance matrices
      for (int i1 = 0; i1 < nPtTrkBins[iPtZ]; i1++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinContent (i1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinContent (i1+1) + event_weight*trks_weights1[0][i1]);
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinError (i1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinError (i1+1) + event_weight*trks_weights1[0][i1]);
        //h_trk_pt_ptz_wgts[iSpc][iPtZ][iCent]->SetBinContent (i1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinContent (i1+1) + event_weight*sqrt(trks_weights2[0][i1]));
        for (int i2 = 0 ; i2 < nPtTrkBins[iPtZ]; i2++)
          h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (i1+1, i2+1, h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (i1+1, i2+1) + event_weight * (trks_weights1[0][i1]) * (trks_weights1[0][i2]));
          //h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (i1+1, i2+1, h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (i1+1, i2+1) + event_weight*(trks_weights1[0][i1])*(trks_weights1[0][i2]));
      } // end loop over i1
      for (int i1 = 0; i1 < nXHZBins[iPtZ]; i1++) {
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinContent (i1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinContent (i1+1) + event_weight*trks_weights1[1][i1]);
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinError (i1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinError (i1+1) + event_weight*trks_weights1[1][i1]);
        //h_trk_xhz_ptz_wgts[iSpc][iPtZ][iCent]->SetBinContent (i1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinContent (i1+1) + event_weight*sqrt (trks_weights2[1][i1]));
        for (int i2 = 0 ; i2 < nXHZBins[iPtZ]; i2++)
          h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (i1+1, i2+1, h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (i1+1, i2+1) + event_weight * (trks_weights1[1][i1]) * (trks_weights1[1][i2]));
          //h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (i1+1, i2+1, h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (i1+1, i2+1) + event_weight*(trks_weights1[1][i1])*(trks_weights1[1][i2]));
      } // end loop over i1

      // reset trk count measurements for next event
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 6; j++) {
          trks_counts[i][j] = 0;
          trks_weights1[i][j] = 0.;
          trks_weights2[i][j] = 0.;
        } // end loop over j
      } // end loop over i

    } // end loop over pp tree

    cout << endl << "Done minbias pp loop." << endl;
  }

  for (int i = 0; i < nZEvts; i++) {
    if (zEventsUsed[i] < mixingFraction) {
      cout << "Warning! Some events used fewer than " << mixingFraction << " times!" << endl;
      break;
    }
  }

  mixedEventsTree->SetDirectory (mixedEventsFile);
  mixedEventsTree->Write ("", TObject :: kOverwrite);

  mixedEventsFile->Close ();
  SaferDelete (mixedEventsFile);

  Delete2DArray (trks_counts, 2, 6);
  Delete2DArray (trks_weights1, 2, 6);
  Delete2DArray (trks_weights2, 2, 6);

  SaveHists (outFileName);

  if (inFile) inFile->Close ();
  SaferDelete (inFile);

  if (mbInFile) mbInFile->Close ();
  SaferDelete (mbInFile);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Generates weights from pre-mixed events.
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: GenerateWeights (const char* inFileName, const char* outFileName) {

  const int nQ2Bins = 20;
  const double* q2Bins = linspace (0, 0.3, nQ2Bins);
  const int nPsi2Bins = 8;
  const double* psi2Bins = linspace (0, pi/2, nPsi2Bins);
  
  //const int nNchBins = 160;
  //const double* nchBins = linspace (-0.5, 160.5, nNchBins);

  TH1D* h_z_q2_dist[3][numFineCentBins][nPtZBins];
  TH1D* h_mixed_q2_dist[3][numFineCentBins][nPtZBins];
  TH1D* h_z_psi2_dist[3][numFineCentBins][nPtZBins];
  TH1D* h_mixed_psi2_dist[3][numFineCentBins][nPtZBins];

  TFile* inFile = new TFile (inFileName, "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbMixedTree");

  for (int iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (int iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
        h_z_q2_dist[iSpc][iFineCent][iPtZ] = new TH1D (Form ("h_z_q2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()), "", nQ2Bins, q2Bins);
        h_mixed_q2_dist[iSpc][iFineCent][iPtZ] = new TH1D (Form ("h_mixed_q2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()), "", nQ2Bins, q2Bins);
        h_z_q2_dist[iSpc][iFineCent][iPtZ]->Sumw2 ();
        h_mixed_q2_dist[iSpc][iFineCent][iPtZ]->Sumw2 ();

        h_z_psi2_dist[iSpc][iFineCent][iPtZ] = new TH1D (Form ("h_z_psi2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()), "", nPsi2Bins, psi2Bins);
        h_mixed_psi2_dist[iSpc][iFineCent][iPtZ] = new TH1D (Form ("h_mixed_psi2_dist_%s_iCent%i_iPtZ%i_%s", spc, iFineCent, iPtZ, name.c_str ()), "", nPsi2Bins, psi2Bins);
        h_z_psi2_dist[iSpc][iFineCent][iPtZ]->Sumw2 ();
        h_mixed_psi2_dist[iSpc][iFineCent][iPtZ]->Sumw2 ();
      }
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    bool isEE = false;
    float fcal_et = 0, z_pt = 0, z_phi = 0, q2 = 0, psi2 = 0, z_q2 = 0, z_psi2 = 0;

    PbPbTree->SetBranchAddress ("isEE",         &isEE);
    PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    PbPbTree->SetBranchAddress ("z_pt",         &z_pt);
    PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
    PbPbTree->SetBranchAddress ("q2",           &q2);
    PbPbTree->SetBranchAddress ("psi2",         &psi2);
    PbPbTree->SetBranchAddress ("z_q2",         &z_q2);
    PbPbTree->SetBranchAddress ("z_psi2",       &z_psi2);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      const short iFineCent = GetFineCentBin (fcal_et);
      if (iFineCent < 1 || iFineCent > numFineCentBins-1)
        continue;

      const short iPtZ = GetPtZBin (z_pt);
      if (iPtZ < 0 || iPtZ > nPtZBins-1)
        continue;

      const short iSpc = (isEE ? 0 : 1);

      h_z_q2_dist[iSpc][iFineCent][iPtZ]->Fill (z_q2);
      h_mixed_q2_dist[iSpc][iFineCent][iPtZ]->Fill (q2);

      float dphi = DeltaPhi (z_phi, z_psi2, false);
      if (dphi > pi/2)
        dphi = pi - dphi;
      h_z_psi2_dist[iSpc][iFineCent][iPtZ]->Fill (dphi);

      dphi = DeltaPhi (z_phi, psi2, false);
      if (dphi > pi/2)
        dphi = pi - dphi;
      h_mixed_psi2_dist[iSpc][iFineCent][iPtZ]->Fill (dphi);
    }
    cout << "Done Pb+Pb loop." << endl;
  }

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (int ix = 1; ix <= h_z_q2_dist[iSpc][iFineCent][iPtZ]->GetNbinsX (); ix++)
          h_z_q2_dist[iSpc][iFineCent][iPtZ]->SetBinError (ix, h_z_q2_dist[iSpc][iFineCent][iPtZ]->GetBinError (ix) * TMath::Sqrt (40));
        for (int ix = 1; ix <= h_z_psi2_dist[iSpc][iFineCent][iPtZ]->GetNbinsX (); ix++)
          h_z_psi2_dist[iSpc][iFineCent][iPtZ]->SetBinError (ix, h_z_psi2_dist[iSpc][iFineCent][iPtZ]->GetBinError (ix) * TMath::Sqrt (40));

        h_z_q2_dist[2][iFineCent][iPtZ]->Add (h_z_q2_dist[iSpc][iFineCent][iPtZ]);
        h_mixed_q2_dist[2][iFineCent][iPtZ]->Add (h_mixed_q2_dist[iSpc][iFineCent][iPtZ]);
        h_z_psi2_dist[2][iFineCent][iPtZ]->Add (h_z_psi2_dist[iSpc][iFineCent][iPtZ]);
        h_mixed_psi2_dist[2][iFineCent][iPtZ]->Add (h_mixed_psi2_dist[iSpc][iFineCent][iPtZ]);
      }
    }
  }

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_q2_dist[iSpc][iFineCent][iPtZ]->Scale (1./h_z_q2_dist[iSpc][iFineCent][iPtZ]->Integral ());
        h_mixed_q2_dist[iSpc][iFineCent][iPtZ]->Scale (1./h_mixed_q2_dist[iSpc][iFineCent][iPtZ]->Integral ());
        h_z_psi2_dist[iSpc][iFineCent][iPtZ]->Scale (1./h_z_psi2_dist[iSpc][iFineCent][iPtZ]->Integral ());
        h_mixed_psi2_dist[iSpc][iFineCent][iPtZ]->Scale (1./h_mixed_psi2_dist[iSpc][iFineCent][iPtZ]->Integral ());
      }
    }
  }

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_q2_dist[iSpc][iFineCent][iPtZ]->Divide (h_mixed_q2_dist[iSpc][iFineCent][iPtZ]);
        h_mixed_q2_dist[iSpc][iFineCent][iPtZ]->Divide (h_mixed_q2_dist[iSpc][iFineCent][iPtZ]);
        h_z_psi2_dist[iSpc][iFineCent][iPtZ]->Divide (h_mixed_psi2_dist[iSpc][iFineCent][iPtZ]);
        h_mixed_psi2_dist[iSpc][iFineCent][iPtZ]->Divide (h_mixed_psi2_dist[iSpc][iFineCent][iPtZ]);
      }
    }
  }

  if (eventWeightsFile && eventWeightsFile->IsOpen ()) {
    eventWeightsFile->Close ();
    eventWeightsLoaded = false;
  }

  eventWeightsFile = new TFile (outFileName, "recreate");
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iFineCent = 0; iFineCent < numFineCentBins; iFineCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        SafeWrite (h_z_q2_dist[iSpc][iFineCent][iPtZ]);
        SafeWrite (h_mixed_q2_dist[iSpc][iFineCent][iPtZ]);
        SafeWrite (h_z_psi2_dist[iSpc][iFineCent][iPtZ]);
        SafeWrite (h_mixed_psi2_dist[iSpc][iFineCent][iPtZ]);
      }
    }
  }
  eventWeightsFile->Close ();
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot centrality distributions for different MB triggers
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: PlotCentralityDists () {
  const char* canvasName = "c_centrality";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    gDirectory->Add (c);
    c->cd ();
  }

  c->cd ();

  c->SetLogy ();

  for (short iMBTrig = 0; iMBTrig < 3; iMBTrig++) {
    h_centrality[iMBTrig]->Scale (1., "width");

    h_centrality[iMBTrig]->SetLineColor (colors[iMBTrig+1]);

    h_centrality[iMBTrig]->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [GeV]");
    h_centrality[iMBTrig]->GetYaxis ()->SetTitle ("dN_{evt} / d#Sigma#it{E}_{T} [GeV^{-1}]");

//    h_centrality[iMBTrig]->GetYaxis ()->SetRangeUser (5e-2, 2e3);

    h_centrality[iMBTrig]->Draw (!canvasExists && iMBTrig == 0 ? "hist" : "same hist");
  }

  myMarkerTextNoLine (0.67, 0.81, colors[1], kFullCircle, "HLT_mb_sptrk_L1ZDC_A_C_VTE50 (PC)", 1.25, 0.04);
  myMarkerTextNoLine (0.67, 0.76, colors[2], kFullCircle, "HLT_noalg_pc_L1TE50_VTE600_0ETA49 (PC)", 1.25, 0.04);
  myMarkerTextNoLine (0.67, 0.71, colors[3], kFullCircle, "HLT_noalg_cc_L1TE600_0ETA49 (CC)", 1.25, 0.04);

  myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
  myText (0.22, 0.81, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.04);

  c->SaveAs (Form ("%s/CentralityDist.pdf", plotPath.Data ()));
}

#endif
