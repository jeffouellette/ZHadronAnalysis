#ifndef __MixedMCAnalysis_h__
#define __MixedMCAnalysis_h__

#include "Params.h"
#include "FullAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class MixedMCAnalysis : public FullAnalysis {

  public:

  bool takeNonTruthTracks = false;
  bool doQ2Mixing = false;
  bool doPsi2Mixing = false;

  int numQ2MixBins = 1;
  double* q2MixBins = nullptr;
  int numPsi2MixBins = 1;
  double* psi2MixBins = nullptr;

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
  

  MixedMCAnalysis (const char* _name = "mc_mixed") : FullAnalysis () {
    name = _name;
    plotFill = true;
    useAltMarker = false;
    isMC = true;
    eventWeightsFileName = "MCAnalysis/Nominal/eventWeightsFile.root";
    plotSignal = false;
    hasBkg = false;
    backgroundSubtracted = true;
    histsUnfolded = true;
    iaaCalculated = true;
    //icpCalculated = true;
  }

  void Execute (const char* inFileName, const char* outFileName) override;
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
////////////////////////////////////////////////////////////////////////////////////////////////
void MixedMCAnalysis :: Execute (const char* inFileName, const char* outFileName) {

  LoadEventWeights ();

  SetupDirectories ("", "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
  if (!inFile || !inFile->IsOpen ()) {
    cout << "Error in Execute: File not open!" << endl;
    return;
  }
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), inFileName) << endl;

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  CreateHists ();

  bool isEE = false;
  int event_number = 0;
  float event_weight = 1;//, fcal_weight = 1;//, q2_weight = 1, psi2_weight = 1;//, vz_weight = 1, nch_weight = 1;
  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0;
  float qx_a = 0, qy_a = 0, qx_c = 0, qy_c = 0;
  float z_pt = 0, z_y = 0, z_phi = 0, z_m = 0;
  float l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0;
  float l1_trk_pt = 0, l1_trk_eta = 0, l1_trk_phi = 0, l2_trk_pt = 0, l2_trk_eta = 0, l2_trk_phi = 0;
  int l1_charge = 0, l2_charge = 0, ntrk = 0;
  float trk_pt[10000], trk_eta[10000], trk_phi[10000];
  bool trk_truth_matched[10000];

  q2MixBins = linspace (0, 0.2, numQ2MixBins);
  psi2MixBins = linspace (-pi/2, pi/2, numPsi2MixBins);


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    int nEvts = PbPbTree->GetEntries ();
    PbPbTree->LoadBaskets (3000000000);
    PbPbTree->SetBranchAddress ("event_number", &event_number);
    PbPbTree->SetBranchAddress ("event_weight", &event_weight);
    PbPbTree->SetBranchAddress ("isEE",         &isEE);
    PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    PbPbTree->SetBranchAddress ("qx_a",         &qx_a);
    PbPbTree->SetBranchAddress ("qy_a",         &qy_a);
    PbPbTree->SetBranchAddress ("qx_c",         &qx_c);
    PbPbTree->SetBranchAddress ("qy_c",         &qy_c);
    PbPbTree->SetBranchAddress ("q2",           &q2);
    PbPbTree->SetBranchAddress ("psi2",         &psi2);
    PbPbTree->SetBranchAddress ("vz",           &vz);
    PbPbTree->SetBranchAddress ("z_pt",         &z_pt);
    PbPbTree->SetBranchAddress ("z_y",          &z_y);
    PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
    PbPbTree->SetBranchAddress ("z_m",          &z_m);
    PbPbTree->SetBranchAddress ("l1_pt",        &l1_pt);
    PbPbTree->SetBranchAddress ("l1_eta",       &l1_eta);
    PbPbTree->SetBranchAddress ("l1_phi",       &l1_phi);
    PbPbTree->SetBranchAddress ("l1_charge",    &l1_charge);
    PbPbTree->SetBranchAddress ("l1_trk_pt",    &l1_trk_pt);
    PbPbTree->SetBranchAddress ("l1_trk_eta",   &l1_trk_eta);
    PbPbTree->SetBranchAddress ("l1_trk_phi",   &l1_trk_phi);
    PbPbTree->SetBranchAddress ("l2_pt",        &l2_pt);
    PbPbTree->SetBranchAddress ("l2_eta",       &l2_eta);
    PbPbTree->SetBranchAddress ("l2_phi",       &l2_phi);
    PbPbTree->SetBranchAddress ("l2_charge",    &l2_charge);
    PbPbTree->SetBranchAddress ("l2_trk_pt",    &l2_trk_pt);
    PbPbTree->SetBranchAddress ("l2_trk_eta",   &l2_trk_eta);
    PbPbTree->SetBranchAddress ("l2_trk_phi",   &l2_trk_phi);
    PbPbTree->SetBranchAddress ("ntrk",         &ntrk);
    PbPbTree->SetBranchAddress ("trk_pt",       &trk_pt);
    PbPbTree->SetBranchAddress ("trk_eta",      &trk_eta);
    PbPbTree->SetBranchAddress ("trk_phi",      &trk_phi);
    PbPbTree->SetBranchAddress ("trk_truth_matched", &trk_truth_matched);

    //cout << "Found " << nEvts << " events, assuming " << nEvts << " / " << mixingFraction << " = " << nEvts / mixingFraction << " are unique" << endl;
    //nEvts = nEvts / mixingFraction;

    //PbPbTree->GetEntry (0);
    //const int init_event_number = event_number;

    //if (nEvts == 0)
    //  cout << "Warning! No Z's to mix with in this run!" << endl;
    //cout << "For this PbPb tree, maximum mixing fraction = " << nMBPbPbEvts / nEvts << endl;

    //bool doShuffle = false;
    //if (mixingFraction * nEvts > nMBPbPbEvts) {
    //  cout << "Warning! Mixing fraction too high, will use " << (float)(nMBPbPbEvts / mixingFraction) / (float)(nEvts) * 100. << "% of Z events" << endl;
    //  nEvts = nMBPbPbEvts / mixingFraction;
    //  doShuffle = true;
    //}

    //std::vector <int> zEventOrder = {};
    //std::vector <int> zEventsUsed = {};
    //for (int i = 0; i < nEvts; i++) {
    //  zEventOrder.push_back (i);
    //  zEventsUsed.push_back (false);
    //}
    ////std::srand (std::time (0));
    //if (doShuffle) std::random_shuffle (zEventOrder.begin (), zEventOrder.end ());

    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;

      PbPbTree->GetEntry (iEvt);
      //PbPbTree->GetEntry (zEventOrder[iEvt % nEvts]);
      //zEventsUsed[iEvt % nEvts]++;

      //if (iEvt > 0 && event_number == init_event_number)
      //  break;

      //if (fabs (vz) > 150)
      //  continue; // vertex cut

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined

      const short iCent = GetCentBin (fcal_et);
      if (iCent < 1 || iCent > numCentBins-1)
        continue;

      const short iFineCent = GetFineCentBin (fcal_et);
      if (iFineCent < 1 || iFineCent > numFineCentBins-1)
        continue;

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin
      if (iPtZ < 0 || iPtZ > nPtZBins-1)
        continue;
      if (iPtZ < 2)
        continue; // no one cares about these events anyways

      //{
      //  float dphi = DeltaPhi (z_phi, psi2, false);
      //  if (dphi > pi/2)
      //    dphi = pi - dphi;
      //  fcal_weight = h_PbPbFCal_weights[iSpc][iPtZ]->GetBinContent (h_PbPbFCal_weights[iSpc][iPtZ]->FindBin (fcal_et));
      //  //q2_weight = h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->FindBin (q2));
      //  //psi2_weight = h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->FindBin (dphi));

      //  event_weight *= fcal_weight;// * psi2_weight;
      //}

      if (event_weight == 0)
        continue;

      // save info about this particular Z boson
      const float _event_weight = event_weight;
      const bool  _isEE         = isEE;
      const float _z_pt         = z_pt;
      const float _z_y          = z_y;
      const float _z_m          = z_m;
      const float _z_phi        = z_phi;
      //const float _l1_pt        = l1_pt;
      //const float _l1_eta       = l1_eta;
      //const float _l1_phi       = l1_phi;
      //const float _l2_pt        = l2_pt;
      //const float _l2_eta       = l2_eta;
      //const float _l2_phi       = l2_phi;
      //const float _l1_trk_pt    = l1_trk_pt;
      //const float _l1_trk_eta   = l1_trk_eta;
      //const float _l1_trk_phi   = l1_trk_phi;
      //const float _l2_trk_pt    = l2_trk_pt;
      //const float _l2_trk_eta   = l2_trk_eta;
      //const float _l2_trk_phi   = l2_trk_phi;

      const short iFCalEt = GetSuperFineCentBin (fcal_et);
      if (iFCalEt < 1 || iFCalEt > numSuperFineCentBins-1)
        continue;
      const short iQ2 = GetQ2MixBin (q2);
      if (doQ2Mixing && (iQ2 < 0 || iQ2 > numQ2MixBins-1)) {
        cout << "Out-of-bounds q2, skipping this Z!" << endl;
        continue;
      }
      const short iPsi2 = GetPsi2MixBin (psi2);
      if (doPsi2Mixing && (iPsi2 < 0 || iPsi2 > numPsi2MixBins-1)) {
        cout << "Out-of-bounds psi2, skipping this Z!" << endl;
        continue;
      }

      for (int iPrime = 1; iPrime <= mixingFraction; iPrime++) {

        if ((iPrime + iEvt) % nEvts == iEvt) {
          cout << "No more events to mix with!" << endl;
          break;
        }

        //cout << "Getting entry " << (iPrime+iEvt) %nEvts << endl;
        PbPbTree->GetEntry ((iPrime+iEvt) % nEvts);
        
        event_weight  = _event_weight;
        isEE          = _isEE;
        z_pt          = _z_pt;
        z_y           = _z_y;
        z_m           = _z_m;
        z_phi         = _z_phi;
        //l1_pt         = _l1_pt;
        //l1_eta        = _l1_eta;
        //l1_phi        = _l1_phi;
        //l2_pt         = _l2_pt;
        //l2_eta        = _l2_eta;
        //l2_phi        = _l2_phi;
        //l1_trk_pt     = _l1_trk_pt;
        //l1_trk_eta    = _l1_trk_eta;
        //l1_trk_phi    = _l1_trk_phi;
        //l2_trk_pt     = _l2_trk_pt;
        //l2_trk_eta    = _l2_trk_eta;
        //l2_trk_phi    = _l2_trk_phi;

        bool goodMixEvent = false;
        {
          const float qx = qx_a + qx_c;
          const float qy = qy_a + qy_c;
          q2 = sqrt (qx*qx + qy*qy) / fcal_et;
          psi2 = 0.5 * atan2 (qy, qx);
        }
        goodMixEvent = (fabs (vz) < 150 && event_weight != 0);
        goodMixEvent = (goodMixEvent && iFCalEt == GetSuperFineCentBin (fcal_et));
        if (doQ2Mixing)
          goodMixEvent = (goodMixEvent && iQ2 == GetQ2MixBin (q2));
        if (doPsi2Mixing)
          goodMixEvent = (goodMixEvent && iPsi2 == GetPsi2MixBin (psi2));

        if (!goodMixEvent) {
          iPrime--;
          continue;
        }

        h_fcal_et->Fill (fcal_et);
        h_fcal_et_reweighted->Fill (fcal_et, event_weight);

        h_q2[iFineCent]->Fill (q2);
        h_q2_reweighted[iFineCent]->Fill (q2, event_weight);
        h_psi2[iFineCent]->Fill (psi2);
        h_psi2_reweighted[iFineCent]->Fill (psi2, event_weight);
        h_PbPb_vz->Fill (vz);
        h_PbPb_vz_reweighted->Fill (vz, event_weight);

        h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
        h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);
        for (int iTrk = 0; iTrk < ntrk; iTrk++) {
          const float trkpt = trk_pt[iTrk];

          if (takeNonTruthTracks && trk_truth_matched[iTrk])
            continue;

          if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
            continue;

          if (doLeptonRejVar && (DeltaR (l1_trk_eta, trk_eta[iTrk], l1_trk_phi, trk_phi[iTrk]) < 0.03 || DeltaR (l2_trk_eta, trk_eta[iTrk], l2_trk_phi, trk_phi[iTrk]) < 0.03))
            continue;

          const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], true);
          const float trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], true);
          if (trkEff == 0 || trkPur == 0)
            continue;
          const float trkWeight = event_weight * trkPur / trkEff;

          // Study correlations (requires dphi in -pi/2 to 3pi/2)
          float dphi = DeltaPhi (z_phi, trk_phi[iTrk], true);
          if (dphi < -pi/2)
            dphi = dphi + 2*pi;

          for (short iPtTrk = 0; iPtTrk < nPtchBins[iPtZ]; iPtTrk++) {
            if (pTchBins[iPtZ][iPtTrk] <= trkpt && trkpt < pTchBins[iPtZ][iPtTrk+1])
              h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, trkWeight);
          }

          // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
          dphi = DeltaPhi (z_phi, trk_phi[iTrk], false);
          for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
            if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
              h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt);
              h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, trkWeight);
              h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt / z_pt, trkWeight);
            }
          }
        } // end loop over tracks

      } // end loop over subsequent events (for mixing)
    } // end loop over Pb+Pb tree
    cout << "Done MC Pb+Pb loop." << endl;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over pp tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (ppTree) {
    int nEvts = ppTree->GetEntries ();
    ppTree->LoadBaskets (4000000000);
    ppTree->SetBranchAddress ("event_number", &event_number);
    ppTree->SetBranchAddress ("event_weight", &event_weight);
    ppTree->SetBranchAddress ("isEE",       &isEE);
    ppTree->SetBranchAddress ("vz",         &vz);
    ppTree->SetBranchAddress ("z_pt",       &z_pt);
    ppTree->SetBranchAddress ("z_y",        &z_y);
    ppTree->SetBranchAddress ("z_phi",      &z_phi);
    ppTree->SetBranchAddress ("z_m",        &z_m);
    ppTree->SetBranchAddress ("l1_pt",      &l1_pt);
    ppTree->SetBranchAddress ("l1_eta",     &l1_eta);
    ppTree->SetBranchAddress ("l1_phi",     &l1_phi);
    ppTree->SetBranchAddress ("l1_charge",  &l1_charge);
    ppTree->SetBranchAddress ("l1_trk_pt",  &l1_trk_pt);
    ppTree->SetBranchAddress ("l1_trk_eta", &l1_trk_eta);
    ppTree->SetBranchAddress ("l1_trk_phi", &l1_trk_phi);
    ppTree->SetBranchAddress ("l2_pt",      &l2_pt);
    ppTree->SetBranchAddress ("l2_eta",     &l2_eta);
    ppTree->SetBranchAddress ("l2_phi",     &l2_phi);
    ppTree->SetBranchAddress ("l2_charge",  &l2_charge);
    ppTree->SetBranchAddress ("l2_trk_pt",  &l2_trk_pt);
    ppTree->SetBranchAddress ("l2_trk_eta", &l2_trk_eta);
    ppTree->SetBranchAddress ("l2_trk_phi", &l2_trk_phi);
    ppTree->SetBranchAddress ("ntrk",       &ntrk);
    ppTree->SetBranchAddress ("trk_pt",     &trk_pt);
    ppTree->SetBranchAddress ("trk_eta",    &trk_eta);
    ppTree->SetBranchAddress ("trk_phi",    &trk_phi);
    ppTree->SetBranchAddress ("trk_truth_matched", &trk_truth_matched);

    //cout << "Found " << nEvts << " events, assuming " << nEvts << " / " << mixingFraction << " = " << nEvts / mixingFraction << " are unique" << endl;
    //nEvts = nEvts / mixingFraction;

    //ppTree->GetEntry (0);
    //const int init_event_number = event_number;

    //if (nEvts == 0)
    //  cout << "Warning! No Z's to mix with in this run!" << endl;
    //cout << "For this pp tree, maximum mixing fraction = " << nMBppEvts / nEvts << endl;

    //bool doShuffle = false;
    //if (mixingFraction * nEvts > nMBppEvts) {
    //  cout << "Warning! Mixing fraction too high, will use " << (float)(nMBppEvts / mixingFraction) / (float)(nEvts) * 100. << "% of Z events" << endl;
    //  nEvts = nMBppEvts / mixingFraction;
    //  doShuffle = true;
    //}

    //std::vector <int> zEventOrder = {};
    //std::vector <int> zEventsUsed = {};
    //for (int i = 0; i < nEvts; i++) {
    //  zEventOrder.push_back (i);
    //  zEventsUsed.push_back (false);
    //}
    ////std::srand (std::time (0));
    //if (doShuffle) std::random_shuffle (zEventOrder.begin (), zEventOrder.end ());

    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;

      ppTree->GetEntry (iEvt);
      //ppTree->GetEntry (zEventOrder[iEvt % nEvts]);
      //zEventsUsed[iEvt % nEvts]++;

      //if (iEvt > 0 && event_number == init_event_number)
      //  break;

      //if (fabs (vz) > 150)
      //  continue; // vertex cut

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin
      if (iPtZ < 0 || iPtZ > nPtZBins-1)
        continue;
      if (iPtZ < 2)
        continue; // no one cares about these events anyways

      //nch_weight = h_ppNch_weights->GetBinContent (h_ppNch_weights->FindBin (ntrk));

      //event_weight = event_weight * vz_weight * nch_weight;
      if (event_weight == 0)
        continue;

      // save info about this particular Z boson
      const float _event_weight = event_weight;
      const float _z_pt         = z_pt;
      const float _z_y          = z_y;
      const float _z_m          = z_m;
      const float _z_phi        = z_phi;
      //const float _l1_pt        = l1_pt;
      //const float _l1_eta       = l1_eta;
      //const float _l1_phi       = l1_phi;
      //const float _l2_pt        = l2_pt;
      //const float _l2_eta       = l2_eta;
      //const float _l2_phi       = l2_phi;
      //const float _l1_trk_pt    = l1_trk_pt;
      //const float _l1_trk_eta   = l1_trk_eta;
      //const float _l1_trk_phi   = l1_trk_phi;
      //const float _l2_trk_pt    = l2_trk_pt;
      //const float _l2_trk_eta   = l2_trk_eta;
      //const float _l2_trk_phi   = l2_trk_phi;

      for (int iPrime = 1; iPrime <= mixingFraction; iPrime++) {
        ppTree->GetEntry ((iPrime+iEvt) % nEvts);
        
        event_weight  = _event_weight;
        z_pt          = _z_pt;
        z_y           = _z_y;
        z_m           = _z_m;
        z_phi         = _z_phi;
        //l1_pt         = _l1_pt;
        //l1_eta        = _l1_eta;
        //l1_phi        = _l1_phi;
        //l2_pt         = _l2_pt;
        //l2_eta        = _l2_eta;
        //l2_phi        = _l2_phi;
        //l1_trk_pt     = _l1_trk_pt;
        //l1_trk_eta    = _l1_trk_eta;
        //l1_trk_phi    = _l1_trk_phi;
        //l2_trk_pt     = _l2_trk_pt;
        //l2_trk_eta    = _l2_trk_eta;
        //l2_trk_phi    = _l2_trk_phi;

        h_pp_vz->Fill (vz);
        h_pp_vz_reweighted->Fill (vz, event_weight);

        h_pp_nch->Fill (ntrk);
        h_pp_nch_reweighted->Fill (ntrk, event_weight);

        h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
        h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);
        for (int iTrk = 0; iTrk < ntrk; iTrk++) {
          const float trkpt = trk_pt[iTrk];

          if (takeNonTruthTracks && trk_truth_matched[iTrk])
            continue;

          if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
            continue;

          if (doLeptonRejVar && (DeltaR (l1_trk_eta, trk_eta[iTrk], l1_trk_phi, trk_phi[iTrk]) < 0.03 || DeltaR (l2_trk_eta, trk_eta[iTrk], l2_trk_phi, trk_phi[iTrk]) < 0.03))
            continue;

          const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], false);
          const float trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], false);
          if (trkEff == 0 || trkPur == 0)
            continue;
          const float trkWeight = event_weight * trkPur / trkEff;

          // Study correlations (requires dphi in -pi/2 to 3pi/2)
          float dphi = DeltaPhi (z_phi, trk_phi[iTrk], true);
          if (dphi < -pi/2)
            dphi = dphi + 2*pi;

          for (short iPtTrk = 0; iPtTrk < nPtchBins[iPtZ]; iPtTrk++) {
            if (pTchBins[iPtZ][iPtTrk] <= trkpt && trkpt < pTchBins[iPtZ][iPtTrk+1])
              h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, trkWeight);
          }

          // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
          dphi = DeltaPhi (z_phi, trk_phi[iTrk], false);
          for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
            if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
              h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt);
              h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, trkWeight);
              h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt / z_pt, trkWeight);
            }
          }
        } // end loop over tracks

      } // end loop over subsequent events (for mixing)
    } // end loop over pp tree
    cout << "Done MC closure pp loop." << endl;
  }

  //CombineHists ();
  //ScaleHists ();

  SaveHists (outFileName);

  inFile->Close ();
  if (inFile) { delete inFile; inFile = nullptr; }

  //Delete2DArray (trkPtProj, numPhiBins, nPtchBins[iPtZ]);
}


#endif
