#ifndef __Cuts_cxx__
#define __Cuts_cxx__

#include "Cuts.h"
#include "Params.h"
#include "TreeVariables.h"

#include <iostream>

namespace ZTrackAnalyzer {

bool IsElectronTrack (TreeVariables* t, const int iTrk, const int iE1, const int iE2) {

  //const float minDR = 0.03;

  //if (iE1 * iE2 < 0) {
  //  cout << "Warning in Cuts.cxx: invalid candidate electron indices! Returning false by default." << endl;
  //  return false;
  //}

  bool isElectron = false;
  //double dR = DeltaR (t->trk_eta[iTrk], t->electron_id_track_eta[iE1], t->trk_phi[iTrk], t->electron_id_track_phi[iE1]);
  isElectron |= (iE1 != -1 && t->trk_pt[iTrk] == t->electron_id_track_pt[iE1] && t->trk_eta[iTrk] == t->electron_id_track_eta[iE1] && t->trk_phi[iTrk] == t->electron_id_track_phi[iE1]);
  //isElectron = isElectron || dR < minDR;
  //dR = DeltaR (t->trk_eta[iTrk], t->electron_id_track_eta[iE2], t->trk_phi[iTrk], t->electron_id_track_phi[iE2]);
  isElectron |= (iE2 != -1 && t->trk_pt[iTrk] == t->electron_id_track_pt[iE2] && t->trk_eta[iTrk] == t->electron_id_track_eta[iE2] && t->trk_phi[iTrk] == t->electron_id_track_phi[iE2]);
  //isElectron = isElectron || dR < minDR;

  return isElectron;

}


bool IsMuonTrack (TreeVariables* t, const int iTrk, const int iM1, const int iM2) {

  //const float minDR = 0.03;

  //if (iM1 * iM2 < 0) {
  //  cout << "Warning in Cuts.cxx: invalid candidate muon indices! Returning false by default." << endl;
  //  return false;
  //}

  bool isMuon = false;
  //double dR = DeltaR (t->trk_eta[iTrk], t->muon_id_track_eta[iM1], t->trk_phi[iTrk], t->muon_id_track_phi[iM1]);
  isMuon |= (iM1 != -1 && t->trk_pt[iTrk] == t->muon_id_track_pt[iM1] && t->trk_eta[iTrk] == t->muon_id_track_eta[iM1] && t->trk_phi[iTrk] == t->muon_id_track_phi[iM1]);
  //isMuon = isMuon || dR < minDR;
  //dR = DeltaR (t->trk_eta[iTrk], t->muon_id_track_eta[iM2], t->trk_phi[iTrk], t->muon_id_track_phi[iM2]);
  isMuon |= (iM2 != -1 && t->trk_pt[iTrk] == t->muon_id_track_pt[iM2] && t->trk_eta[iTrk] == t->muon_id_track_eta[iM2] && t->trk_phi[iTrk] == t->muon_id_track_phi[iM2]);
  //isMuon = isMuon || dR < minDR;

  return isMuon;
}

} // end namespace

#endif
