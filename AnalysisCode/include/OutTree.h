#ifndef __OutTree_h__
#define __OutTree_h__

#include <TTree.h>
#include <vector>

using namespace std;

namespace ZTrackAnalyzer {

bool isMC = false;
float event_weight = 0;

unsigned int event_number = 0;
unsigned int run_number = 0;

float ip = 0;
float eventPlane = 0;

float fcal_et = 0;
float zdcEnergy = 0;
float q2x_a = 0;
float q2y_a = 0;
float q2x_c = 0;
float q2y_c = 0;
float q2 = 0;
float psi2 = 0;
float q3 = 0;
float psi3 = 0;
float q4 = 0;
float psi4 = 0;
float vz = 0;

bool isEE = false;
bool isSameSign = false;
bool isQCDBkg = false;
float z_pt = 0;
float z_y = 0;
float z_phi = 0;
float z_m = 0;

float phi_transmin = 0;
float phi_transmax = 0;

float truth_z_pt = 0;
float truth_z_y = 0;
float truth_z_phi = 0;
float truth_z_m = 0;

float l1_pt = 0;
float l1_eta = 0;
float l1_phi = 0;
int   l1_charge = 0;
float l1_trk_pt = 0;
float l1_trk_eta = 0;
float l1_trk_phi = 0;
float l1_pt_sys = 0;
float l1_eta_sys = 0;
float l1_phi_sys = 0;
float l1_d0sig = 0;
float l1_z0 = 0;

float l2_pt = 0;
float l2_eta = 0;
float l2_phi = 0;
float l2_trk_pt = 0;
float l2_trk_eta = 0;
float l2_trk_phi = 0;
int   l2_charge = 0;
float l2_pt_sys = 0;
float l2_eta_sys = 0;
float l2_phi_sys = 0;
float l2_d0sig = 0;
float l2_z0 = 0;

int ntrk_perp = 0;
int ntrk = 0;
float trk_pt[10000];
float trk_eta[10000];
float trk_phi[10000];
float trk_charge[10000];
float trk_d0[10000];
float trk_z0[10000];
bool  trk_truth_matched[10000];

int truth_ntrk = 0;
float truth_trk_pt[10000];
float truth_trk_eta[10000];
float truth_trk_phi[10000];
float truth_trk_charge[10000];

int njet = 0;
float jet_pt[40];
float jet_eta[40];
float jet_phi[40];
float jet_e[40];

int truth_jet_n = 0;
float truth_jet_pt[40];
float truth_jet_eta[40];
float truth_jet_phi[40];
float truth_jet_e[40];

struct OutTree {
  private:
  bool branchEventInfo = false;
  bool branchLeptons = false;
  bool branchZs = false;
  bool branchTruthZs = false;
  bool branchJets = false;
  bool branchTracks = false;
  bool branchTruthTracks = false;
  bool branchTruthJets = false;

  public:
  TTree* tree = nullptr;

  OutTree () {
    tree = nullptr;
  }
  OutTree (const char* name, TFile* file) {
    tree = new TTree (name, name);
    tree->SetDirectory (file);
  }

  void SetBranchEventInfo   (const bool _branchEventInfo = true)    { branchEventInfo = _branchEventInfo; }
  void SetBranchLeptons     (const bool _branchLeptons = true)      { branchLeptons = _branchLeptons; }
  void SetBranchZs          (const bool _branchZs = true)           { branchZs = _branchZs; }
  void SetBranchTruthZs     (const bool _branchTruthZs = true)      { branchTruthZs = _branchTruthZs; }
  void SetBranchJets        (const bool _branchJets = true)         { branchJets = _branchJets; }
  void SetBranchTracks      (const bool _branchTracks = true)       { branchTracks = _branchTracks; }
  void SetBranchTruthTracks (const bool _branchTruthTracks = true)  { branchTruthTracks = _branchTruthTracks; }
  void SetBranchTruthJets   (const bool _branchTruthJets = true)    { branchTruthJets = _branchTruthJets; }

  void SetBranches () {
    if (branchEventInfo) {
      tree->Branch ("event_number",  &event_number,  "event_number/I");
      tree->Branch ("run_number",    &run_number,    "run_number/I");
      tree->Branch ("event_weight",  &event_weight,  "event_weight/F");
      tree->Branch ("ntrk_perp",     &ntrk_perp,     "ntrk_perp/I");
      tree->Branch ("fcal_et",       &fcal_et,       "fcal_et/F");
      tree->Branch ("zdcEnergy",     &zdcEnergy,     "zdcEnergy/F");
      tree->Branch ("q2x_a",         &q2x_a,         "q2x_a/F");
      tree->Branch ("q2y_a",         &q2y_a,         "q2y_a/F");
      tree->Branch ("q2x_c",         &q2x_c,         "q2x_c/F");
      tree->Branch ("q2y_c",         &q2y_c,         "q2y_c/F");
      tree->Branch ("q2",            &q2,            "q2/F");
      tree->Branch ("psi2",          &psi2,          "psi2/F");
      tree->Branch ("q3",            &q3,            "q3/F");
      tree->Branch ("psi3",          &psi3,          "psi3/F");
      tree->Branch ("q4",            &q4,            "q4/F");
      tree->Branch ("psi4",          &psi4,          "psi4/F");
      tree->Branch ("vz",            &vz,            "vz/F");
    }
    if (branchLeptons) {
      tree->Branch ("l1_pt",         &l1_pt,         "l1_pt/F");
      tree->Branch ("l1_eta",        &l1_eta,        "l1_eta/F");
      tree->Branch ("l1_phi",        &l1_phi,        "l1_phi/F");
      tree->Branch ("l1_charge",     &l1_charge,     "l1_charge/I");
      tree->Branch ("l1_trk_pt",     &l1_trk_pt,     "l1_trk_pt/F");
      tree->Branch ("l1_trk_eta",    &l1_trk_eta,    "l1_trk_eta/F");
      tree->Branch ("l1_trk_phi",    &l1_trk_phi,    "l1_trk_phi/F");
      tree->Branch ("l1_pt_sys",     &l1_pt_sys,     "l1_pt_sys/F");
      tree->Branch ("l1_eta_sys",    &l1_eta_sys,    "l1_eta_sys/F");
      tree->Branch ("l1_phi_sys",    &l1_phi_sys,    "l1_phi_sys/F");
      tree->Branch ("l1_d0sig",      &l1_d0sig,      "l1_d0sig/F");
      tree->Branch ("l1_z0",         &l1_z0,         "l1_z0/F");

      tree->Branch ("l2_pt",         &l2_pt,         "l2_pt/F");
      tree->Branch ("l2_eta",        &l2_eta,        "l2_eta/F");
      tree->Branch ("l2_phi",        &l2_phi,        "l2_phi/F");
      tree->Branch ("l2_charge",     &l2_charge,     "l2_charge/I");
      tree->Branch ("l2_trk_pt",     &l2_trk_pt,     "l2_trk_pt/F");
      tree->Branch ("l2_trk_eta",    &l2_trk_eta,    "l2_trk_eta/F");
      tree->Branch ("l2_trk_phi",    &l2_trk_phi,    "l2_trk_phi/F");
      tree->Branch ("l2_pt_sys",     &l2_pt_sys,     "l2_pt_sys/F");
      tree->Branch ("l2_eta_sys",    &l2_eta_sys,    "l2_eta_sys/F");
      tree->Branch ("l2_phi_sys",    &l2_phi_sys,    "l2_phi_sys/F");
      tree->Branch ("l2_d0sig",      &l2_d0sig,      "l2_d0sig/F");
      tree->Branch ("l2_z0",         &l2_z0,         "l2_z0/F");
    }
    if (branchZs) {
      tree->Branch ("isEE",          &isEE,          "isEE/O");
      tree->Branch ("isSameSign",    &isSameSign,    "isSameSign/O");
      tree->Branch ("isQCDBkg",      &isQCDBkg,      "isQCDBkg/O");
      tree->Branch ("z_pt",          &z_pt,          "z_pt/F");
      tree->Branch ("z_y",           &z_y,           "z_y/F");
      tree->Branch ("z_phi",         &z_phi,         "z_phi/F");
      tree->Branch ("z_m",           &z_m,           "z_m/F");
      tree->Branch ("phi_transmin",  &phi_transmin,  "phi_transmin/F");
      tree->Branch ("phi_transmax",  &phi_transmax,  "phi_transmax/F");
    }

    if (branchTruthZs) {
      tree->Branch ("truth_z_pt",    &truth_z_pt,    "truth_z_pt/F");
      tree->Branch ("truth_z_y",     &truth_z_y,     "truth_z_y/F");
      tree->Branch ("truth_z_phi",   &truth_z_phi,   "truth_z_phi/F");
      tree->Branch ("truth_z_m",     &truth_z_m,     "truth_z_m/F");
    }

    if (branchJets) {
      tree->Branch ("njet",     &njet);
      tree->Branch ("jet_pt",   &jet_pt,  "jet_pt[njet]/F");
      tree->Branch ("jet_eta",  &jet_eta, "jet_eta[njet]/F");
      tree->Branch ("jet_phi",  &jet_phi, "jet_phi[njet]/F");
      tree->Branch ("jet_e",    &jet_e,   "jet_e[njet]/F");
    }

    if (branchTracks) {
      tree->Branch ("ntrk",               &ntrk,              "ntrk/I");
      tree->Branch ("trk_pt",             &trk_pt,            "trk_pt[ntrk]/F");
      tree->Branch ("trk_eta",            &trk_eta,           "trk_eta[ntrk]/F");
      tree->Branch ("trk_phi",            &trk_phi,           "trk_phi[ntrk]/F");
      tree->Branch ("trk_charge",         &trk_charge,        "trk_charge[ntrk]/F");
      tree->Branch ("trk_d0",             &trk_d0,            "trk_d0[ntrk]/F");
      tree->Branch ("trk_z0",             &trk_z0,            "trk_z0[ntrk]/F");
      tree->Branch ("trk_truth_matched",  &trk_truth_matched, "trk_truth_matched[ntrk]/O");
    }

    if (branchTruthTracks) {
      tree->Branch ("truth_ntrk",         &truth_ntrk,        "truth_ntrk/I");
      tree->Branch ("truth_trk_pt",       &truth_trk_pt,      "truth_trk_pt[truth_ntrk]/F");
      tree->Branch ("truth_trk_eta",      &truth_trk_eta,     "truth_trk_eta[truth_ntrk]/F");
      tree->Branch ("truth_trk_phi",      &truth_trk_phi,     "truth_trk_phi[truth_ntrk]/F");
      tree->Branch ("truth_trk_charge",   &truth_trk_charge,  "truth_trk_charge[truth_ntrk]/F");
    }

    if (branchTruthJets) {
      tree->Branch ("truth_jet_n",    &truth_jet_n,   "truth_jet_n/I");
      tree->Branch ("truth_jet_pt",   &truth_jet_pt,  "truth_jet_pt[truth_jet_n]/F");
      tree->Branch ("truth_jet_eta",  &truth_jet_eta, "truth_jet_eta[truth_jet_n]/F");
      tree->Branch ("truth_jet_phi",  &truth_jet_phi, "truth_jet_phi[truth_jet_n]/F");
      tree->Branch ("truth_jet_e",    &truth_jet_e,   "truth_jet_e[truth_jet_n]/F");
    }
    return;
  }

  void Fill () {
    tree->Fill ();
    return;
  }
};

}

#endif
