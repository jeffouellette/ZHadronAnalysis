#ifndef __OutTree_h__
#define __OutTree_h__

#include <TTree.h>
#include <vector>

using namespace std;

namespace ZTrackAnalyzer {

bool isMC = false;
bool isEE = false;
float event_weight = 0;

unsigned int event_number = 0;
unsigned int run_number = 0;

float fcal_et = 0;
float zdcEnergy = 0;
float q2 = 0;
float psi2 = 0;
float vz = 0;

float z_pt = 0;
float z_y = 0;
float z_phi = 0;
float z_m = 0;

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

int ntrk_all = 0;
int ntrk = 0;
vector<float> trk_pt (0);
vector<float> trk_eta (0);
vector<float> trk_phi (0);
vector<float> trk_charge (0);
vector<float> trk_z0 (0);
//vector<bool>  trk_hiloose (0);
//vector<bool>  trk_hitight (0);

int njet = 0;
vector<float> jet_pt (0);
vector<float> jet_eta (0);
vector<float> jet_phi (0);
vector<float> jet_e (0);

int truth_jet_n = 0;
vector<float> truth_jet_pt (0);
vector<float> truth_jet_eta (0);
vector<float> truth_jet_phi (0);
vector<float> truth_jet_e (0);

struct OutTree {
  private:
  bool branchEventInfo = false;
  bool branchLeptons = false;
  bool branchZs = false;
  bool branchTruthZs = false;
  bool branchJets = false;
  bool branchTracks = false;
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

  void SetBranchEventInfo (const bool _branchEventInfo = true) { branchEventInfo = _branchEventInfo; }
  void SetBranchLeptons   (const bool _branchLeptons = true)   { branchLeptons = _branchLeptons; }
  void SetBranchZs        (const bool _branchZs = true)        { branchZs = _branchZs; }
  void SetBranchTruthZs   (const bool _branchTruthZs = true)   { branchTruthZs = _branchTruthZs; }
  void SetBranchJets      (const bool _branchJets = true)      { branchJets = _branchJets; }
  void SetBranchTracks    (const bool _branchTracks = true)    { branchTracks = _branchTracks; }
  void SetBranchTruthJets (const bool _branchTruthJets = true) { branchTruthJets = _branchTruthJets; }

  void SetBranches () {
    if (branchEventInfo) {
      tree->Branch ("event_number",  &event_number,  "event_number/I");
      tree->Branch ("run_number",    &run_number,    "run_number/I");
      tree->Branch ("isEE",          &isEE,          "isEE/O");
      tree->Branch ("event_weight",  &event_weight,  "event_weight/F");
      tree->Branch ("ntrk_all",      &ntrk_all,      "ntrk_all/I");
      tree->Branch ("fcal_et",       &fcal_et,       "fcal_et/F");
      tree->Branch ("zdcEnergy",     &zdcEnergy,     "zdcEnergy/F");
      tree->Branch ("q2",            &q2,            "q2/F");
      tree->Branch ("psi2",          &psi2,          "psi2/F");
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
    }
    if (branchZs) {
      tree->Branch ("z_pt",          &z_pt,          "z_pt/F");
      tree->Branch ("z_y",           &z_y,           "z_y/F");
      tree->Branch ("z_phi",         &z_phi,         "z_phi/F");
      tree->Branch ("z_m",           &z_m,           "z_m/F");
    }
    if (branchTruthZs) {
      tree->Branch ("truth_z_pt",    &truth_z_pt,    "truth_z_pt/F");
      tree->Branch ("truth_z_y",     &truth_z_y,     "truth_z_y/F");
      tree->Branch ("truth_z_phi",   &truth_z_phi,   "truth_z_phi/F");
      tree->Branch ("truth_z_m",     &truth_z_m,     "truth_z_m/F");
    }
    if (branchJets) {
      tree->Branch ("njet",          &njet);
      tree->Branch ("jet_pt",        &jet_pt);
      tree->Branch ("jet_eta",       &jet_eta);
      tree->Branch ("jet_phi",       &jet_phi);
      tree->Branch ("jet_e",         &jet_e);
    }

    if (branchTracks) {
      tree->Branch ("ntrk",          &ntrk,          "ntrk/I");
      tree->Branch ("trk_pt",        &trk_pt);
      tree->Branch ("trk_eta",       &trk_eta);
      tree->Branch ("trk_phi",       &trk_phi);
      tree->Branch ("trk_charge",    &trk_charge);
      //tree->Branch ("trk_hiloose",   &trk_hiloose);
      //tree->Branch ("trk_hitight",   &trk_hitight);
    }
    if (branchTruthJets) {
      tree->Branch ("truth_jet_n",   &truth_jet_n,   "truth_jet_n/I");
      tree->Branch ("truth_jet_pt",  &truth_jet_pt);
      tree->Branch ("truth_jet_eta", &truth_jet_eta);
      tree->Branch ("truth_jet_phi", &truth_jet_phi);
      tree->Branch ("truth_jet_e",   &truth_jet_e);
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
