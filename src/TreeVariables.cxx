#ifndef __TreeVariables_cxx__
#define __TreeVariables_cxx__

#include "TreeVariables.h"
#include "ZTrackUtilities.h"

#include <iostream>
#include <iomanip>

namespace ZTrackAnalyzer {

TreeVariables :: TreeVariables (TTree* t, const bool _isMC) {
  tree = t;
  isMC = _isMC;
}


void TreeVariables :: SetBranchAddresses () {

  tree->SetBranchAddress ("run_number",     &run_number);
  tree->SetBranchAddress ("event_number",   &event_number);
  tree->SetBranchAddress ("lumi_block",     &lumi_block);
  tree->SetBranchAddress ("isOOTPU",        &isOOTPU);
  tree->SetBranchAddress ("BlayerDesyn",    &BlayerDesyn);
  tree->SetBranchAddress ("passes_toroid",  &passes_toroid);

  if (getCollisionRateInfo) {
    tree->SetBranchAddress ("actualInteractionsPerCrossing",  &actualInteractionsPerCrossing);
    tree->SetBranchAddress ("averageInteractionsPerCrossing", &averageInteractionsPerCrossing);
  } else {
    tree->SetBranchStatus ("actualInteractionsPerCrossing",  0);
    tree->SetBranchStatus ("averageInteractionsPerCrossing", 0);
  }

  if (isMC) {
    tree->SetBranchAddress ("mcEventWeights", &mcEventWeights);
  }

  if (getVertices) {
    tree->SetBranchAddress ("nvert",     &nvert);
    tree->SetBranchAddress ("vert_x",    vert_x);
    tree->SetBranchAddress ("vert_y",    vert_y);
    tree->SetBranchAddress ("vert_z",    vert_z);
    tree->SetBranchAddress ("vert_ntrk", vert_ntrk);
    tree->SetBranchAddress ("vert_type", vert_type);
  } else {
    tree->SetBranchStatus ("nvert",     0);
    tree->SetBranchStatus ("vert_x",    0);
    tree->SetBranchStatus ("vert_y",    0);
    tree->SetBranchStatus ("vert_z",    0);
    tree->SetBranchStatus ("vert_ntrk", 0);
    tree->SetBranchStatus ("vert_type", 0);
  }

  if (getFCals) {
    tree->SetBranchAddress ("fcalA_et",     &fcalA_et);
    tree->SetBranchAddress ("fcalC_et",     &fcalC_et);
    tree->SetBranchAddress ("fcalA_et_Cos", &fcalA_et_Cos);
    tree->SetBranchAddress ("fcalC_et_Cos", &fcalC_et_Cos);
    tree->SetBranchAddress ("fcalA_et_Sin", &fcalA_et_Sin);
    tree->SetBranchAddress ("fcalC_et_Sin", &fcalC_et_Sin);
  } else {
    tree->SetBranchStatus ("fcalA_et",     0);
    tree->SetBranchStatus ("fcalC_et",     0);
    tree->SetBranchStatus ("fcalA_et_Cos", 0);
    tree->SetBranchStatus ("fcalC_et_Cos", 0);
    tree->SetBranchStatus ("fcalA_et_Sin", 0);
    tree->SetBranchStatus ("fcalC_et_Sin", 0);
  }

  if (getZdc) {
    tree->SetBranchAddress ("ZdcCalibEnergy_A",  &ZdcCalibEnergy_A);
    tree->SetBranchAddress ("ZdcCalibEnergy_C",  &ZdcCalibEnergy_C);
    tree->SetBranchAddress ("L1_ZDC_A",          &L1_ZDC_A);
    tree->SetBranchAddress ("L1_ZDC_A_tbp",      &L1_ZDC_A_tbp);
    tree->SetBranchAddress ("L1_ZDC_A_tap",      &L1_ZDC_A_tap);
    tree->SetBranchAddress ("L1_ZDC_A_tav",      &L1_ZDC_A_tav);
    tree->SetBranchAddress ("L1_ZDC_A_prescale", &L1_ZDC_A_prescale);
    tree->SetBranchAddress ("L1_ZDC_C",          &L1_ZDC_C);
    tree->SetBranchAddress ("L1_ZDC_C_tbp",      &L1_ZDC_C_tbp);
    tree->SetBranchAddress ("L1_ZDC_C_tap",      &L1_ZDC_C_tap);
    tree->SetBranchAddress ("L1_ZDC_C_tav",      &L1_ZDC_C_tav);
    tree->SetBranchAddress ("L1_ZDC_C_prescale", &L1_ZDC_C_prescale);
  } else {
    tree->SetBranchStatus ("ZdcCalibEnergy_A",  0);
    tree->SetBranchStatus ("ZdcCalibEnergy_C",  0);
    tree->SetBranchStatus ("L1_ZDC_A",          0);
    tree->SetBranchStatus ("L1_ZDC_A_tbp",      0);
    tree->SetBranchStatus ("L1_ZDC_A_tap",      0);
    tree->SetBranchStatus ("L1_ZDC_A_tav",      0);
    tree->SetBranchStatus ("L1_ZDC_A_prescale", 0);
    tree->SetBranchStatus ("L1_ZDC_C",          0);
    tree->SetBranchStatus ("L1_ZDC_C_tbp",      0);
    tree->SetBranchStatus ("L1_ZDC_C_tap",      0);
    tree->SetBranchStatus ("L1_ZDC_C_tav",      0);
    tree->SetBranchStatus ("L1_ZDC_C_prescale", 0);
  }

  if (getTracks) {
    tree->SetBranchAddress ("ntrk",               &ntrk);
    tree->SetBranchAddress ("trk_pt",             &trk_pt);
    tree->SetBranchAddress ("trk_eta",            &trk_eta);
    tree->SetBranchAddress ("trk_phi",            &trk_phi);
    tree->SetBranchAddress ("trk_charge",         &trk_charge);
    tree->SetBranchAddress ("trk_HItight",        &trk_HItight);
    tree->SetBranchAddress ("trk_HIloose",        &trk_HIloose);
    tree->SetBranchAddress ("trk_d0",             &trk_d0);
    tree->SetBranchAddress ("trk_z0",             &trk_z0);
    tree->SetBranchAddress ("trk_vz",             &trk_vz);
    tree->SetBranchAddress ("trk_theta",          &trk_theta);
    if (isMC) {
      tree->SetBranchAddress ("trk_prob_truth",     &trk_prob_truth);
      tree->SetBranchAddress ("trk_truth_pt",       &trk_truth_pt);
      tree->SetBranchAddress ("trk_truth_eta",      &trk_truth_eta);
      tree->SetBranchAddress ("trk_truth_phi",      &trk_truth_phi);
      tree->SetBranchAddress ("trk_truth_charge",   &trk_truth_charge);
      tree->SetBranchAddress ("trk_truth_type",     &trk_truth_type);
      tree->SetBranchAddress ("trk_truth_orig",     &trk_truth_orig);
      tree->SetBranchAddress ("trk_truth_barcode",  &trk_truth_barcode);
      tree->SetBranchAddress ("trk_truth_pdgid",    &trk_truth_pdgid);
      tree->SetBranchAddress ("trk_truth_vz",       &trk_truth_vz);
      tree->SetBranchAddress ("trk_truth_nIn",      &trk_truth_nIn);
      tree->SetBranchAddress ("trk_truth_isHadron", &trk_truth_isHadron);
    }
  } else {
    tree->SetBranchStatus ("ntrk",                0);
    tree->SetBranchStatus ("trk_pt",              0);
    tree->SetBranchStatus ("trk_eta",             0);
    tree->SetBranchStatus ("trk_phi",             0);
    tree->SetBranchStatus ("trk_charge",          0);
    tree->SetBranchStatus ("trk_HItight",         0);
    tree->SetBranchStatus ("trk_HIloose",         0);
    tree->SetBranchStatus ("trk_d0",              0);
    tree->SetBranchStatus ("trk_z0",              0);
    tree->SetBranchStatus ("trk_vz",              0);
    tree->SetBranchStatus ("trk_theta",           0);
    if (isMC) {
      tree->SetBranchStatus ("trk_prob_truth",      0);
      tree->SetBranchStatus ("trk_truth_pt",        0);
      tree->SetBranchStatus ("trk_truth_eta",       0);
      tree->SetBranchStatus ("trk_truth_phi",       0);
      tree->SetBranchStatus ("trk_truth_charge",    0);
      tree->SetBranchStatus ("trk_truth_type",      0);
      tree->SetBranchStatus ("trk_truth_orig",      0);
      tree->SetBranchStatus ("trk_truth_barcode",   0);
      tree->SetBranchStatus ("trk_truth_pdgid",     0);
      tree->SetBranchStatus ("trk_truth_vz",        0);
      tree->SetBranchStatus ("trk_truth_nIn",       0);
      tree->SetBranchStatus ("trk_truth_isHadron",  0);
    }
  }

  if (getElectrons) {
    tree->SetBranchAddress ("electron_n",                   &electron_n);
    tree->SetBranchAddress ("electron_pt_precalib",         &electron_pt_precalib);
    tree->SetBranchAddress ("electron_pt",                  &electron_pt);
    tree->SetBranchAddress ("electron_eta",                 &electron_eta);
    tree->SetBranchAddress ("electron_phi",                 &electron_phi);
    tree->SetBranchAddress ("electron_charge",              &electron_charge);
    tree->SetBranchAddress ("electron_lhtight",             &electron_lhtight);
    tree->SetBranchAddress ("electron_lhmedium",            &electron_lhmedium);
    tree->SetBranchAddress ("electron_lhloose",             &electron_lhloose);
    tree->SetBranchAddress ("electron_lhmedium_hi",         &electron_lhmedium_hi);
    tree->SetBranchAddress ("electron_lhloose_hi",          &electron_lhloose_hi);
    tree->SetBranchAddress ("electron_matched",             &electron_matched);
    tree->SetBranchAddress ("electron_etcone20",            &electron_etcone20);
    tree->SetBranchAddress ("electron_etcone30",            &electron_etcone30);
    tree->SetBranchAddress ("electron_etcone40",            &electron_etcone40);
    tree->SetBranchAddress ("electron_topoetcone20",        &electron_topoetcone20);
    tree->SetBranchAddress ("electron_topoetcone30",        &electron_topoetcone30);
    tree->SetBranchAddress ("electron_topoetcone40",        &electron_topoetcone40);
    tree->SetBranchAddress ("electron_id_track_pt",         &electron_id_track_pt);
    tree->SetBranchAddress ("electron_id_track_eta",        &electron_id_track_eta);
    tree->SetBranchAddress ("electron_id_track_phi",        &electron_id_track_phi);
    tree->SetBranchAddress ("electron_id_track_charge",     &electron_id_track_charge);
    tree->SetBranchAddress ("electron_id_track_d0sig",      &electron_id_track_d0sig);
    tree->SetBranchAddress ("electron_id_track_z0",         &electron_id_track_z0);
    tree->SetBranchAddress ("electron_id_track_vz",         &electron_id_track_vz);
    tree->SetBranchAddress ("electron_id_track_theta",      &electron_id_track_theta);
    tree->SetBranchAddress ("electron_pt_sys",              &electron_pt_sys);
    tree->SetBranchAddress ("electron_eta_sys",             &electron_eta_sys);
    tree->SetBranchAddress ("electron_phi_sys",             &electron_phi_sys);
  } else {
    tree->SetBranchStatus ("electron_n",                  0);
    tree->SetBranchStatus ("electron_pt_precalib",        0);
    tree->SetBranchStatus ("electron_pt",                 0);
    tree->SetBranchStatus ("electron_eta",                0);
    tree->SetBranchStatus ("electron_phi",                0);
    tree->SetBranchStatus ("electron_charge",             0);
    tree->SetBranchStatus ("electron_lhtight",            0);
    tree->SetBranchStatus ("electron_lhmedium",           0);
    tree->SetBranchStatus ("electron_lhloose",            0);
    tree->SetBranchStatus ("electron_lhmedium_hi",        0);
    tree->SetBranchStatus ("electron_lhloose_hi",         0);
    tree->SetBranchStatus ("electron_matched",            0);
    tree->SetBranchStatus ("electron_etcone20",           0);
    tree->SetBranchStatus ("electron_etcone30",           0);
    tree->SetBranchStatus ("electron_etcone40",           0);
    tree->SetBranchStatus ("electron_topoetcone20",       0);
    tree->SetBranchStatus ("electron_topoetcone30",       0);
    tree->SetBranchStatus ("electron_topoetcone40",       0);
    tree->SetBranchStatus ("electron_ntrk",               0);
    tree->SetBranchStatus ("electron_id_track_pt",        0);
    tree->SetBranchStatus ("electron_id_track_eta",       0);
    tree->SetBranchStatus ("electron_id_track_phi",       0);
    tree->SetBranchStatus ("electron_id_track_charge",    0);
    tree->SetBranchStatus ("electron_id_track_d0sig",     0);
    tree->SetBranchStatus ("electron_id_track_z0",        0);
    tree->SetBranchStatus ("electron_id_track_vz",        0);
    tree->SetBranchStatus ("electron_id_track_theta",     0);
    tree->SetBranchStatus ("electron_pt_sys",             0);
    tree->SetBranchStatus ("electron_eta_sys",            0);
    tree->SetBranchStatus ("electron_phi_sys",            0);
  }

  if (getClusters) {
    tree->SetBranchAddress ("Cluster_n",        &Cluster_n);
    tree->SetBranchAddress ("Cluster_pt",       &Cluster_pt);
    tree->SetBranchAddress ("Cluster_et",       &Cluster_et);
    tree->SetBranchAddress ("Cluster_eta",      &Cluster_eta);
    tree->SetBranchAddress ("Cluster_phi",      &Cluster_phi);
    tree->SetBranchAddress ("Cluster_energyBE", &Cluster_energyBE);
    tree->SetBranchAddress ("Cluster_etaBE",    &Cluster_etaBE);
    tree->SetBranchAddress ("Cluster_phiBE",    &Cluster_phiBE);
    tree->SetBranchAddress ("Cluster_calE",     &Cluster_calE);
    tree->SetBranchAddress ("Cluster_calEta",   &Cluster_calEta);
    tree->SetBranchAddress ("Cluster_calPhi",   &Cluster_calPhi);
    tree->SetBranchAddress ("Cluster_size",     &Cluster_size);
    tree->SetBranchAddress ("Cluster_status",   &Cluster_status);
  } else {
    tree->SetBranchStatus ("Cluster_n",         0);
    tree->SetBranchStatus ("Cluster_pt",        0);
    tree->SetBranchStatus ("Cluster_et",        0);
    tree->SetBranchStatus ("Cluster_eta",       0);
    tree->SetBranchStatus ("Cluster_phi",       0);
    tree->SetBranchStatus ("Cluster_energyBE",  0);
    tree->SetBranchStatus ("Cluster_etaBE",     0);
    tree->SetBranchStatus ("Cluster_phiBE",     0);
    tree->SetBranchStatus ("Cluster_calE",      0);
    tree->SetBranchStatus ("Cluster_calEta",    0);
    tree->SetBranchStatus ("Cluster_calPhi",    0);
    tree->SetBranchStatus ("Cluster_size",      0);
    tree->SetBranchStatus ("Cluster_status",    0);
  }

  if (getMuons) {
    tree->SetBranchAddress ("muon_n",                       &muon_n);
    tree->SetBranchAddress ("muon_pt_precalib",             &muon_pt_precalib);
    tree->SetBranchAddress ("muon_pt",                      &muon_pt);
    tree->SetBranchAddress ("muon_ms_pt_precalib",          &muon_ms_pt_precalib);
    tree->SetBranchAddress ("muon_ms_pt",                   &muon_ms_pt);
    tree->SetBranchAddress ("muon_eta",                     &muon_eta);
    tree->SetBranchAddress ("muon_phi",                     &muon_phi);
    tree->SetBranchAddress ("muon_charge",                  &muon_charge);
    tree->SetBranchAddress ("muon_tight",                   &muon_tight);
    tree->SetBranchAddress ("muon_medium",                  &muon_medium);
    tree->SetBranchAddress ("muon_loose",                   &muon_loose);
    tree->SetBranchAddress ("muon_matched",                 &muon_matched);
    tree->SetBranchAddress ("muon_etcone20",                &muon_etcone20);
    tree->SetBranchAddress ("muon_etcone30",                &muon_etcone30);
    tree->SetBranchAddress ("muon_etcone40",                &muon_etcone40);
    tree->SetBranchAddress ("muon_topoetcone20",            &muon_topoetcone20);
    tree->SetBranchAddress ("muon_topoetcone30",            &muon_topoetcone30);
    tree->SetBranchAddress ("muon_topoetcone40",            &muon_topoetcone40);
    tree->SetBranchAddress ("muon_id_track_pt",             &muon_id_track_pt);
    tree->SetBranchAddress ("muon_id_track_eta",            &muon_id_track_eta);
    tree->SetBranchAddress ("muon_id_track_phi",            &muon_id_track_phi);
    tree->SetBranchAddress ("muon_id_track_charge",         &muon_id_track_charge);
    tree->SetBranchAddress ("muon_id_track_d0sig",          &muon_id_track_d0sig);
    tree->SetBranchAddress ("muon_id_track_z0",             &muon_id_track_z0);
    tree->SetBranchAddress ("muon_id_track_vz",             &muon_id_track_vz);
    tree->SetBranchAddress ("muon_id_track_theta",          &muon_id_track_theta);
    tree->SetBranchAddress ("muon_id_track_tightprimary",   &muon_id_track_tightprimary);
    tree->SetBranchAddress ("muon_id_track_hiloose",        &muon_id_track_hiloose);
    tree->SetBranchAddress ("muon_id_track_hitight",        &muon_id_track_hitight);
    tree->SetBranchAddress ("muon_ms_track_pt",             &muon_ms_track_pt);
    tree->SetBranchAddress ("muon_ms_track_eta",            &muon_ms_track_eta);
    tree->SetBranchAddress ("muon_ms_track_phi",            &muon_ms_track_phi);
    tree->SetBranchAddress ("muon_ms_track_charge",         &muon_ms_track_charge);
    tree->SetBranchAddress ("muon_ms_track_d0sig",          &muon_ms_track_d0sig);
    tree->SetBranchAddress ("muon_ms_track_z0",             &muon_ms_track_z0);
    tree->SetBranchAddress ("muon_ms_track_vz",             &muon_ms_track_vz);
    tree->SetBranchAddress ("muon_ms_track_theta",          &muon_ms_track_theta);
    tree->SetBranchAddress ("muon_pt_sys",                  &muon_pt_sys);
    tree->SetBranchAddress ("muon_eta_sys",                 &muon_eta_sys);
    tree->SetBranchAddress ("muon_phi_sys",                 &muon_phi_sys);
  } else {
    tree->SetBranchStatus ("muon_n",                      0);
    tree->SetBranchStatus ("muon_pt_precalib",            0);
    tree->SetBranchStatus ("muon_pt",                     0);
    tree->SetBranchStatus ("muon_ms_pt_precalib",         0);
    tree->SetBranchStatus ("muon_ms_pt",                  0);
    tree->SetBranchStatus ("muon_eta",                    0);
    tree->SetBranchStatus ("muon_phi",                    0);
    tree->SetBranchStatus ("muon_charge",                 0);
    tree->SetBranchStatus ("muon_tight",                  0);
    tree->SetBranchStatus ("muon_medium",                 0);
    tree->SetBranchStatus ("muon_loose",                  0);
    tree->SetBranchStatus ("muon_matched",                0);
    tree->SetBranchStatus ("muon_etcone20",               0);
    tree->SetBranchStatus ("muon_etcone30",               0);
    tree->SetBranchStatus ("muon_etcone40",               0);
    tree->SetBranchStatus ("muon_topoetcone20",           0);
    tree->SetBranchStatus ("muon_topoetcone30",           0);
    tree->SetBranchStatus ("muon_topoetcone40",           0);
    tree->SetBranchStatus ("muon_id_track_pt",            0);
    tree->SetBranchStatus ("muon_id_track_eta",           0);
    tree->SetBranchStatus ("muon_id_track_phi",           0);
    tree->SetBranchStatus ("muon_id_track_charge",        0);
    tree->SetBranchStatus ("muon_id_track_d0",            0);
    tree->SetBranchStatus ("muon_id_track_z0",            0);
    tree->SetBranchStatus ("muon_id_track_vz",            0);
    tree->SetBranchStatus ("muon_id_track_theta",         0);
    tree->SetBranchStatus ("muon_id_track_tightprimary",  0);
    tree->SetBranchStatus ("muon_id_track_hiloose",       0);
    tree->SetBranchStatus ("muon_id_track_hitight",       0);
    tree->SetBranchStatus ("muon_ms_track_pt",            0);
    tree->SetBranchStatus ("muon_ms_track_eta",           0);
    tree->SetBranchStatus ("muon_ms_track_phi",           0);
    tree->SetBranchStatus ("muon_ms_track_charge",        0);
    tree->SetBranchStatus ("muon_ms_track_d0",            0);
    tree->SetBranchStatus ("muon_ms_track_z0",            0);
    tree->SetBranchStatus ("muon_ms_track_vz",            0);
    tree->SetBranchStatus ("muon_ms_track_theta",         0);
    tree->SetBranchStatus ("muon_pt_sys",                 0);
    tree->SetBranchStatus ("muon_eta_sys",                0);
    tree->SetBranchStatus ("muon_phi_sys",                0);
  }

  if (getJets) {
    tree->SetBranchAddress ("akt4emtopo_jet_n",   &akt4emtopo_jet_n);
    tree->SetBranchAddress ("akt4emtopo_jet_pt",  &akt4emtopo_jet_pt);
    tree->SetBranchAddress ("akt4emtopo_jet_eta", &akt4emtopo_jet_eta);
    tree->SetBranchAddress ("akt4emtopo_jet_phi", &akt4emtopo_jet_phi);
    tree->SetBranchAddress ("akt4emtopo_jet_e",   &akt4emtopo_jet_e);
  } else {
    tree->SetBranchStatus ("akt4emtopo_jet_n",    0);
    tree->SetBranchStatus ("akt4emtopo_jet_pt",   0);
    tree->SetBranchStatus ("akt4emtopo_jet_eta",  0);
    tree->SetBranchStatus ("akt4emtopo_jet_phi",  0);
    tree->SetBranchStatus ("akt4emtopo_jet_e",    0);
  }

  if (isMC) {

    if (getTruthElectrons) {
      tree->SetBranchAddress ("truth_electron_n",       &truth_electron_n);
      tree->SetBranchAddress ("truth_electron_pt",      &truth_electron_pt);
      tree->SetBranchAddress ("truth_electron_eta",     &truth_electron_eta);
      tree->SetBranchAddress ("truth_electron_phi",     &truth_electron_phi);
      tree->SetBranchAddress ("truth_electron_charge",  &truth_electron_charge);
      tree->SetBranchAddress ("truth_electron_barcode", &truth_electron_barcode);
    } else {
      tree->SetBranchStatus ("truth_electron_n",        0);
      tree->SetBranchStatus ("truth_electron_pt",       0);
      tree->SetBranchStatus ("truth_electron_eta",      0);
      tree->SetBranchStatus ("truth_electron_phi",      0);
      tree->SetBranchStatus ("truth_electron_charge",   0);
      tree->SetBranchStatus ("truth_electron_barcode",  0);
    }

    if (getTruthMuons) {
      tree->SetBranchAddress ("truth_muon_n",       &truth_muon_n);
      tree->SetBranchAddress ("truth_muon_pt",      &truth_muon_pt);
      tree->SetBranchAddress ("truth_muon_eta",     &truth_muon_eta);
      tree->SetBranchAddress ("truth_muon_phi",     &truth_muon_phi);
      tree->SetBranchAddress ("truth_muon_charge",  &truth_muon_charge);
      tree->SetBranchAddress ("truth_muon_barcode", &truth_muon_barcode);
    } else {
      tree->SetBranchStatus ("truth_muon_n",        0);
      tree->SetBranchStatus ("truth_muon_pt",       0);
      tree->SetBranchStatus ("truth_muon_eta",      0);
      tree->SetBranchStatus ("truth_muon_phi",      0);
      tree->SetBranchStatus ("truth_muon_charge",   0);
      tree->SetBranchStatus ("truth_muon_barcode",  0);
    }

    if (getTruthTracks) {
      tree->SetBranchAddress ("truth_trk_n",        &truth_trk_n);
      tree->SetBranchAddress ("truth_trk_pt",       &truth_trk_pt);
      tree->SetBranchAddress ("truth_trk_eta",      &truth_trk_eta);
      tree->SetBranchAddress ("truth_trk_phi",      &truth_trk_phi);
      tree->SetBranchAddress ("truth_trk_charge",   &truth_trk_charge);
      tree->SetBranchAddress ("truth_trk_pdgid",    &truth_trk_pdgid);
      tree->SetBranchAddress ("truth_trk_barcode",  &truth_trk_barcode);
      tree->SetBranchAddress ("truth_trk_isHadron", &truth_trk_isHadron);
    } else {
      tree->SetBranchStatus ("truth_trk_n",         0);
      tree->SetBranchStatus ("truth_trk_pt",        0);
      tree->SetBranchStatus ("truth_trk_eta",       0);
      tree->SetBranchStatus ("truth_trk_phi",       0);
      tree->SetBranchStatus ("truth_trk_charge",    0);
      tree->SetBranchStatus ("truth_trk_pdgid",     0);
      tree->SetBranchStatus ("truth_trk_barcode",   0);
      tree->SetBranchStatus ("truth_trk_isHadron",  0);
    }

    if (getTruthJets) {
      tree->SetBranchAddress ("akt4_truth_jet_n",   &truth_jet_n);
      tree->SetBranchAddress ("akt4_truth_jet_pt",  &truth_jet_pt);
      tree->SetBranchAddress ("akt4_truth_jet_eta", &truth_jet_eta);
      tree->SetBranchAddress ("akt4_truth_jet_phi", &truth_jet_phi);
      tree->SetBranchAddress ("akt4_truth_jet_e",   &truth_jet_e);
    } else {
      tree->SetBranchStatus ("akt4_truth_jet_n",   0);
      tree->SetBranchStatus ("akt4_truth_jet_pt",  0);
      tree->SetBranchStatus ("akt4_truth_jet_eta", 0);
      tree->SetBranchStatus ("akt4_truth_jet_phi", 0);
      tree->SetBranchStatus ("akt4_truth_jet_e",   0);
    }
  }

  return;
}

void TreeVariables :: PrintAll (const long long entry) {
  GetEntry (entry);
  cout << endl << "////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "// Getting event " << entry << "..." << endl;
  cout << "////////////////////////////////////////////////////////////////////////////////" << endl;

  if (getElectrons) {
    cout << endl << "Electrons:" << endl;
    cout << setw (2) << "#"
         << setw (12) << "Pt"
         << setw (12) << "Eta"
         << setw (12) << "Phi"
         << setw (12) << "Charge" << endl;
    for (int j = 0; j < electron_n; j++) {
      cout << setw (2) << j
           << setw (12) << electron_pt[j]
           << setw (12) << electron_eta[j]
           << setw (12) << electron_phi[j]
           << setw (12) << electron_charge[j] << endl;
    }
  }

  if (isMC && getTruthElectrons) {
    cout << endl << "Truth electrons:" << endl;
    cout << setw (2) << "#"
         << setw (12) << "Pt"
         << setw (12) << "Eta"
         << setw (12) << "Phi"
         << setw (12) << "Charge" << endl;
    for (int j = 0; j < truth_electron_n; j++) {
      cout << setw (2) << j
           << setw (12) << truth_electron_pt[j]
           << setw (12) << truth_electron_eta[j]
           << setw (12) << truth_electron_phi[j]
           << setw (12) << truth_electron_charge[j] << endl;
    }
  }

  if (getMuons) {
    cout << endl << "Muons:" << endl;
    cout << setw (2) << "#"
         << setw (12) << "Pt"
         << setw (12) << "Eta"
         << setw (12) << "Phi"
         << setw (12) << "Charge" << endl;
    for (int j = 0; j < muon_n; j++) {
      cout << setw (2) << j
           << setw (12) << muon_pt[j]
           << setw (12) << muon_eta[j]
           << setw (12) << muon_phi[j]
           << setw (12) << muon_charge[j] << endl;
    }
  }

  if (isMC && getTruthMuons) {
    cout << endl << "Truth muons:" << endl;
    cout << setw (2) << "#"
         << setw (12) << "Pt"
         << setw (12) << "Eta"
         << setw (12) << "Phi"
         << setw (12) << "Charge" << endl;
    for (int j = 0; j < truth_muon_n; j++) {
      cout << setw (2) << j
           << setw (12) << truth_muon_pt[j]
           << setw (12) << truth_muon_eta[j]
           << setw (12) << truth_muon_phi[j]
           << setw (12) << truth_muon_charge[j] << endl;
    }
  }

  if (getTracks) {
    cout << endl << "Tracks:" << endl;
    cout << setw (3) << "#"
         << setw (12) << "HI Tight"
         << setw (12) << "HI Loose"
         << setw (12) << "Pt"
         << setw (12) << "Eta"
         << setw (12) << "Phi";
    if (isMC) { 
      cout << setw (12) << "truth prob."
           << setw (12) << "pdgid";
    }
    cout << endl;

    for (int t = 0; t < ntrk; t++) {
      float minDR = 1;
      for (int iETrk = 0; iETrk < electron_n; iETrk++) {
          minDR = std::fmin (minDR, (float)DeltaR (electron_id_track_eta[iETrk], trk_eta[t], electron_id_track_phi[iETrk], trk_phi[t]));
      }
      for (int iMTrk = 0; iMTrk < muon_n; iMTrk++) {
          minDR = std::fmin (minDR, (float)DeltaR (muon_id_track_eta[iMTrk], trk_eta[t], muon_id_track_phi[iMTrk], trk_phi[t]));
      }
      if (!(0.05 < minDR && minDR < 0.1))
        continue;

      cout << setw (3) << t
           << setw (12) << trk_HIloose[t]
           << setw (12) << trk_HItight[t]
           << setw (12) << trk_pt[t]
           << setw (12) << trk_eta[t]
           << setw (12) << trk_phi[t];
      if (isMC) {
        cout << setw (12) << trk_prob_truth[t]
             << setw (12) << trk_truth_pdgid[t];
      }
      cout << endl;
    }
  }

  return;
}

/** Setter functions **/
void TreeVariables :: SetGetCollisionRateInfo (const bool _getCollisionRateInfo)  { getCollisionRateInfo  = _getCollisionRateInfo;  }
void TreeVariables :: SetGetVertices          (const bool _getVertices)           { getVertices           = _getVertices;           }
void TreeVariables :: SetGetFCals             (const bool _getFCals)              { getFCals              = _getFCals;              }
void TreeVariables :: SetGetZdc               (const bool _getZdc)                { getZdc                = _getZdc;                }
void TreeVariables :: SetGetTracks            (const bool _getTracks)             { getTracks             = _getTracks;             }
void TreeVariables :: SetGetElectrons         (const bool _getElectrons)          { getElectrons          = _getElectrons;          }
void TreeVariables :: SetGetClusters          (const bool _getClusters)           { getClusters           = _getClusters;           }
void TreeVariables :: SetGetMuons             (const bool _getMuons)              { getMuons              = _getMuons;              }
void TreeVariables :: SetGetJets              (const bool _getJets)               { getJets               = _getJets;               }
void TreeVariables :: SetGetTruthTracks       (const bool _getTruthTracks)        { getTruthTracks        = _getTruthTracks;        }
void TreeVariables :: SetGetTruthElectrons    (const bool _getTruthElectrons)     { getTruthElectrons     = _getTruthElectrons;     }
void TreeVariables :: SetGetTruthMuons        (const bool _getTruthMuons)         { getTruthMuons         = _getTruthMuons;         }
void TreeVariables :: SetGetTruthJets         (const bool _getTruthJets)          { getTruthJets          = _getTruthJets;          }

} // end namespace

#endif
