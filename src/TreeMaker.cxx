#include "TreeMaker.h"
#include "Params.h"
#include "TreeVariables.h"
#include "OutTree.h"
#include "Trigger.h"
#include "ZTrackUtilities.h"
#include "Cuts.h"

#include <TChain.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH2F.h>
#include <TLorentzVector.h>

#include <iostream>

using namespace std;

namespace ZTrackAnalyzer {

bool TreeMaker (const char* directory,
                const int dataSet,
                const char* inFileName) {
 
  cout << "Info: In TreeMaker.cxx: Entered TreeMaker routine." << endl;
  cout << "Info: In TreeMaker.cxx: Printing systematic onfiguration:";
  cout << "\n\tdoElectronPtUpVar: " << doElectronPtUpVar;
  cout << "\n\tdoElectronPtDownVar: " << doElectronPtDownVar;
  cout << "\n\tdoMuonPtUpVar: " << doMuonPtUpVar;
  cout << "\n\tdoMuonPtDownVar: " << doMuonPtDownVar;
  cout << "\n\tdoElectronLHMediumVar: " << doElectronLHMediumVar;
  cout << "\n\tdoMuonTightVar: " << doMuonTightVar;
  cout << "\n\tdoHITightVar: " << doHITightVar;
  cout << endl;

  if (isMC)
    SetupDirectories ("MCAnalysis");
  else
    SetupDirectories ("DataAnalysis");

  const bool isOverlayMC = (isMC && strstr (inFileName, "PbPb") != NULL && strstr (inFileName, "Hijing") == NULL);
  if (isOverlayMC)
    cout << "Info: In TreeMaker.cxx: Running over data overlay, will check data conditions" << endl;
    

  const TString identifier = GetIdentifier (dataSet, directory, inFileName);
  cout << "Info: In TreeMaker.cxx: File Identifier: " << identifier << endl;
  cout << "Info: In TreeMaker.cxx: Saving output to " << rootPath << endl;

  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (!isMC) {
      if (dataSet == 0)
        fileIdentifier = "PhysCont.AOD.";
      else
        fileIdentifier = to_string (dataSet);
    }
    else {
      cout << "Error: In TreeMaker.C: Cannot identify this MC file! Quitting." << endl;
      return false;
    }
  }
  else
    fileIdentifier = inFileName;

  TChain* tree = new TChain ("bush", "bush");
  TString pattern = "*.root";
  auto dir = gSystem->OpenDirectory (dataPath + directory);
  while (const char* f = gSystem->GetDirEntry (dir)) {
    TString file = TString (f);
    if (!file.Contains (fileIdentifier))
      continue;
    cout << "Adding " << dataPath + directory + "/" + file + "/*.root" << " to TChain" << endl;
    tree->Add (dataPath + directory + "/" + file + "/*.root");
    break;
  }
  cout << "Chain has " << tree->GetListOfFiles ()->GetEntries () << " files, " << tree->GetEntries () << " entries" << endl;

  if (tree == nullptr) {
    cout << "Error: In TagAndProbe.cxx: TTree not obtained for given data set. Quitting." << endl;
    return false;
  }


  TH1D* h_zdcCuts = nullptr;
  if (isPbPb) {
    h_zdcCuts = GetZdcCuts ();
    if (h_zdcCuts == nullptr) {
      cout << "Error: In TreeMaker.cxx: Zdc in-time pile-up cuts not found. Quitting." << endl;
      return false;
    }
  }

  TFile* f_muonTrigEff = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/TagAndProbe/muontrigger_sf_2017_mc16d_v01.root", "read");
  TH2F* h_muonTrigEff_barrel = (TH2F*) f_muonTrigEff->Get ("Medium/PeriodK/HLT_mu14/eff_etaphi_fine_barrel_data_nominal");
  TH2F* h_muonTrigEff_endcap = (TH2F*) f_muonTrigEff->Get ("Medium/PeriodK/HLT_mu14/eff_etaphi_fine_endcap_data_nominal");

  TFile* f_muonIDEff = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/TagAndProbe/Reco_Medium_Z.root", "read");
  TH3D* h_muonIDEff = (TH3D*) f_muonIDEff->Get ("Eff_2017");

  TFile* f_electronEff = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/TagAndProbe/Nominal/outFile.root", "read");
  TH2D* h_electronTrigEff_pt_eta_pp = (TH2D*) f_electronEff->Get ("h2_electronTrigEffNum_pt_eta_pp")->Clone ("h_electronTrigEff_pt_eta_pp");
  h_electronTrigEff_pt_eta_pp->Divide ((TH2D*) f_electronEff->Get ("h2_electronTrigEffDen_pt_eta_pp"));
  TH2D* h_electronIDEff_pt_eta_pp = (TH2D*) f_electronEff->Get ("h2_electronIDEffNum_pt_eta_pp")->Clone ("h_electronIDEff_pt_eta_pp");
  h_electronIDEff_pt_eta_pp->Divide ((TH2D*) f_electronEff->Get ("h2_electronIDEffDen_pt_eta_pp"));


  TreeVariables* t = new TreeVariables (tree, isMC);
  t->SetGetFCals ();
  t->SetGetVertices ();
  t->SetGetElectrons ();
  t->SetGetMuons ();
  t->SetGetTracks ();
  if (!isPbPb) t->SetGetJets ();
  else if (isPbPb) t->SetGetZdc ();
  t->SetBranchAddresses ();

  Trigger* electronTrigger = nullptr, *muonTrigger = nullptr;
  std::string electronTrigName, muonTrigName;
  if (collisionSystem == PbPb15) {
    electronTrigName = "HLT_e15_loose_ion_L1EM12";
    muonTrigName = "HLT_mu14";
  }
  else if (collisionSystem == pp17) {
    electronTrigName = "HLT_e15_lhloose_L1EM12";
    muonTrigName = "HLT_mu14";
  }
  else if (collisionSystem == PbPb18) {
    electronTrigName = "HLT_e15_lhloose_ion_L1EM12";
    muonTrigName = "HLT_mu14";
  }
  else {
    cout << "Error: In TagAndProbe.cxx: Invalid collision system, quitting." << endl;
    return false;
  }
  if (!isMC) {
    for (int iTrig = 0; iTrig < nElectronTrig; iTrig++) {
      if (electronTrigNames[iTrig] != electronTrigName)
        continue;
      Trigger* trig = new Trigger (electronTrigNames[iTrig].Data (), electronTrigMinPtCuts[iTrig], electronTrigMinEtaCuts[iTrig], electronTrigMaxEtaCuts[iTrig]);
      trig->minPt = electronTrigMinPtCuts[iTrig];
      trig->maxPt = electronTrigMaxPtCuts[iTrig];
      tree->SetBranchAddress (electronTrigNames[iTrig], &(trig->trigBool));
      tree->SetBranchAddress (electronTrigNames[iTrig]+"_prescale", &(trig->trigPrescale));
      electronTrigger = trig;
      break;
    }
    for (int iTrig = 0; iTrig < nMuonTrig; iTrig++) {
      if (muonTrigNames[iTrig] != muonTrigName)
        continue;
      Trigger* trig = new Trigger (muonTrigNames[iTrig].Data (), muonTrigMinPtCuts[iTrig], muonTrigMinEtaCuts[iTrig], muonTrigMaxEtaCuts[iTrig]);
      trig->minPt = muonTrigMinPtCuts[iTrig];
      trig->maxPt = muonTrigMaxPtCuts[iTrig];
      tree->SetBranchAddress (muonTrigNames[iTrig], &(trig->trigBool));
      tree->SetBranchAddress (muonTrigNames[iTrig]+"_prescale", &(trig->trigPrescale));
      muonTrigger = trig;
      break;
    }
  }

  TFile* outFile = new TFile (Form ("%s/%s.root", rootPath.Data (), identifier.Data ()), "recreate");

  const char* outTreeName = isPbPb ? "PbPbZTrackTree" : "ppZTrackTree";
  outFile->Delete (Form ("%s;*", outTreeName));
  OutTree* outTree = new OutTree (outTreeName, outFile);
  outTree->SetBranchEventInfo ();
  outTree->SetBranchLeptons ();
  outTree->SetBranchZs ();
  if (!isPbPb) outTree->SetBranchJets ();
  outTree->SetBranchTracks ();
  outTree->SetBranches ();

  const int nEvts = tree->GetEntries ();
  TLorentzVector l1, l2;

  // First loop over events looking for Z->ee candidates
  isEE = true;
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TreeMaker.cxx: Electrons " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    t->GetEntry (iEvt);

    if (collisionSystem == PbPb18 && t->BlayerDesyn)
      continue;
    if (isPbPb && t->isOOTPU)
      continue;

    event_number = t->event_number;
    run_number = t->run_number;

    bool hasPrimaryVert = false;
    bool hasPileupVert = false;
    for (int iVert = 0; iVert < t->nvert; iVert++) {
      if (t->vert_type[iVert] == 1) {
        if (hasPrimaryVert) {
          hasPrimaryVert = false;
          break;
        }
        hasPrimaryVert = (t->vert_type[iVert] == 1);
        vz = t->vert_z[iVert];
      }
      else if (t->vert_type[iVert] == 3) {
        hasPileupVert = true;
      }
    }
    if (!hasPrimaryVert)// || hasPileupVert)
      continue;

    event_weight = -1;
    if (!isMC) {
      if (electronTrigger->trigBool)
        event_weight = electronTrigger->trigPrescale;
    }
    else {
      event_weight = t->mcEventWeights->at (0) * mcFilterEfficiency / (mcNumberEvents);
    }
    if (event_weight == -1)
      continue; // trigger requirement

    if (isPbPb) {
      fcal_et = t->fcalA_et + t->fcalC_et;
      if (fcal_et < 50)
        continue; // loose cut on "noise" events
      const double qx = t->fcalA_et_Cos + t->fcalC_et_Cos;
      const double qy = t->fcalA_et_Sin + t->fcalC_et_Sin;
      if (fcal_et > 0)
        q2 = sqrt (qx*qx + qy*qy) / fcal_et;
      else
        q2 = 0;
      psi2 = 0.5 * atan2 (qy, qx);

      zdcEnergy = t->ZdcCalibEnergy_A + t->ZdcCalibEnergy_C; // gets zdc energy in TeV
      const float nNeutrons = (zdcEnergy) / (2.51);
      const int bin = h_zdcCuts->FindFixBin (fcal_et * 1e-3); // gets x-axis bin corresponding to Fcal Sum Et in TeV
      if (bin < 1 || h_zdcCuts->GetNbinsX () < bin || (is2015data ? nNeutrons : zdcEnergy) > h_zdcCuts->GetBinContent (bin))
        continue; // Zdc-based in-time pile-up cut
    }

    //const float _event_weight = event_weight; // for storage, in case of multiple Z's in one event

    for (int iE1 = 0; iE1 < t->electron_n; iE1++) {
      if (!doElectronPtDownVar && !doElectronPtUpVar) l1_pt = t->electron_pt[iE1];
      else if (doElectronPtDownVar) l1_pt = t->electron_pt[iE1] - t->electron_pt_sys[iE1];
      else if (doElectronPtUpVar) l1_pt = t->electron_pt[iE1] + t->electron_pt_sys[iE1];
      l1_eta = t->electron_eta[iE1];
      l1_phi = t->electron_phi[iE1];
      l1_charge = t->electron_charge[iE1];
      l1_trk_pt = t->electron_id_track_pt[iE1];
      l1_trk_eta = t->electron_id_track_eta[iE1];
      l1_trk_phi = t->electron_id_track_phi[iE1];
      l1_pt_sys = t->electron_pt_sys[iE1];
      l1_eta_sys = t->electron_eta_sys[iE1];
      l1_phi_sys = t->electron_phi_sys[iE1];

      if (l1_pt < electron_pt_cut)
        continue; // basic electron pT cut
      if (!InEMCal (l1_eta))
        continue; // reject electrons reconstructed outside the EMCal
      //if (t->electron_topoetcone20[iE1] / l1_pt > 0.2)
      //  continue; // isolation cut for electrons

      if (!doElectronLHMediumVar) {
        if (!isPbPb && !t->electron_lhloose[iE1])
          continue; // reject non-loose electrons
        else if (isPbPb && !t->electron_lhloose_hi[iE1])
          continue; // alternative for HI collisions
      }
      else {
        if (!isPbPb && !t->electron_lhmedium[iE1])
          continue; // reject non-medium electrons
        else if (isPbPb && !t->electron_lhmedium_hi[iE1])
          continue; // alternative for HI collisions
      }

      if (fabs (t->electron_id_track_d0sig[iE1]) > 5)
        continue; // electron d0 significance vertex compatibility cut

      if (fabs ( (t->electron_id_track_z0[iE1] + t->electron_id_track_vz[iE1] - vz) * sin (t->electron_id_track_theta[iE1])) > 0.5)
        continue; // electron z0 vertex compatibility cut

      l1.SetPtEtaPhiM (l1_pt, l1_eta, l1_phi, electron_mass);

      for (int iE2 = 0; iE2 < iE1; iE2++) {
        if (!doElectronPtDownVar && !doElectronPtUpVar) l2_pt = t->electron_pt[iE2];
        else if (doElectronPtDownVar) l2_pt = t->electron_pt[iE2] - t->electron_pt_sys[iE2];
        else if (doElectronPtUpVar) l2_pt = t->electron_pt[iE2] + t->electron_pt_sys[iE2];
        l2_eta = t->electron_eta[iE2];
        l2_phi = t->electron_phi[iE2];
        l2_charge = t->electron_charge[iE2];
        l2_trk_pt = t->electron_id_track_pt[iE2];
        l2_trk_eta = t->electron_id_track_eta[iE2];
        l2_trk_phi = t->electron_id_track_phi[iE2];
        l2_pt_sys = t->electron_pt_sys[iE2];
        l2_eta_sys = t->electron_eta_sys[iE2];
        l2_phi_sys = t->electron_phi_sys[iE2];

        if (l2_pt < electron_pt_cut)
          continue; // basic electron pT cut
        if (!InEMCal (l2_eta))
          continue; // reject electrons reconstructed outside the EMCal
        //if (t->electron_topoetcone20[iE2] / l2_pt > 0.2)
        //  continue; // isolation cut for electrons

        if (!doElectronLHMediumVar) {
          if (!isPbPb && !t->electron_lhloose[iE2])
            continue; // reject non-loose electrons
          else if (isPbPb && !t->electron_lhloose_hi[iE2])
            continue; // alternative for HI collisions
        }
        else {
          if (!isPbPb && !t->electron_lhmedium[iE2])
            continue; // reject non-medium electrons
          else if (isPbPb && !t->electron_lhmedium_hi[iE2])
            continue; // alternative for HI collisions
        }

        if (fabs (t->electron_id_track_d0sig[iE2]) > 5)
          continue; // electron d0 significance vertex compatibility cut

        if (fabs ( (t->electron_id_track_z0[iE2] + t->electron_id_track_vz[iE2] - vz) * sin (t->electron_id_track_theta[iE2])) > 0.5)
          continue; // electron z0 vertex compatibility cut

        l2.SetPtEtaPhiM (l2_pt, l2_eta, l2_phi, electron_mass);

        z_pt = (l1+l2).Pt ();
        z_y = (l1+l2).Rapidity ();
        z_phi = (l1+l2).Phi ();
        z_m = (l1+l2).M ();

        if (l1_charge == l2_charge) 
          continue; // require oppositely charged electrons
        if (z_m < 76 || 106 < z_m)
          continue; // require Z to be in mass window

        //const float trigEff1 = (!isMC ? (h_electronTrigEff_pt_eta_pp->GetBinContent (h_electronTrigEff_pt_eta_pp->FindBin (l1_eta, l1_pt))) : 1);
        //const float trigEff2 = (!isMC ? (h_electronTrigEff_pt_eta_pp->GetBinContent (h_electronTrigEff_pt_eta_pp->FindBin (l2_eta, l2_pt))) : 1);
        //const float recoEff1 = (h_electronIDEff_pt_eta_pp->GetBinContent (h_electronIDEff_pt_eta_pp->FindBin (l1_eta, l1_pt)));
        //const float recoEff2 = (h_electronIDEff_pt_eta_pp->GetBinContent (h_electronIDEff_pt_eta_pp->FindBin (l2_eta, l2_pt)));
        //const float combEff = (1.-(1.-trigEff1)*(1.-trigEff2))*recoEff1*recoEff2;

        //if (combEff == 0)
        //  continue;
        //event_weight = _event_weight / combEff;

        ntrk_all = 0;
        ntrk = 0;
        trk_pt.clear ();
        trk_eta.clear ();
        trk_phi.clear ();
        trk_charge.clear ();
        //trk_hiloose.clear ();
        //trk_hitight.clear ();

        for (int iTrk = 0; iTrk < t->ntrk; iTrk++) {
          //if (!t->trk_HItight[iTrk] && !t->trk_HIloose[iTrk])
          //  continue;
          if (doHITightVar && !t->trk_HItight[iTrk])
            continue;
          else if (!doHITightVar && !t->trk_HIloose[iTrk])
            continue;

          ntrk_all++;

          if (t->trk_pt[iTrk] < trk_pt_cut)
            continue; // track minimum pT

          //if (doTrkD0Var && t->trk_d0[iTrk] > 1.0)
          //  continue;
          //if (doTrkZ0Var && fabs (t->trk_z0[iTrk] * sin (t->trk_theta[iTrk])) > 1.0) 
          //  continue;

          if (IsElectronTrack (t, iTrk, iE1, iE2))
            continue;

          trk_pt.push_back (t->trk_pt[iTrk]);
          trk_eta.push_back (t->trk_eta[iTrk]);
          trk_phi.push_back (t->trk_phi[iTrk]);
          trk_charge.push_back (t->trk_charge[iTrk]);
          //trk_hitight.push_back (t->trk_HItight[iTrk]);
          //trk_hiloose.push_back (t->trk_HIloose[iTrk]);
          ntrk++;
        } // end loop over tracks

        if (!isPbPb) {
          njet = 0;
          jet_pt.clear ();
          jet_eta.clear ();
          jet_phi.clear ();
          jet_e.clear ();
          for (int iJet = 0; iJet < t->akt4emtopo_jet_n; iJet++) {
            float _jpt = t->akt4emtopo_jet_pt[iJet];
            float _jeta = t->akt4emtopo_jet_eta[iJet];
            float _jphi = t->akt4emtopo_jet_phi[iJet];
            float _je = t->akt4emtopo_jet_e[iJet];

            if (_jpt < 15)
              continue;
            if (DeltaR (_jeta, l1_eta, _jphi, l1_phi) < 0.2)
              continue;
            if (DeltaR (_jeta, l2_eta, _jphi, l2_phi) < 0.2)
              continue;

            jet_pt.push_back (_jpt);
            jet_eta.push_back (_jeta);
            jet_phi.push_back (_jphi);
            jet_e.push_back (_je);
            njet++;
          }
        }

        outTree->Fill ();
        
      } // end iE2 loop
    } // end iE1 loop

  } // end electron selection
  cout << endl << "Info: In TreeMaker.cxx: Finished processing electrons." << endl;


  // Now loop over events looking for Z->mumu candidates
  isEE = false;
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << "Info: In TreeMaker.cxx: Muons " << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    t->GetEntry (iEvt);

    event_number = t->event_number;
    run_number = t->run_number;

    if (collisionSystem == PbPb18 && t->BlayerDesyn)
      continue;
    if (isPbPb && t->isOOTPU)
      continue;
    if ((!isMC || isOverlayMC) && !t->passes_toroid)
      continue; // additional check for muon quality

    bool hasPrimaryVert = false;
    bool hasPileupVert = false;
    for (int iVert = 0; iVert < t->nvert; iVert++) {
      if (t->vert_type[iVert] == 1) {
        if (hasPrimaryVert) {
          hasPrimaryVert = false;
          break;
        }
        hasPrimaryVert = (t->vert_type[iVert] == 1);
        vz = t->vert_z[iVert];
      }
      else if (t->vert_type[iVert] == 3) {
        hasPileupVert = true;
      }
    }
    if (!hasPrimaryVert)// || hasPileupVert)
      continue; // check for primary vertex

    event_weight = -1;
    if (!isMC) {
      if (muonTrigger->trigBool)
        event_weight = muonTrigger->trigPrescale;
    }
    else {
      event_weight = t->mcEventWeights->at (0) * mcFilterEfficiency / (mcNumberEvents);
    }
    if (event_weight == -1)
      continue; // trigger requirement

    if (isPbPb) {
      fcal_et = t->fcalA_et + t->fcalC_et;
      const double qx = t->fcalA_et_Cos + t->fcalC_et_Cos;
      const double qy = t->fcalA_et_Sin + t->fcalC_et_Sin;
      if (fcal_et > 0)
        q2 = sqrt (qx*qx + qy*qy) / fcal_et;
      else
        q2 = 0;
      psi2 = 0.5 * atan2 (qy, qx);

      zdcEnergy = t->ZdcCalibEnergy_A + t->ZdcCalibEnergy_C; // gets zdc energy in TeV
      const float nNeutrons = (zdcEnergy) / (2.51);
      const int bin = h_zdcCuts->FindFixBin (fcal_et * 1e-3); // gets x-axis bin corresponding to Fcal Sum Et in TeV
      if (bin < 1 || h_zdcCuts->GetNbinsX () < bin || (is2015data ? nNeutrons : zdcEnergy) > h_zdcCuts->GetBinContent (bin))
        continue; // Zdc-based in-time pile-up cut
    }

    //const float _event_weight = event_weight; // for storage, in case of multiple Z's in one event

    for (int iM1 = 0; iM1 < t->muon_n; iM1++) {
      if (!doMuonPtUpVar && !doMuonPtDownVar) l1_pt = t->muon_pt[iM1];
      else if (doMuonPtDownVar) l1_pt = t->muon_pt[iM1] - t->muon_pt_sys[iM1];
      else if (doMuonPtUpVar) l1_pt = t->muon_pt[iM1] + t->muon_pt_sys[iM1];
      l1_eta = t->muon_eta[iM1];
      l1_phi = t->muon_phi[iM1];
      l1_charge = t->muon_charge[iM1];
      l1_trk_pt = t->muon_id_track_pt[iM1];
      l1_trk_eta = t->muon_id_track_eta[iM1];
      l1_trk_phi = t->muon_id_track_phi[iM1];
      l1_pt_sys = t->muon_pt_sys[iM1];
      l1_eta_sys = t->muon_eta_sys[iM1];
      l1_phi_sys = t->muon_phi_sys[iM1];

      if (l1_pt < muon_pt_cut)
        continue; // basic muon pT cut
      if (fabs (l1_eta) > 2.5)
        continue; // reject muons reconstructed outside the MS
      if (!doMuonTightVar && !t->muon_medium[iM1])
        continue; // reject non-tight muons
      else if (doMuonTightVar && !t->muon_tight[iM1])
        continue; // reject non-loose muons

      if (fabs (t->muon_id_track_d0sig[iM1]) > 3)
        continue; // muon d0 cut
      if (fabs ( (t->muon_id_track_z0[iM1] + t->muon_id_track_vz[iM1] - vz) * sin (t->muon_id_track_theta[iM1])) > 0.5)
        continue; // muon z0 to vertex cut

      l1.SetPtEtaPhiM (l1_pt, l1_eta, l1_phi, muon_mass);

      for (int iM2 = 0; iM2 < iM1; iM2++) {
        if (!doMuonPtUpVar && !doMuonPtDownVar) l2_pt = t->muon_pt[iM2];
        else if (doMuonPtDownVar) l2_pt = t->muon_pt[iM2] - t->muon_pt_sys[iM2];
        else if (doMuonPtUpVar) l2_pt = t->muon_pt[iM2] + t->muon_pt_sys[iM2];
        l2_eta = t->muon_eta[iM2];
        l2_phi = t->muon_phi[iM2];
        l2_charge = t->muon_charge[iM2];
        l2_trk_pt = t->muon_id_track_pt[iM2];
        l2_trk_eta = t->muon_id_track_eta[iM2];
        l2_trk_phi = t->muon_id_track_phi[iM2];
        l2_pt_sys = t->muon_pt_sys[iM2];
        l2_eta_sys = t->muon_eta_sys[iM2];
        l2_phi_sys = t->muon_phi_sys[iM2];

        if (l2_pt < muon_pt_cut)
          continue; // basic muon pT cut
        if (fabs (l2_eta) > 2.5)
          continue; // reject muons reconstructed outside the EMCal
        if (!doMuonTightVar && !t->muon_medium[iM2])
          continue; // reject non-tight muons
        else if (doMuonTightVar && !t->muon_tight[iM2])
          continue; // reject non-loose muons

        if (fabs (t->muon_id_track_d0sig[iM2]) > 3)
          continue; // muon d0 significance vertex compatibility cut

        if (fabs ( (t->muon_id_track_z0[iM2] + t->muon_id_track_vz[iM2] - vz) * sin (t->muon_id_track_theta[iM2])) > 0.5)
          continue; // muon z0 vertex compatibility cut

        l2.SetPtEtaPhiM (l2_pt, l2_eta, l2_phi, muon_mass);

        z_pt = (l1+l2).Pt ();
        z_y = (l1+l2).Rapidity ();
        z_phi = (l1+l2).Phi ();
        z_m = (l1+l2).M ();

        if (l1_charge == l2_charge) 
          continue; // require oppositely charged muons
        if (z_m < 76 || 106 < z_m)
          continue; // require Z to be in mass window

        //const float _l1_phi = (InTwoPi (l1_phi) > pi ? InTwoPi (l1_phi) - 2*pi : InTwoPi (l1_phi));
        //const float _l2_phi = (InTwoPi (l2_phi) > pi ? InTwoPi (l2_phi) - 2*pi : InTwoPi (l2_phi));
        //const float trigEff1 = (!isMC ? (fabs (l1_eta) < 1.05 ? h_muonTrigEff_barrel->GetBinContent (h_muonTrigEff_barrel->FindBin (l1_eta, _l1_phi)) : h_muonTrigEff_endcap->GetBinContent (h_muonTrigEff_endcap->FindBin (l1_eta, _l1_phi))) : 1);
        //const float trigEff2 = (!isMC ? (fabs (l2_eta) < 1.05 ? h_muonTrigEff_barrel->GetBinContent (h_muonTrigEff_barrel->FindBin (l2_eta, _l2_phi)) : h_muonTrigEff_endcap->GetBinContent (h_muonTrigEff_endcap->FindBin (l2_eta, _l2_phi))) : 1);
        //const float recoEff1 = h_muonIDEff->GetBinContent (h_muonIDEff->FindBin (l1_eta, _l1_phi, l1_pt));
        //const float recoEff2 = h_muonIDEff->GetBinContent (h_muonIDEff->FindBin (l2_eta, _l2_phi, l2_pt));
        //const float combEff = (1-(1-trigEff1)*(1-trigEff2))*recoEff1*recoEff2;

        //if (combEff == 0)
        //  continue;
        //event_weight = _event_weight / combEff;

        ntrk_all = 0;
        ntrk = 0;
        trk_pt.clear ();
        trk_eta.clear ();
        trk_phi.clear ();
        trk_charge.clear ();
        //trk_hiloose.clear ();
        //trk_hitight.clear ();

        for (int iTrk = 0; iTrk < t->ntrk; iTrk++) {
          //if (!t->trk_HItight[iTrk] && !t->trk_HIloose[iTrk])
          //  continue;
          if (doHITightVar && !t->trk_HItight[iTrk])
            continue;
          else if (!doHITightVar && !t->trk_HIloose[iTrk])
            continue;

          ntrk_all++;

          if (t->trk_pt[iTrk] < trk_pt_cut)
            continue; // track minimum pT

          //if (doTrkD0Var && t->trk_d0[iTrk] > 1.0)
          //  continue;
          //if (doTrkZ0Var && fabs (t->trk_z0[iTrk] * sin (t->trk_theta[iTrk])) > 1.0) 
          //  continue;

          if (IsMuonTrack (t, iTrk, iM1, iM2))
            continue;

          trk_pt.push_back (t->trk_pt[iTrk]);
          trk_eta.push_back (t->trk_eta[iTrk]);
          trk_phi.push_back (t->trk_phi[iTrk]);
          trk_charge.push_back (t->trk_charge[iTrk]);
          //trk_hitight.push_back (t->trk_HItight[iTrk]);
          //trk_hiloose.push_back (t->trk_HIloose[iTrk]);
          ntrk++;
        } // end loop over tracks

        if (!isPbPb) {
          njet = 0;
          jet_pt.clear ();
          jet_eta.clear ();
          jet_phi.clear ();
          jet_e.clear ();
          for (int iJet = 0; iJet < t->akt4emtopo_jet_n; iJet++) {
            float _jpt = t->akt4emtopo_jet_pt[iJet];
            float _jeta = t->akt4emtopo_jet_eta[iJet];
            float _jphi = t->akt4emtopo_jet_phi[iJet];
            float _je = t->akt4emtopo_jet_e[iJet];

            if (_jpt < 15)
              continue;
            if (DeltaR (_jeta, l1_eta, _jphi, l1_phi) < 0.2)
              continue;
            if (DeltaR (_jeta, l2_eta, _jphi, l2_phi) < 0.2)
              continue;

            jet_pt.push_back (_jpt);
            jet_eta.push_back (_jeta);
            jet_phi.push_back (_jphi);
            jet_e.push_back (_je);
            njet++;
          }
        }

        outTree->Fill ();
        
      } // end iM2 loop
    } // end iM1 loop

  } // end muon selection
  cout << endl << "Info: In TreeMaker.cxx: Finished processing muons." << endl;

  SaferDelete (electronTrigger);
  SaferDelete (muonTrigger);

  SaferDelete (t);

  outFile->Write (0, TObject::kOverwrite);
  outFile->Close ();
  SaferDelete (outFile);

  SaferDelete (h_zdcCuts);

  return true;
}

} // end namespace
