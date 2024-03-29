#ifndef __FullAnalysis_h__
#define __FullAnalysis_h__

#include "Params.h"
#include "PhysicsAnalysis.h"

#include <ArrayTemplates.h>

#include <TEfficiency.h>
#include <TClass.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TLine.h>
#include <TF1.h>
#include <TLorentzVector.h>

#include <iostream>
#include <string>

using namespace std;

class FullAnalysis : public PhysicsAnalysis {

  public:

  // A full analysis includes checks of Z boson distributions
  TH1D*** h_z_phi           = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH1D*** h_z_pt            = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH1D*** h_z_pt_ratio      = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH2D**** h_z_y_phi        = Get3DArray <TH2D*> (numCentBins, 3, nPtZBins+1);   // iCent, iSpc, iPtZ + integrated bin
  TH1D**** h_z_eta          = Get3DArray <TH1D*> (numCentBins, 3, nPtZBins+1);   // iCent, iSpc, iPtZ + integrated bin
  TH1D**** h_z_eta_ratio    = Get3DArray <TH1D*> (numCentBins, 3, nPtZBins+1);   // iCent, iSpc, iPtZ + integrated bin
  TH1D**** h_z_y            = Get3DArray <TH1D*> (numCentBins, 3, nPtZBins+1);   // iCent, iSpc, iPtZ + integrated bin
  TH1D**** h_z_y_ratio      = Get3DArray <TH1D*> (numCentBins, 3, nPtZBins+1);   // iCent, iSpc, iPtZ + integrated bin
  TH1D**** h_z_m            = Get3DArray <TH1D*> (numCentBins, 3, 3);          // iCent, iSpc, iReg
  TH1D**** h_z_m_ratio      = Get3DArray <TH1D*> (numCentBins, 3, 3);          // iCent, iSpc, iReg
  TH1D*** h_z_lepton_dphi   = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH1D*** h_lepton_pt       = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH1D*** h_lepton_pt_ratio = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH1D*** h_lepton_eta      = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH1D*** h_lepton_trk_pt   = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH1D*** h_trk_pt          = Get2DArray <TH1D*> (numCentBins, 3);             // iCent, iSpc
  TH2D*** h_lepton_trk_dr   = Get2DArray <TH2D*> (numCentBins, 3);             // iCent, iSpc

  TGraph**** g_trk_pt_ptz             = Get3DArray <TGraph*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TGraph**** g_trk_xhz_ptz            = Get3DArray <TGraph*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent


  FullAnalysis () { }

  FullAnalysis (const char* _name) {//, const char* subDir) {
    name = _name;
    //directory = Form ("%s/", subDir);
    plotFill = false;
  }

  virtual ~FullAnalysis () {

    Delete2DArray (&h_z_phi,         numCentBins, 3);
    Delete2DArray (&h_z_pt,          numCentBins, 3);
    Delete2DArray (&h_z_pt_ratio,    numCentBins, 3);
    Delete3DArray (&h_z_y_phi,       numCentBins, 3, nPtZBins+1);
    Delete3DArray (&h_z_eta,         numCentBins, 3, nPtZBins+1);
    Delete3DArray (&h_z_eta_ratio,   numCentBins, 3, nPtZBins+1);
    Delete3DArray (&h_z_y,           numCentBins, 3, nPtZBins+1);
    Delete3DArray (&h_z_y_ratio,     numCentBins, 3, nPtZBins+1);
    Delete3DArray (&h_z_m,           numCentBins, 3, 3);
    Delete3DArray (&h_z_m_ratio,     numCentBins, 3, 3);
    Delete2DArray (&h_z_lepton_dphi, numCentBins, 3);
    Delete2DArray (&h_lepton_pt,     numCentBins, 3);
    Delete2DArray (&h_lepton_eta,    numCentBins, 3);
    Delete2DArray (&h_lepton_trk_pt, numCentBins, 3);
    Delete2DArray (&h_trk_pt,        numCentBins, 3);
    Delete2DArray (&h_lepton_trk_dr, numCentBins, 3);

    Delete3DArray (&g_trk_pt_ptz,    3, nPtZBins, numCentBins);
    Delete3DArray (&g_trk_xhz_ptz,   3, nPtZBins, numCentBins);
  }


  protected:
  void LabelZMassSpectra (const short iSpc, const short iCent, const short iReg);

  public:

  virtual void CreateHists () override;
  virtual void CopyAnalysis (FullAnalysis* a, const bool copyBkgs = false);
  virtual void ClearHists () override;
  virtual void CombineFullHists ();
  virtual void LoadHists (const char* histFileName = "savedHists.root", const bool _finishHists = true) override;
  virtual void SaveHists (const char* histFileName = "savedHists.root") override;
  virtual void ScaleHists () override;
  virtual void Execute (const char* inFileName, const char* outFileName) override;

  //void PrintZYields ();

  void PlotLeptonPtSpectra (FullAnalysis* a = nullptr);
  void PlotLeptonEtaSpcComp (FullAnalysis* a = nullptr);
  void PlotLeptonTrackPtSpectra ();
  void PlotLeptonTrackDR ();
  void PlotLeptonTrackDRProjX ();
  void PlotZPtSpectra (FullAnalysis* a = nullptr);
  void PlotZPtSpectraChannelCompare ();
  void CalculateZPtMeans ();
  void PlotZYPhiMap ();
  void PlotZEtaMap (FullAnalysis* a = nullptr);
  void PlotZYMap (FullAnalysis* a = nullptr);
  void PlotZYMapSpcComp (const short pPtZ, FullAnalysis* a = nullptr);
  void PlotZMassSpectra (FullAnalysis* a = nullptr);
  void PlotZPhiYield (const short pSpc = 2);
  void PlotZLeptonDPhi ();
  void PlotAllYields_Scatter_dPtZ (const bool useTrkPt = true, const short pSpc = 2);

};


////////////////////////////////////////////////////////////////////////////////////////////////
// Create new histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: CreateHists () {
  PhysicsAnalysis :: CreateHists ();

  const int n_z_pt_bins = 300;
  const double* z_pt_bins = logspace (1, 300, n_z_pt_bins);
  //const double z_pt_bins[401] = {1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5, 4.55, 4.6, 4.65, 4.7, 4.75, 4.8, 4.85, 4.9, 4.95, 5, 5.125, 5.25, 5.375, 5.5, 5.625, 5.75, 5.875, 6, 6.125, 6.25, 6.375, 6.5, 6.625, 6.75, 6.875, 7, 7.125, 7.25, 7.375, 7.5, 7.625, 7.75, 7.875, 8, 8.125, 8.25, 8.375, 8.5, 8.625, 8.75, 8.875, 9, 9.125, 9.25, 9.375, 9.5, 9.625, 9.75, 9.875, 10, 10.125, 10.25, 10.375, 10.5, 10.625, 10.75, 10.875, 11, 11.125, 11.25, 11.375, 11.5, 11.625, 11.75, 11.875, 12, 12.125, 12.25, 12.375, 12.5, 12.625, 12.75, 12.875, 13, 13.125, 13.25, 13.375, 13.5, 13.625, 13.75, 13.875, 14, 14.125, 14.25, 14.375, 14.5, 14.625, 14.75, 14.875, 15, 15.1875, 15.375, 15.5625, 15.75, 15.9375, 16.125, 16.3125, 16.5, 16.6875, 16.875, 17.0625, 17.25, 17.4375, 17.625, 17.8125, 18, 18.1875, 18.375, 18.5625, 18.75, 18.9375, 19.125, 19.3125, 19.5, 19.6875, 19.875, 20.0625, 20.25, 20.4375, 20.625, 20.8125, 21, 21.1875, 21.375, 21.5625, 21.75, 21.9375, 22.125, 22.3125, 22.5, 22.6875, 22.875, 23.0625, 23.25, 23.4375, 23.625, 23.8125, 24, 24.1875, 24.375, 24.5625, 24.75, 24.9375, 25.125, 25.3125, 25.5, 25.6875, 25.875, 26.0625, 26.25, 26.4375, 26.625, 26.8125, 27, 27.1875, 27.375, 27.5625, 27.75, 27.9375, 28.125, 28.3125, 28.5, 28.6875, 28.875, 29.0625, 29.25, 29.4375, 29.625, 29.8125, 30, 30.375, 30.75, 31.125, 31.5, 31.875, 32.25, 32.625, 33, 33.375, 33.75, 34.125, 34.5, 34.875, 35.25, 35.625, 36, 36.375, 36.75, 37.125, 37.5, 37.875, 38.25, 38.625, 39, 39.375, 39.75, 40.125, 40.5, 40.875, 41.25, 41.625, 42, 42.375, 42.75, 43.125, 43.5, 43.875, 44.25, 44.625, 45, 45.375, 45.75, 46.125, 46.5, 46.875, 47.25, 47.625, 48, 48.375, 48.75, 49.125, 49.5, 49.875, 50.25, 50.625, 51, 51.375, 51.75, 52.125, 52.5, 52.875, 53.25, 53.625, 54, 54.375, 54.75, 55.125, 55.5, 55.875, 56.25, 56.625, 57, 57.375, 57.75, 58.125, 58.5, 58.875, 59.25, 59.625, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 96, 99, 102, 105, 108, 111, 114, 117, 120, 123, 126, 129, 132, 135, 138, 141, 144, 147, 150, 153, 156, 159, 162, 165, 168, 171, 174, 177, 180, 183, 186, 189, 192, 195, 198, 201, 204, 207, 210, 213, 216, 219, 222, 225, 228, 231, 234, 237, 240, 243, 246, 249, 252, 255, 258, 261, 264, 267, 270, 273, 276, 279, 282, 285, 288, 291, 294, 297, 300};
  const double* lepton_pt_bins = logspace (20, 250, 230);

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      h_z_phi[iCent][iSpc]          = new TH1D (Form ("h_z_phi_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 80, 0, pi);
      h_z_pt[iCent][iSpc]           = new TH1D (Form ("h_z_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", n_z_pt_bins, z_pt_bins);
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_y_phi[iCent][iSpc][iPtZ]  = new TH2D (Form ("h_z_y_phi_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()), "", 20, -2.5, 2.5, 20, 0, 2*pi);
        h_z_eta[iCent][iSpc][iPtZ]    = new TH1D (Form ("h_z_eta_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()), "", 16, -5.0, 5.0);
        h_z_y[iCent][iSpc][iPtZ]      = new TH1D (Form ("h_z_y_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()), "", 16, -2.5, 2.5);

        g_trk_pt_ptz[iSpc][iPtZ][iCent] = new TGraph ();
        g_trk_pt_ptz[iSpc][iPtZ][iCent]->SetName (Form ("g_trk_pt_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        g_trk_xhz_ptz[iSpc][iPtZ][iCent] = new TGraph ();
        g_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetName (Form ("g_trk_xhz_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
      } // end loop over iPtZ
      for (short iReg = 0; iReg < 3; iReg++) {
        h_z_m[iCent][iSpc][iReg]    = new TH1D (Form ("h_z_m_%s_iCent%i_iReg%i_%s", spc, iCent, iReg, name.c_str ()), "", 40, 76, 106);
        h_z_m[iCent][iSpc][iReg]->Sumw2 ();
      } // end loop over iReg
      h_z_lepton_dphi[iCent][iSpc]  = new TH1D (Form ("h_z_lepton_dphi_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 45, 0, pi);
      h_lepton_pt[iCent][iSpc]      = new TH1D (Form ("h_lepton_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 230, lepton_pt_bins);
      h_lepton_eta[iCent][iSpc]     = new TH1D (Form ("h_lepton_eta_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 50, -2.5, 2.5);
      h_lepton_trk_pt[iCent][iSpc]  = new TH1D (Form ("h_lepton_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 230, lepton_pt_bins);
      h_trk_pt[iCent][iSpc]         = new TH1D (Form ("h_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 100, 0, 100);
      h_lepton_trk_dr[iCent][iSpc]  = new TH2D (Form ("h_lepton_trk_dr_%s_iCent%i_%s", spc, iCent, name.c_str ()), "", 100, 0, 1, 40, 0, 2);
      
      h_z_phi[iCent][iSpc]->Sumw2 ();
      h_z_pt[iCent][iSpc]->Sumw2 ();
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_y_phi[iCent][iSpc][iPtZ]->Sumw2 ();
        h_z_eta[iCent][iSpc][iPtZ]->Sumw2 ();
        h_z_y[iCent][iSpc][iPtZ]->Sumw2 ();
      }
      h_z_lepton_dphi[iCent][iSpc]->Sumw2 ();
      h_lepton_pt[iCent][iSpc]->Sumw2 ();
      h_lepton_eta[iCent][iSpc]->Sumw2 ();
      h_lepton_trk_pt[iCent][iSpc]->Sumw2 ();
      h_trk_pt[iCent][iSpc]->Sumw2 ();
      h_lepton_trk_dr[iCent][iSpc]->Sumw2 ();
    } // end loop over iSpc
  } // end loop over iCent

  //delete[] z_pt_bins;
  delete[] lepton_pt_bins;

  histsLoaded = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Create new histograms as clones from another analysis, where appropriate
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: CopyAnalysis (FullAnalysis* a, const bool copyBkgs) {
  if (name == "")
    cout << "Warning in FullAnalysis :: CopyAnalysis: name of analysis not set!" << endl;

  PhysicsAnalysis :: CopyAnalysis ((PhysicsAnalysis*)a, copyBkgs);

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      h_z_phi[iCent][iSpc]          = (TH1D*) a->h_z_phi[iCent][iSpc]->Clone (Form ("h_z_phi_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_z_pt[iCent][iSpc]           = (TH1D*) a->h_z_pt[iCent][iSpc]->Clone (Form ("h_z_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      if (a->h_z_pt_ratio[iCent][iSpc])
        h_z_pt_ratio[iCent][iSpc]   = (TH1D*) a->h_z_pt_ratio[iCent][iSpc]->Clone (Form ("h_z_pt_ratio_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_y_phi[iCent][iSpc][iPtZ]        = (TH2D*) a->h_z_y_phi[iCent][iSpc][iPtZ]->Clone (Form ("h_z_y_phi_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        h_z_eta[iCent][iSpc][iPtZ]          = (TH1D*) a->h_z_eta[iCent][iSpc][iPtZ]->Clone (Form ("h_z_eta_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        if (a->h_z_eta_ratio[iCent][iSpc][iPtZ])
          h_z_eta_ratio[iCent][iSpc][iPtZ]  = (TH1D*) a->h_z_eta_ratio[iCent][iSpc][iPtZ]->Clone (Form ("h_z_eta_ratio_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        h_z_y[iCent][iSpc][iPtZ]            = (TH1D*) a->h_z_y[iCent][iSpc][iPtZ]->Clone (Form ("h_z_y_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        if (a->h_z_y_ratio[iCent][iSpc][iPtZ])
          h_z_y_ratio[iCent][iSpc][iPtZ]    = (TH1D*) a->h_z_y_ratio[iCent][iSpc][iPtZ]->Clone (Form ("h_z_y_ratio_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        g_trk_pt_ptz[iSpc][iPtZ][iCent]     = (TGraph*) a->g_trk_pt_ptz[iSpc][iPtZ][iCent]->Clone (Form ("g_trk_pt_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt
        g_trk_xhz_ptz[iSpc][iPtZ][iCent]    = (TGraph*) a->g_trk_xhz_ptz[iSpc][iPtZ][iCent]->Clone (Form ("g_trk_xhz_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh
      } // end loop over iPtZ

      for (short iReg = 0; iReg < 3; iReg++) {
        h_z_m[iCent][iSpc][iReg]    = (TH1D*) a->h_z_m[iCent][iSpc][iReg]->Clone (Form ("h_z_m_%s_iCent%i_iReg%i_%s", spc, iCent, iReg, name.c_str ()));
        if (a->h_z_m_ratio[iCent][iSpc][iReg])
          h_z_m_ratio[iCent][iSpc][iReg] = (TH1D*) a->h_z_m_ratio[iCent][iSpc][iReg]->Clone (Form ("h_z_m_ratio_%s_iCent%i_iReg%i_%s", spc, iCent, iReg, name.c_str ()));
      } // end loop over iReg

      if (a->h_z_lepton_dphi[iCent][iSpc])  h_z_lepton_dphi[iCent][iSpc]  = (TH1D*) a->h_z_lepton_dphi[iCent][iSpc]->Clone (Form ("h_z_lepton_dphi_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      if (a->h_lepton_pt[iCent][iSpc])      h_lepton_pt[iCent][iSpc]      = (TH1D*) a->h_lepton_pt[iCent][iSpc]->Clone (Form ("h_lepton_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      if (a->h_lepton_eta[iCent][iSpc])     h_lepton_eta[iCent][iSpc]     = (TH1D*) a->h_lepton_eta[iCent][iSpc]->Clone (Form ("h_lepton_eta_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      if (a->h_lepton_trk_pt[iCent][iSpc])  h_lepton_trk_pt[iCent][iSpc]  = (TH1D*) a->h_lepton_trk_pt[iCent][iSpc]->Clone (Form ("h_lepton_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      if (a->h_trk_pt[iCent][iSpc])         h_trk_pt[iCent][iSpc]         = (TH1D*) a->h_trk_pt[iCent][iSpc]->Clone (Form ("h_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      if (a->h_lepton_trk_dr[iCent][iSpc])  h_lepton_trk_dr[iCent][iSpc]  = (TH2D*) a->h_lepton_trk_dr[iCent][iSpc]->Clone (Form ("h_lepton_trk_dr_%s_iCent%i_%s", spc, iCent, name.c_str ()));

    } // end loop over iSpc
  } // end loop over iCent

  histsLoaded = true;
  histsScaled = true;
  return;    
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Clears histograms from memory (allows them to be overwritten).
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: ClearHists () {

  PhysicsAnalysis :: ClearHists ();

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        SaferDelete (&g_trk_pt_ptz[iSpc][iPtZ][iCent]);
        SaferDelete (&g_trk_xhz_ptz[iSpc][iPtZ][iCent]);
      } // end loop over iPtZ
    } // end loop over iSpc
  } // end loop over iCent
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Load pre-filled histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: LoadHists (const char* histFileName, const bool _finishHists) {
  if (histsLoaded)
    return;

  PhysicsAnalysis :: LoadHists (histFileName, false);

  TDirectory* _gDirectory = gDirectory;
  if (!histFile->IsOpen ()) {
    cout << "Error in FullAnalysis :: LoadHists: histFile not open after calling parent function, exiting." << endl;
    return;
  }

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      h_z_phi[iCent][iSpc]          = (TH1D*) histFile->Get (Form ("h_z_phi_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_z_pt[iCent][iSpc]           = (TH1D*) histFile->Get (Form ("h_z_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_z_y_phi[iCent][iSpc][nPtZBins]    = new TH2D (Form ("h_z_y_phi_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()), "", 20, -2.5, 2.5, 20, 0, 2*pi);
        h_z_y_phi[iCent][iSpc][iPtZ]        = (TH2D*) histFile->Get (Form ("h_z_y_phi_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        h_z_eta[iCent][iSpc][nPtZBins]      = new TH1D (Form ("h_z_eta_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()), "", 16, -5.0, 5.0);
        h_z_eta[iCent][iSpc][iPtZ]          = (TH1D*) histFile->Get (Form ("h_z_eta_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        h_z_y[iCent][iSpc][nPtZBins]        = new TH1D (Form ("h_z_y_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()), "", 16, -2.5, 2.5);
        h_z_y[iCent][iSpc][iPtZ]            = (TH1D*) histFile->Get (Form ("h_z_y_%s_iCent%i_iPtZ%i_%s", spc, iCent, iPtZ, name.c_str ()));
        g_trk_pt_ptz[iSpc][iPtZ][iCent]       = (TGraph*) histFile->Get (Form ("g_trk_pt_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zpt
        g_trk_xhz_ptz[iSpc][iPtZ][iCent]      = (TGraph*) histFile->Get (Form ("g_trk_xhz_ptz_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ())); // old name: h_z_trk_zxzh
      } // end loop over iPtZ

      for (short iReg = 0; iReg < 3; iReg++) {
        h_z_m[iCent][iSpc][iReg]        = (TH1D*) histFile->Get (Form ("h_z_m_%s_iCent%i_iReg%i_%s", spc, iCent, iReg, name.c_str ()));
      } // end loop over iReg
      h_z_lepton_dphi[iCent][iSpc]  = (TH1D*) histFile->Get (Form ("h_z_lepton_dphi_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_lepton_pt[iCent][iSpc]      = (TH1D*) histFile->Get (Form ("h_lepton_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_lepton_eta[iCent][iSpc]     = (TH1D*) histFile->Get (Form ("h_lepton_eta_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_lepton_trk_pt[iCent][iSpc]  = (TH1D*) histFile->Get (Form ("h_lepton_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_trk_pt[iCent][iSpc]         = (TH1D*) histFile->Get (Form ("h_trk_pt_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h_lepton_trk_dr[iCent][iSpc]  = (TH2D*) histFile->Get (Form ("h_lepton_trk_dr_%s_iCent%i_%s", spc, iCent, name.c_str ()));

    } // end loop over iSpc
  } // end loop over iCent
  
  histsLoaded = true;

  if (_finishHists) {
    CombineFullHists ();
    ScaleHists ();
    SetVariances ();
    MakeSummaryCorrelations ();
  }

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Save histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: SaveHists (const char* histFileName) {
  if (!histsLoaded)
    return;

  PhysicsAnalysis :: SaveHists (histFileName);

  TDirectory* _gDirectory = gDirectory;
  if (!histFile) {
    histFile = new TFile (Form ("%s/%s", rootPath.Data (), histFileName), "update");
    histFile->cd ();
  }
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      SafeWrite (h_z_phi[iCent][iSpc]);
      SafeWrite (h_z_pt[iCent][iSpc]);
      SafeWrite (h_z_pt_ratio[iCent][iSpc]);
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        SafeWrite (h_z_y_phi[iCent][iSpc][iPtZ]);
        SafeWrite (h_z_eta[iCent][iSpc][iPtZ]);
        SafeWrite (h_z_eta_ratio[iCent][iSpc][iPtZ]);
        SafeWrite (h_z_y[iCent][iSpc][iPtZ]);
        SafeWrite (h_z_y_ratio[iCent][iSpc][iPtZ]);
        SafeWrite (g_trk_pt_ptz[iSpc][iPtZ][iCent]);
        SafeWrite (g_trk_xhz_ptz[iSpc][iPtZ][iCent]);
      } // end loop over iPtZ
      for (short iReg = 0; iReg < 3; iReg++) {
        SafeWrite (h_z_m[iCent][iSpc][iReg]);
        SafeWrite (h_z_m_ratio[iCent][iSpc][iReg]);
      } // end loop over iReg
      SafeWrite (h_z_lepton_dphi[iCent][iSpc]);
      SafeWrite (h_lepton_pt[iCent][iSpc]);
      SafeWrite (h_lepton_eta[iCent][iSpc]);
      SafeWrite (h_lepton_trk_pt[iCent][iSpc]);
      SafeWrite (h_trk_pt[iCent][iSpc]);
      SafeWrite (h_lepton_trk_dr[iCent][iSpc]);
    } // end loop over iSpc
  } // end loop over iCent
  
  histFile->Close ();
  histFile = nullptr;
  histsLoaded = false;

  _gDirectory->cd ();
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Fill combined species histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: CombineFullHists () {

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 2; iSpc++) {
      if (h_z_phi[iCent][iSpc]) h_z_phi[iCent][2]->Add (h_z_phi[iCent][iSpc]);
      if (h_z_pt[iCent][iSpc]) h_z_pt[iCent][2]->Add (h_z_pt[iCent][iSpc]);
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        if (h_z_y_phi[iCent][iSpc][iPtZ]) {
          h_z_y_phi[iCent][2][iPtZ]->Add (h_z_y_phi[iCent][iSpc][iPtZ]);
          if (iPtZ >= 2) {
            h_z_y_phi[iCent][iSpc][nPtZBins]->Add (h_z_y_phi[iCent][iSpc][iPtZ]);
            h_z_y_phi[iCent][2][nPtZBins]->Add (h_z_y_phi[iCent][iSpc][iPtZ]);
          }
        }
        if (h_z_eta[iCent][iSpc][iPtZ]) {
          h_z_eta[iCent][2][iPtZ]->Add (h_z_eta[iCent][iSpc][iPtZ]);
          if (iPtZ >= 2) {
            h_z_eta[iCent][iSpc][nPtZBins]->Add (h_z_eta[iCent][iSpc][iPtZ]);
            h_z_eta[iCent][2][nPtZBins]->Add (h_z_eta[iCent][iSpc][iPtZ]);
          }
        }
        if (h_z_y[iCent][iSpc][iPtZ]) {
          h_z_y[iCent][2][iPtZ]->Add (h_z_y[iCent][iSpc][iPtZ]);
          if (iPtZ >= 2) {
            h_z_y[iCent][iSpc][nPtZBins]->Add (h_z_y[iCent][iSpc][iPtZ]);
            h_z_y[iCent][2][nPtZBins]->Add (h_z_y[iCent][iSpc][iPtZ]);
          }
        }

        double x, y;
        for (int i = 0; i < g_trk_pt_ptz[iSpc][iPtZ][iCent]->GetN (); i++) {
          g_trk_pt_ptz[iSpc][iPtZ][iCent]->GetPoint (i, x, y);
          g_trk_pt_ptz[2][iPtZ][iCent]->SetPoint (g_trk_pt_ptz[2][iPtZ][iCent]->GetN (), x, y);
        } // end loop over i
        for (int i = 0; i < g_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetN (); i++) {
          g_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetPoint (i, x, y);
          g_trk_xhz_ptz[2][iPtZ][iCent]->SetPoint (g_trk_xhz_ptz[2][iPtZ][iCent]->GetN (), x, y);
        } // end loop over i
      } // end loop over iPtZ
      for (short iReg = 0; iReg < 2; iReg++) {
        if (h_z_m[iCent][iSpc][iReg]) h_z_m[iCent][2][iReg]->Add (h_z_m[iCent][iSpc][iReg]);
        if (h_z_m[iCent][iSpc][iReg]) h_z_m[iCent][iSpc][2]->Add (h_z_m[iCent][iSpc][iReg]);
        if (h_z_m[iCent][iSpc][iReg]) h_z_m[iCent][2][2]->Add (h_z_m[iCent][iSpc][iReg]);
      } // end loop over iReg
      if (h_lepton_pt[iCent][iSpc]) h_lepton_pt[iCent][2]->Add (h_lepton_pt[iCent][iSpc]);
      if (h_lepton_eta[iCent][iSpc]) h_lepton_eta[iCent][2]->Add (h_lepton_eta[iCent][iSpc]);
      if (h_lepton_pt[iCent][iSpc]) h_lepton_trk_pt[iCent][2]->Add (h_lepton_pt[iCent][iSpc]);
      if (h_trk_pt[iCent][iSpc]) h_trk_pt[iCent][2]->Add (h_trk_pt[iCent][iSpc]);
      if (h_lepton_trk_dr[iCent][iSpc]) h_lepton_trk_dr[iCent][2]->Add (h_lepton_trk_dr[iCent][iSpc]);
    } // end loop over iSpc
  } // end loop over iCent

  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Scale histograms for plotting, calculating signals, etc.
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: ScaleHists () {
  if (histsScaled || !histsLoaded)
    return;

  PhysicsAnalysis :: ScaleHists ();

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {

      double normFactor = 0;
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++)
        normFactor += h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (2);

      if (h_lepton_pt[iCent][iSpc]) {
        h_lepton_pt[iCent][iSpc]->Rebin (10);
        if (normFactor > 0)
          h_lepton_pt[iCent][iSpc]->Scale (1 / normFactor, "width");
      }

      if (h_lepton_eta[iCent][iSpc]) {
        h_lepton_eta[iCent][iSpc]->Rebin (2);
        if (normFactor > 0)
          h_lepton_eta[iCent][iSpc]->Scale (1 / normFactor, "width");
      }

      if (h_lepton_trk_pt[iCent][iSpc]) {
        h_lepton_trk_pt[iCent][iSpc]->Rebin (5);
        if (normFactor > 0)
          h_lepton_trk_pt[iCent][iSpc]->Scale (1 / normFactor, "width");
      }

      if (h_trk_pt[iCent][iSpc]) {
        h_trk_pt[iCent][iSpc]->Rebin (5);
        if (normFactor > 0)
          h_trk_pt[iCent][iSpc]->Scale (1 / normFactor, "width");
      }

      if (h_z_pt[iCent][iSpc]) {
        h_z_pt[iCent][iSpc]->Rebin (3);
        //h_z_pt[iCent][iSpc]->Scale (1./10.);
        if (h_z_pt[iCent][iSpc]->Integral () > 0)
          h_z_pt[iCent][iSpc]->Scale (1 / h_z_pt[iCent][iSpc]->Integral (), "width");
      }

      for (short iReg = 0; iReg < 3; iReg++) {
        //if (h_z_m[iCent][iSpc]->GetMaximum () > 0)
        //  h_z_m[iCent][iSpc]->Scale (1 / (h_z_m[iCent][iSpc]->GetMaximum ()));
        if (h_z_m[iCent][iSpc][iReg]) {
          if (h_z_m[iCent][iSpc][iReg]->Integral () > 0)
            h_z_m[iCent][iSpc][iReg]->Scale (1 / h_z_m[iCent][iSpc][iReg]->Integral ());
        }
      } // end loop over iReg

      for (short iPtZ = 0; iPtZ <= nPtZBins; iPtZ++) {
        if (h_z_eta[iCent][iSpc][iPtZ]) {
          if (h_z_eta[iCent][iSpc][iPtZ]->Integral () > 0)
            h_z_eta[iCent][iSpc][iPtZ]->Scale (1 / h_z_eta[iCent][iSpc][iPtZ]->Integral (), "width");
        }

        if (h_z_y[iCent][iSpc][iPtZ]) {
          //h_z_y[iCent][iSpc][iPtZ]->Rebin (2);
          if (h_z_y[iCent][iSpc][iPtZ]->Integral () > 0)
            h_z_y[iCent][iSpc][iPtZ]->Scale (1 / h_z_y[iCent][iSpc][iPtZ]->Integral (), "width");
        }
      } // end loop over iPtZ

      if (h_z_phi[iCent][iSpc]) {
        h_z_phi[iCent][iSpc]->Rebin (8);
        if (h_z_phi[iCent][iSpc]->Integral () > 0)
          h_z_phi[iCent][iSpc]->Scale (1 / h_z_phi[iCent][iSpc]->Integral (), "width");
      }
    } // end loop over iCent
  } // end loop over iSpc

  histsScaled = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
// Designed to be overloaded. The default here is for analyzing data.
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: Execute (const char* inFileName, const char* outFileName) {

  cout << "Arguments provided: " << endl;
  cout << "inFileName = " << inFileName << endl;
  cout << "outFileName = " << outFileName << endl;

  LoadEventWeights ();
  //eventPlaneCalibrator = EventPlaneCalibrator (Form ("%s/FCalCalibration/Nominal/data18hi.root", rootPath.Data ()));

  TFile* inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), inFileName) << endl;

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  CreateHists ();

  bool isEE = false;
  float event_weight = 1, fcal_weight = 1, q2_weight = 1, psi2_weight = 1;
  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0;
  float z_pt = 0, z_eta = 0, z_y = 0, z_phi = 0, z_m = 0;
  float l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0;
  float l1_trk_pt = 0, l1_trk_eta = 0, l1_trk_phi = 0, l2_trk_pt = 0, l2_trk_eta = 0, l2_trk_phi = 0;
  int l1_charge = 0, l2_charge = 0, ntrk = 0;
  float trk_pt[10000], trk_eta[10000], trk_phi[10000];

  int***    trks_counts   = Get3DArray <int> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  float***  trks_weights1 = Get3DArray <float> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  float***  trks_weights2 = Get3DArray <float> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  int**     trks_counts_inPhi   = Get2DArray <int> (maxNPtchBins, 40);
  float**   trks_weights1_inPhi = Get2DArray <float> (maxNPtchBins, 40);
  float**   trks_weights2_inPhi = Get2DArray <float> (maxNPtchBins, 40);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("event_weight", &event_weight);
    PbPbTree->SetBranchAddress ("isEE",         &isEE);
    PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
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

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      if (fabs (vz) > 150) continue;

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

      // do a reweighting procedure
      {
        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2)  dphi = pi - dphi;
        if (useCentWgts)  fcal_weight = h_PbPbFCal_weights[iSpc][iPtZ]->GetBinContent (h_PbPbFCal_weights[iSpc][iPtZ]->FindBin (fcal_et));
        if (useQ2Wgts)    q2_weight   = h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->FindBin (q2));
        if (usePsi2Wgts)  psi2_weight = h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->FindBin (dphi));

        event_weight *= fcal_weight * q2_weight * psi2_weight;
      }

      if (event_weight == 0) continue;

      TLorentzVector zvec;
      zvec.SetPxPyPzE (z_pt*cos(z_phi), z_pt*sin(z_phi), sqrt(z_pt*z_pt+z_m*z_m)*sinh(z_y), sqrt(z_pt*z_pt+z_m*z_m)*cosh(z_y));
      z_eta = zvec.Eta ();

      h_z_pt[iCent][iSpc]->Fill (z_pt, event_weight);
      h_z_y_phi[iCent][iSpc][iPtZ]->Fill (z_y, InTwoPi (z_phi), event_weight);
      h_z_eta[iCent][iSpc][iPtZ]->Fill (z_eta, event_weight);
      h_z_y[iCent][iSpc][iPtZ]->Fill (z_y, event_weight);
      int iReg = (fabs (z_y) > 1.00 ? 1 : 0); // barrel vs. endcaps
      h_z_m[iCent][iSpc][iReg]->Fill (z_m, event_weight);

      h_lepton_pt[iCent][iSpc]->Fill (l1_pt, event_weight);
      h_lepton_pt[iCent][iSpc]->Fill (l2_pt, event_weight);
      h_lepton_eta[iCent][iSpc]->Fill (l1_eta, event_weight);
      h_lepton_eta[iCent][iSpc]->Fill (l2_eta, event_weight);
      h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l1_phi), event_weight);
      h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l2_phi), event_weight);

      if (z_pt > 5) {
        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2) dphi = pi - dphi;
        h_z_phi[iCent][iSpc]->Fill (2*dphi, event_weight);
      }

      h_lepton_trk_pt[iCent][iSpc]->Fill (l1_trk_pt, event_weight);
      h_lepton_trk_pt[iCent][iSpc]->Fill (l2_trk_pt, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (2.5, pow (event_weight, 2));

      if (iPtZ < 2) continue;

      h_fcal_et->Fill (fcal_et);
      h_fcal_et_reweighted->Fill (fcal_et, event_weight);

      h_q2[iFineCent]->Fill (q2);
      h_q2_reweighted[iFineCent]->Fill (q2, event_weight);
      h_psi2[iFineCent]->Fill (psi2);
      h_psi2_reweighted[iFineCent]->Fill (psi2, event_weight);
      h_PbPb_vz->Fill (vz);
      h_PbPb_vz_reweighted->Fill (vz, event_weight);

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];
        const float xhz = trkpt / z_pt;

        if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
          continue;
        //if (trk_eta[iTrk] < 0)
        //  continue;

        if (doLeptonRejVar && (DeltaR (l1_trk_eta, trk_eta[iTrk], l1_trk_phi, trk_phi[iTrk]) < 0.02 || DeltaR (l2_trk_eta, trk_eta[iTrk], l2_trk_phi, trk_phi[iTrk]) < 0.02))
          continue;

        {
          float mindr = pi;
          float ptdiff = 0;
          float dr = DeltaR (trk_eta[iTrk], l1_trk_eta, trk_phi[iTrk], l1_trk_phi);
          if (dr < mindr) {
            mindr = dr;
            ptdiff = 2. * fabs (trkpt - l1_trk_pt) / (trkpt + l1_trk_pt);
          }
          dr = DeltaR (trk_eta[iTrk], l2_trk_eta, trk_phi[iTrk], l2_trk_phi);
          if (dr < mindr) {
            mindr = dr;
            ptdiff = 2. * fabs (trkpt - l2_trk_pt) / (trkpt + l2_trk_pt);
          }
          h_lepton_trk_dr[iCent][iSpc]->Fill (mindr, ptdiff);
        }

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], true);
        const double trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], true);
        if (trkEff == 0 || trkPur == 0)
          continue;
        const float trkWeight = trkPur / trkEff;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight * trkWeight);

        // Study correlations (requires dPhi in -pi/2 to 3pi/2)
        float dPhi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dPhi < -pi/2) dPhi = dPhi + 2*pi;

        short iPtch = -1;
        if (allPtchBins[0] <= trkpt) {
          iPtch = 0;
          while (iPtch < maxNPtchBins && allPtchBins[iPtch+1] < trkpt) iPtch++;
        }

        if (iPtch != -1 && iPtch < maxNPtchBins) {
          short idPhi = 0;
          while (idPhi < GetNdPhiBins (0, 0) && (-pi/2.)+(2.*pi/GetNdPhiBins (0, 0))*(idPhi+1) < dPhi) idPhi++;

          trks_counts_inPhi[iPtch][idPhi]   += 1;
          trks_weights1_inPhi[iPtch][idPhi] += trkWeight;
          trks_weights2_inPhi[iPtch][idPhi] += pow (trkWeight, 2);
        }

        short iXhZ = -1;
        if (allXhZBins[0] <= xhz) {
          iXhZ = 0;
          while (iXhZ < maxNXhZBins && allXhZBins[iXhZ+1] < xhz) iXhZ++;
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dPhi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dPhi && dPhi <= phiHighBins[idPhi]) {
            if (iPtch != -1 && iPtch < maxNPtchBins) {
              trks_counts[0][iPtch][idPhi]    += 1;
              trks_weights1[0][iPtch][idPhi]  += trkWeight;
              trks_weights2[0][iPtch][idPhi]  += pow (trkWeight, 2);
            }
            if (iXhZ != -1 && iXhZ < maxNXhZBins) {
              trks_counts[1][iXhZ][idPhi]   += 1;
              trks_weights1[1][iXhZ][idPhi] += trkWeight;
              trks_weights2[1][iXhZ][idPhi] += pow (trkWeight, 2);
            }
          }
        } // end loop over idPhi
        if (3*pi/4 <= dPhi) {
          if (iPtch != -1 && iPtch < maxNPtchBins) {
            trks_counts[0][iPtch][numPhiBins]   += 1;
            trks_weights1[0][iPtch][numPhiBins] += trkWeight;
            trks_weights2[0][iPtch][numPhiBins] += pow (trkWeight, 2);
          }
          if (iXhZ != -1 && iXhZ < maxNXhZBins) {
            trks_counts[1][iXhZ][numPhiBins]    += 1;
            trks_weights1[1][iXhZ][numPhiBins]  += trkWeight;
            trks_weights2[1][iXhZ][numPhiBins]  += pow (trkWeight, 2);
          }
        }
      } // end loop over tracks

      // fill phi correlation histograms and covariance matrices
      for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
        for (int idPhi1 = 0; idPhi1 < GetNdPhiBins (0, 0); idPhi1++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]));
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinError (idPhi1+1, sqrt (pow (h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinError (idPhi1+1), 2) + pow (event_weight, 2) * (trks_weights2_inPhi[iPtch][idPhi1])));
          for (int idPhi2 = 0; idPhi2 < GetNdPhiBins (0, 0); idPhi2++)
            h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, idPhi2+1, h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1, idPhi2+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]) * (trks_weights1_inPhi[iPtch][idPhi2]));
        } // end loop over iPtch
      } // end loop over idPhi1

      // fill yield histograms binned in dPhi and covariance matrices
      for (int idPhi = 0; idPhi < numPhiBins; idPhi++) {
        for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + trks_counts[0][iPtch1][idPhi]);
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]));
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][idPhi])));
          for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
            h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]) * (trks_weights1[0][iPtch2][idPhi]));
        } // end loop over iPtch1
        for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]));
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][idPhi])));
          for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
            h2_trk_xhz_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]) * (trks_weights1[1][iXhZ2][idPhi]));
        } // end loop over iXhZ1
      } // end loop over idPhi

      // fill yield histograms and covariance matrices (for dPhi integrated yield)
      for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]));
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][numPhiBins])));
        for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
          h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]) * (trks_weights1[0][iPtch2][numPhiBins]));
      } // end loop over iPtch1
      for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]));
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][numPhiBins])));
        for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
          h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]) * (trks_weights1[1][iXhZ2][numPhiBins]));
      } // end loop over iXhZ1

      // reset trk count measurements for next event
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < max (maxNPtchBins, maxNXhZBins); j++) {
          for (int k = 0; k <= numPhiBins; k++) {
            trks_counts[i][j][k] = 0;
            trks_weights1[i][j][k] = 0;
            trks_weights2[i][j][k] = 0;
          } // end loop over k
        } // end loop over j
      } // end loop over i
      for (int i = 0; i < maxNPtchBins; i++) {
        for (int j = 0; j < 40; j++) {
          trks_counts_inPhi[i][j] = 0;
          trks_weights1_inPhi[i][j] = 0;
          trks_weights2_inPhi[i][j] = 0;
        } // end loop over j
      } // end loop over i

    } // end loop over Pb+Pb tree
    cout << "Done primary Pb+Pb loop." << endl;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over pp tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (ppTree) {
    ppTree->SetBranchAddress ("event_weight", &event_weight);
    ppTree->SetBranchAddress ("isEE",         &isEE);
    ppTree->SetBranchAddress ("vz",           &vz);
    ppTree->SetBranchAddress ("z_pt",         &z_pt);
    ppTree->SetBranchAddress ("z_y",          &z_y);
    ppTree->SetBranchAddress ("z_phi",        &z_phi);
    ppTree->SetBranchAddress ("z_m",          &z_m);
    ppTree->SetBranchAddress ("l1_pt",        &l1_pt);
    ppTree->SetBranchAddress ("l1_eta",       &l1_eta);
    ppTree->SetBranchAddress ("l1_phi",       &l1_phi);
    ppTree->SetBranchAddress ("l1_charge",    &l1_charge);
    ppTree->SetBranchAddress ("l1_trk_pt",    &l1_trk_pt);
    ppTree->SetBranchAddress ("l1_trk_eta",   &l1_trk_eta);
    ppTree->SetBranchAddress ("l1_trk_phi",   &l1_trk_phi);
    ppTree->SetBranchAddress ("l2_pt",        &l2_pt);
    ppTree->SetBranchAddress ("l2_eta",       &l2_eta);
    ppTree->SetBranchAddress ("l2_phi",       &l2_phi);
    ppTree->SetBranchAddress ("l2_charge",    &l2_charge);
    ppTree->SetBranchAddress ("l2_trk_pt",    &l2_trk_pt);
    ppTree->SetBranchAddress ("l2_trk_eta",   &l2_trk_eta);
    ppTree->SetBranchAddress ("l2_trk_phi",   &l2_trk_phi);
    ppTree->SetBranchAddress ("ntrk",         &ntrk);
    ppTree->SetBranchAddress ("trk_pt",       &trk_pt);
    ppTree->SetBranchAddress ("trk_eta",      &trk_eta);
    ppTree->SetBranchAddress ("trk_phi",      &trk_phi);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      if (fabs (vz) > 150) continue;

      if (event_weight == 0) continue;

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp
      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin
      if (iPtZ < 0 || iPtZ > nPtZBins-1) continue;

      h_pp_nch->Fill (ntrk);
      h_pp_nch_reweighted->Fill (ntrk, event_weight);

      h_pp_vz->Fill (vz);
      h_pp_vz_reweighted->Fill (vz, event_weight);

      TLorentzVector zvec;
      zvec.SetPxPyPzE (z_pt*cos(z_phi), z_pt*sin(z_phi), sqrt(z_pt*z_pt+z_m*z_m)*sinh(z_y), sqrt(z_pt*z_pt+z_m*z_m)*cosh(z_y));
      z_eta = zvec.Eta ();

      h_z_pt[iCent][iSpc]->Fill (z_pt, event_weight);
      h_z_y_phi[iCent][iSpc][iPtZ]->Fill (z_y, InTwoPi (z_phi), event_weight);
      h_z_eta[iCent][iSpc][iPtZ]->Fill (z_eta, event_weight);
      h_z_y[iCent][iSpc][iPtZ]->Fill (z_y, event_weight);
      int iReg = (fabs (z_y) > 1.00 ? 1 : 0); // barrel vs. endcaps
      h_z_m[iCent][iSpc][iReg]->Fill (z_m, event_weight);

      h_lepton_pt[iCent][iSpc]->Fill (l1_pt, event_weight);
      h_lepton_pt[iCent][iSpc]->Fill (l2_pt, event_weight);
      h_lepton_eta[iCent][iSpc]->Fill (l1_eta, event_weight);
      h_lepton_eta[iCent][iSpc]->Fill (l2_eta, event_weight);
      h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l1_phi), event_weight);
      h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l2_phi), event_weight);

      if (z_pt > 5) {
        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2) dphi = pi - dphi;
        h_z_phi[iCent][iSpc]->Fill (2*dphi, event_weight);
      }

      h_lepton_trk_pt[iCent][iSpc]->Fill (l1_trk_pt, event_weight);
      h_lepton_trk_pt[iCent][iSpc]->Fill (l2_trk_pt, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (2.5, pow (event_weight, 2));

      if (iPtZ < 2) continue;

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];
        const float xhz = trkpt / z_pt;

        if (trkpt < trk_min_pt) continue;

        if (doLeptonRejVar && (DeltaR (l1_trk_eta, trk_eta[iTrk], l1_trk_phi, trk_phi[iTrk]) < 0.02 || DeltaR (l2_trk_eta, trk_eta[iTrk], l2_trk_phi, trk_phi[iTrk]) < 0.02)) continue;

        {
          float mindr = pi;
          float ptdiff = 0;
          float dr = DeltaR (trk_eta[iTrk], l1_trk_eta, trk_phi[iTrk], l1_trk_phi);
          if (dr < mindr) {
            mindr = dr;
            ptdiff = 2. * fabs (trkpt - l1_trk_pt) / (trkpt + l1_trk_pt);
          }
          dr = DeltaR (trk_eta[iTrk], l2_trk_eta, trk_phi[iTrk], l2_trk_phi);
          if (dr < mindr) {
            mindr = dr;
            ptdiff = 2. * fabs (trkpt - l2_trk_pt) / (trkpt + l2_trk_pt);
          }
          h_lepton_trk_dr[iCent][iSpc]->Fill (mindr, ptdiff);
        }

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], false);
        const double trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], false);
        if (trkEff == 0 || trkPur == 0) continue;
        const float trkWeight = trkPur / trkEff;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight * trkWeight);

        // Study correlations (requires dPhi in -pi/2 to 3pi/2)
        float dPhi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dPhi < -pi/2) dPhi = dPhi + 2*pi;

        short iPtch = -1;
        if (allPtchBins[0] <= trkpt) {
          iPtch = 0;
          while (iPtch < maxNPtchBins && allPtchBins[iPtch+1] < trkpt) iPtch++;
        }

        if (iPtch != -1 && iPtch < maxNPtchBins) {
          short idPhi = 0;
          while (idPhi < GetNdPhiBins (0, 0) && (-pi/2.)+(2.*pi/GetNdPhiBins (0, 0))*(idPhi+1) < dPhi) idPhi++;

          trks_counts_inPhi[iPtch][idPhi]   += 1;
          trks_weights1_inPhi[iPtch][idPhi] += trkWeight;
          trks_weights2_inPhi[iPtch][idPhi] += pow (trkWeight, 2);
        }

        short iXhZ = -1;
        if (allXhZBins[0] <= xhz) {
          iXhZ = 0;
          while (iXhZ < maxNXhZBins && allXhZBins[iXhZ+1] < xhz) iXhZ++;
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dPhi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dPhi && dPhi <= phiHighBins[idPhi]) {
            if (iPtch != -1 && iPtch < maxNPtchBins) {
              trks_counts[0][iPtch][idPhi]    += 1;
              trks_weights1[0][iPtch][idPhi]  += trkWeight;
              trks_weights2[0][iPtch][idPhi]  += pow (trkWeight, 2);
            }
            if (iXhZ != -1 && iXhZ < maxNXhZBins) {
              trks_counts[1][iXhZ][idPhi]   += 1;
              trks_weights1[1][iXhZ][idPhi] += trkWeight;
              trks_weights2[1][iXhZ][idPhi] += pow (trkWeight, 2);
            }
          }
        } // end loop over idPhi
        if (3*pi/4 <= dPhi) {
          if (iPtch != -1 && iPtch < maxNPtchBins) {
            trks_counts[0][iPtch][numPhiBins]   += 1;
            trks_weights1[0][iPtch][numPhiBins] += trkWeight;
            trks_weights2[0][iPtch][numPhiBins] += pow (trkWeight, 2);
          }
          if (iXhZ != -1 && iXhZ < maxNXhZBins) {
            trks_counts[1][iXhZ][numPhiBins]    += 1;
            trks_weights1[1][iXhZ][numPhiBins]  += trkWeight;
            trks_weights2[1][iXhZ][numPhiBins]  += pow (trkWeight, 2);
          }
        }
      } // end loop over tracks

      // fill phi correlation histograms and covariance matrices
      for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
        for (int idPhi1 = 0; idPhi1 < GetNdPhiBins (0, 0); idPhi1++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]));
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinError (idPhi1+1, sqrt (pow (h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinError (idPhi1+1), 2) + pow (event_weight, 2) * (trks_weights2_inPhi[iPtch][idPhi1])));
          for (int idPhi2 = 0; idPhi2 < GetNdPhiBins (0, 0); idPhi2++)
            h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, idPhi2+1, h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1, idPhi2+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]) * (trks_weights1_inPhi[iPtch][idPhi2]));
        } // end loop over iPtch
      } // end loop over idPhi1

      // fill yield histograms binned in dPhi and covariance matrices
      for (int idPhi = 0; idPhi < numPhiBins; idPhi++) {
        for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + trks_counts[0][iPtch1][idPhi]);
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]));
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][idPhi])));
          for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
            h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]) * (trks_weights1[0][iPtch2][idPhi]));
        } // end loop over iPtch1
        for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]));
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][idPhi])));
          for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
            h2_trk_xhz_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]) * (trks_weights1[1][iXhZ2][idPhi]));
        } // end loop over iXhZ1
      } // end loop over idPhi

      // fill yield histograms and covariance matrices (for dPhi integrated yield)
      for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]));
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][numPhiBins])));
        for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
          h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]) * (trks_weights1[0][iPtch2][numPhiBins]));
      } // end loop over iPtch1
      for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]));
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][numPhiBins])));
        for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
          h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]) * (trks_weights1[1][iXhZ2][numPhiBins]));
      } // end loop over iXhZ1

      // reset trk count measurements for next event
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < max (maxNPtchBins, maxNXhZBins); j++) {
          for (int k = 0; k <= numPhiBins; k++) {
            trks_counts[i][j][k] = 0;
            trks_weights1[i][j][k] = 0;
            trks_weights2[i][j][k] = 0;
          } // end loop over k
        } // end loop over j
      } // end loop over i
      for (int i = 0; i < maxNPtchBins; i++) {
        for (int j = 0; j < 40; j++) {
          trks_counts_inPhi[i][j] = 0;
          trks_weights1_inPhi[i][j] = 0;
          trks_weights2_inPhi[i][j] = 0;
        } // end loop over j
      } // end loop over i

    } // end loop over pp tree
    cout << "Done primary pp loop." << endl;
  }

  Delete3DArray (&trks_counts, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete3DArray (&trks_weights1, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete3DArray (&trks_weights2, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete2DArray (&trks_counts_inPhi, maxNPtchBins, 40);
  Delete2DArray (&trks_weights1_inPhi, maxNPtchBins, 40);
  Delete2DArray (&trks_weights2_inPhi, maxNPtchBins, 40);

  SaveHists (outFileName);

  if (inFile) inFile->Close ();
  SaferDelete (&inFile);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot lepton Pt spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotLeptonPtSpectra (FullAnalysis* a) {

  if (a)
    a->PlotLeptonPtSpectra ();

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    const char* spc = iSpc == 0 ? "e" : "#mu";

    //c->cd (iSpc+1);
    //gPad->SetLogy ();

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      const char* canvasName = Form ("c_%s_pt_iCent%i", iSpc == 0 ? "electron" : (iSpc == 1 ? "muon" : "lepton"), iCent);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      TPad* uPad = nullptr, *dPad = nullptr;
      if (canvasExists) {
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
        uPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_uPad", canvasName)));
        dPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_dPad", canvasName)));
      }
      else {
        c = new TCanvas (canvasName, "", 600, 600);
        c->cd ();
        uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 0.4, 1.0, 1.0);
        dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, 0.4);
        uPad->SetBottomMargin (0);
        dPad->SetTopMargin (0);
        dPad->SetBottomMargin (0.25);
        uPad->Draw ();
        dPad->Draw ();
        gDirectory->Add (c);
        gDirectory->Add (uPad);
        gDirectory->Add (dPad);
        c->cd ();
      }

      uPad->cd ();
      uPad->SetLogx ();
      uPad->SetLogy ();

      TH1D* h = (TH1D*) h_lepton_pt[iCent][iSpc];

      if (!canvasExists) {
        TH1D* htemp = (TH1D*) h->Clone ("htemp");
        htemp->Reset ();
        htemp->GetYaxis ()->SetRangeUser (2e-6, 0.50);
        htemp->GetXaxis ()->SetMoreLogLabels ();
        htemp->GetXaxis ()->SetTitle (Form ("#it{p}_{T}^{ %s} [GeV]", spc));
        htemp->GetYaxis ()->SetTitle (Form ("(1/N_{Z}) (dN_{%s}/d#it{p}_{T}) [GeV^{-1}]", spc));
        htemp->GetXaxis ()->SetTitleSize (0);
        htemp->GetYaxis ()->SetTitleSize (0.04/0.6);
        htemp->GetXaxis ()->SetLabelSize (0);
        htemp->GetYaxis ()->SetLabelSize (0.04/0.6);
        htemp->GetYaxis ()->SetTitleOffset (1.1);
        htemp->DrawCopy ("hist");
        SaferDelete (&htemp);
      }

      if (plotFill) {
        //h->SetFillColorAlpha (fillColors[iCent], fillAlpha);
        h->SetFillColorAlpha (kAzure+10, fillAlpha);
        h->SetLineColor (kBlack);
        h->SetMarkerSize (0);
        h->SetLineWidth (0);

        h->DrawCopy ("bar same");
        h->SetLineWidth (1);
        h->Draw ("hist same");

        gPad->RedrawAxis ();
      } else {
        TGraphAsymmErrors* g = GetTGAE (h);
        ResetXErrors (g);

        const int markerStyle = kFullCircle;
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (1.25);
        g->SetLineWidth (2);
        g->SetLineColor (kBlack);
        g->SetMarkerColor (kBlack);
        //g->SetLineColor (colors[iCent]);
        //g->SetMarkerColor (colors[iCent]);

        g->Draw ("P");
      }

      if (!a) continue;

      h = (TH1D*) h_lepton_pt[iCent][iSpc]->Clone (Form ("h_lepton_pt_ratio_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h->Divide (a->h_lepton_pt[iCent][iSpc]);
        
      dPad->cd ();
      dPad->SetLogx ();
      //dPad->SetLogy ();

      TH1D* htemp = (TH1D*) h->Clone ("htemp");
      htemp->Reset ();
      htemp->GetXaxis ()->SetMoreLogLabels ();
      htemp->GetYaxis ()->SetRangeUser (0.1, 1.9);
      htemp->GetXaxis ()->SetTitle (Form ("#it{p}_{T}^{ %s} [GeV]", spc));
      htemp->GetYaxis ()->SetTitle ("Data / MC");
      htemp->GetXaxis ()->SetTitleSize (0.04/0.4);
      htemp->GetYaxis ()->SetTitleSize (0.04/0.4);
      htemp->GetXaxis ()->SetLabelSize (0.04/0.4);
      htemp->GetYaxis ()->SetLabelSize (0.04/0.4);
      htemp->GetXaxis ()->SetTitleOffset (1.2);
      htemp->GetYaxis ()->SetTitleOffset (1.1);
      htemp->GetYaxis ()->CenterTitle ();
      htemp->DrawCopy ("hist");
      SaferDelete (&htemp);

      TLine* l = new TLine (20, 1, 250, 1);
      l->SetLineColor (46);
      l->SetLineWidth (2);
      l->SetLineStyle (5);
      l->Draw ("same");

      if (h) {
        TGraphAsymmErrors* g = GetTGAE (h);
        ResetXErrors (g);

        const int markerStyle = kFullCircle;
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (1.25);
        g->SetLineWidth (2);
        g->SetLineColor (kBlack);
        g->SetMarkerColor (kBlack);

        g->Draw ("P");
      }
      else {
        cout << "Warning in FullAnalysis :: PlotLeptonPtSpectra: Lepton pT spectra ratio needs to be calculated!" << endl;
      }

      uPad->cd ();
      myText (0.22, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.045/0.6);
      if (iCent == 0) myText (0.68, 0.85, kBlack, "#it{pp}, 5.02 TeV", 0.04/0.6);
      else            myText (0.68, 0.85, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04/0.6);
      const char* zstr = iSpc == 0 ? "#it{Z} #rightarrow #it{e}#it{e}" : (iSpc == 1 ? "#it{Z} #rightarrow #it{#mu}#it{#mu}" : "#it{Z} #rightarrow #it{l}#it{l}");
      myText (0.71, 0.76, kBlack, zstr, 0.04/0.6);

      //myMarkerText (0.753, 0.67, kBlack, kFullCircle, "Data", 1.25, 0.04/0.6);
      myMarkerText (0.78, 0.67, kBlack, kFullCircle, "Data", 1.25, 0.04/0.6);
      //myOnlyBoxText (0.76, 0.58, 1.2, kAzure+10, kBlack, 1, "MC", 0.04/0.6, 1001, 1);
      myOnlyBoxText (0.78, 0.58, 1.2, kAzure+10, kBlack, 1, "MC", 0.04/0.6, 1001, 1);

      c->SaveAs (Form ("%s/LeptonPtSpectra/%sPtSpectra_iCent%i.pdf", plotPath.Data (), iSpc == 0 ? "Electron" : "Muon", iCent));
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot lepton eta spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotLeptonEtaSpcComp (FullAnalysis* a) {
  const char* canvasName = "c_lepton_eta_SpcComp";
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  TPad* uPad = nullptr, *dPad = nullptr;
  if (canvasExists) {
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName)); 
    uPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_uPad", canvasName)));
    dPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_dPad", canvasName)));
  }
  else {
    c = new TCanvas (canvasName, "", 800, 800);
    c->cd ();
    uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 0.4, 1.0, 1.0);
    dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, 0.4);
    uPad->SetBottomMargin (0);
    dPad->SetTopMargin (0);
    dPad->SetBottomMargin (0.25);
    uPad->Draw ();
    dPad->Draw ();
    gDirectory->Add (c);
    gDirectory->Add (uPad);
    gDirectory->Add (dPad);
    c->cd ();
  }

  const short iCent = 0;

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    TH1D* h = (TH1D*) h_lepton_eta[iCent][iSpc]->Clone ();

    uPad->cd ();
    if (plotFill) {
      h->SetFillColorAlpha (fillColors[iSpc+1], 0.2);
      h->SetLineColor (kBlack);
      h->SetMarkerSize (0);
      h->SetLineWidth (0);
      //h->GetYaxis ()->SetRangeUser (0, 1.3);
      h->GetYaxis ()->SetRangeUser (0, 1.0);

      h->GetXaxis ()->SetTitle ("Lepton #eta");
      //h->GetYaxis ()->SetTitle ("Arb. Units");
      h->GetYaxis ()->SetTitle ("1/N_{Z} dN_{l}/d#eta");
      h->GetXaxis ()->SetTitleSize (0.04/0.6);
      h->GetYaxis ()->SetTitleSize (0.04/0.6);
      h->GetXaxis ()->SetLabelSize (0.04/0.6);
      h->GetYaxis ()->SetLabelSize (0.04/0.6);
      h->GetXaxis ()->SetTitleOffset (1.5*0.6);
      h->GetYaxis ()->SetTitleOffset (1.5*0.6);

      h->DrawCopy (iSpc == 0 && !canvasExists ? "bar" : "bar same");
      h->SetLineWidth (1);
      h->Draw ("hist same");

      gPad->RedrawAxis ();
    }
    else {
      TGraphAsymmErrors* g = GetTGAE (h);
      ResetXErrors (g);
      //deltaize (g, 0.1*(-1.5+iCent));

      const int markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
      g->SetMarkerStyle (markerStyle);
      g->SetMarkerSize (1);
      g->SetLineWidth (1);
      g->SetLineColor (colors[iSpc+1]);
      g->SetMarkerColor (colors[iSpc+1]);
      g->GetYaxis ()->SetRangeUser (0, 1.0);

      g->GetXaxis ()->SetTitle ("Lepton #eta");
      g->GetYaxis ()->SetTitle ("1/N_{Z} dN_{l}/d#eta");
      g->GetXaxis ()->SetTitleSize (0.04/0.6);
      g->GetYaxis ()->SetTitleSize (0.04/0.6);
      g->GetXaxis ()->SetLabelSize (0.04/0.6);
      g->GetYaxis ()->SetLabelSize (0.04/0.6);
      g->GetXaxis ()->SetTitleOffset (1.5*0.6);
      g->GetYaxis ()->SetTitleOffset (1.5*0.6);
      g->Draw (iSpc == 0 && !canvasExists ? "AP" : "P");
    }

    if (a) {
      dPad->cd ();
      //while (a->h_lepton_eta[iCent][iSpc]->GetNbinsX () >= 2*h->GetNbinsX ()) {
      //  a->h_lepton_eta[iCent][iSpc]->Rebin (2);
      //  a->h_lepton_eta[iCent][iSpc]->Scale (0.5);
      //}
      h->Divide (a->h_lepton_eta[iCent][iSpc]);
      if (h) {
        TGraphAsymmErrors* g = GetTGAE (h);
        ResetXErrors (g);

        const int markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (1);
        g->SetLineWidth (1);
        g->SetLineColor (colors[iSpc+1]);
        g->SetMarkerColor (colors[iSpc+1]);
        //g->GetYaxis ()->SetRangeUser (0.76, 1.24);
        g->GetYaxis ()->SetRangeUser (0.9, 1.1);

        g->GetXaxis ()->SetTitle ("Lepton #eta");
        g->GetYaxis ()->SetTitle ("Data / MC Reco");
        g->GetXaxis ()->SetTitleSize (0.04/0.4);
        g->GetYaxis ()->SetTitleSize (0.04/0.4);
        g->GetXaxis ()->SetLabelSize (0.04/0.4);
        g->GetYaxis ()->SetLabelSize (0.04/0.4);
        g->GetXaxis ()->SetTitleOffset (2.5*0.4);
        g->GetYaxis ()->SetTitleOffset (1.5*0.4);

        g->GetYaxis ()->CenterTitle ();
        g->Draw (iSpc == 0 && !canvasExists ? "AP" : "P");

        if (iSpc == 0) {
          TLine* l = new TLine (-2.5, 1, 2.5, 1);
          l->SetLineColor (46);
          l->SetLineWidth (2);
          l->SetLineStyle (5);
          l->Draw ("same");
        }
      }
    }
  }
    
  uPad->cd ();
  myText (0.22, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.045/0.6);
  if (iCent == 0)
    myText (0.22, 0.80, kBlack, Form ("#it{pp}, 5.02 TeV"), 0.04/0.6);
  else
    myText (0.22, 0.83-iCent*0.06, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04/0.6);

  myText (0.58, 0.90, kBlack, "Data", 0.032/0.6);
  myText (0.66, 0.90, kBlack, "MC Reco.", 0.032/0.6);
  for (int iSpc = 0; iSpc < 2; iSpc++) {
    const char* spcLabel = iSpc == 0 ? "#it{Z} #rightarrow #it{e}#it{e}" : (iSpc == 1 ? "#it{Z} #rightarrow #it{#mu}#it{#mu}" : "#it{Z} #rightarrow #it{l}#it{l}");
    myMarkerTextNoLine (0.64, 0.823-iSpc*0.06, colors[iSpc+1], kFullCircle, "", 1.8, 0.04/0.6);
    //myBoxText (0.72, 0.823-iSpc*0.06, colors[iSpc+1], kOpenCircle, "", 1.5, 0.04/0.6);
    myOnlyBoxText (0.75, 0.823-iSpc*0.06, 1.2, fillColors[iSpc+1], kBlack, 1, "", 0.06, 1001, 1.);
    myText (0.79, 0.82-iSpc*0.06, colors[iSpc+1], spcLabel, 0.04/0.6);
  }

  c->SaveAs (Form ("%s/LeptonEtaSpcComp.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot track Pt spectra for each lepton species
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotLeptonTrackPtSpectra () {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* canvasName = Form ("c_%s_trk_pt", iSpc == 0 ? "electron" : (iSpc == 1 ? "muon" : "lepton"));
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 800, 600);
      gDirectory->Add (c);
    }
    c->cd ();
    c->SetLogy ();

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      //TH1D* h = (TH1D*)h_lepton_trk_pt[iCent][iSpc]->Clone (Form ("h_lepton_trk_pt_iSpc%i_iCent%i", iSpc, iCent));
      TH1D* h = h_lepton_trk_pt[iCent][iSpc];

      //h->Rebin (5);
      //h->Scale (1./h_z_pt[iCent][iSpc]->Integral (h_z_pt[iCent][iSpc]->GetXaxis ()->FindBin (5), h_z_pt[iCent][iSpc]->GetNbinsX ()), "width");
      //cout << iSpc << ", " << iCent << ", " << h_z_pt[iCent][iSpc]->Integral (h_z_pt[iCent][iSpc]->GetXaxis ()->FindBin (5), h_z_pt[iCent][iSpc]->GetNbinsX ()) << endl;

      h->GetYaxis ()->SetRangeUser (6e-6, 450);

      h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
      if (iSpc == 0) {
        h->GetYaxis ()->SetTitle ("(1/N_{Z#rightarrow#it{ee}}) (dN_{ch}/d#it{p}_{T}) [GeV^{-1}]");
      }
      else if (iSpc == 1) {
        h->GetYaxis ()->SetTitle ("(1/N_{Z#rightarrow#it{#mu#mu}}) (dN_{ch}/d#it{p}_{T}) [GeV^{-1}]");
      }
      else {
        h->GetYaxis ()->SetTitle ("(1/N_{Z#rightarrow#it{ll}}) (dN_{ch}/d#it{p}_{T}) [GeV^{-1}]");
      }

      h->SetLineColor (colors[iCent]);
      h->SetMarkerColor (colors[iCent]);
      h->SetMarkerStyle (kFullCircle);
      h->SetMarkerSize (0.75);

      h->Draw (iCent == 0 ? "e1" : "e1 same");
    }
    gPad->RedrawAxis ();

    myText (0.66, 0.82, colors[0], "#it{pp}", 0.04);
    for (short iCent = 1; iCent < numCentBins; iCent++) {
      myText (0.66, 0.82-0.06*iCent, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);
    }

    myText (0.66, 0.88, kBlack, iSpc == 0 ? "#it{Z} #rightarrow #it{e}#it{e}" : (iSpc == 1 ?"#it{Z} #rightarrow #it{#mu}#it{#mu}" : "#it{Z} #rightarrow #it{l}#it{l}"), 0.04);
    c->SaveAs (Form ("%s/LeptonTrackPtSpectra/%sTrackPtSpectra.pdf", plotPath.Data (), iSpc == 0 ? "Electron" : (iSpc == 1 ? "Muon" : "Comb")));
  }

}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot dR between leptons and tracks
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotLeptonTrackDR () {
  
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      const char* canvasName = Form ("c_%s_trk_dr", iSpc == 0 ? "electron" : (iSpc == 1 ? "muon" : "lepton"));
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 800, 600);
        gDirectory->Add (c);
        FormatTH2Canvas (c, true);
        c->SetLogz ();
      }
      c->cd ();

      TH2D* h2 = h_lepton_trk_dr[iCent][iSpc];
      //h2->RebinX (2);
      //h2->RebinY (2);
      h2->GetXaxis ()->SetTitle (Form ("min (#DeltaR (track, %s))", iSpc == 0 ? "electrons" : (iSpc == 1 ? "muons" : "leptons")));
      //h2->GetYaxis ()->SetTitle ("#it{p}_{T}^{h} [GeV]");
      h2->GetYaxis ()->SetTitle ("|#Delta#it{p}_{T}| / <#it{p}_{T}>");
      //h2->GetYaxis ()->SetTitle ("#Delta#phi");
      h2->GetZaxis ()->SetTitle ("Counts");

      h2->GetYaxis ()->SetTitleOffset (1.1);

      h2->Draw ("colz");

      myText (0.56, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);

      if (iCent != 0) {
        myText (0.56, 0.82, kBlack, "Pb+Pb, 5.02 TeV", 0.04);
        myText (0.56, 0.76, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);
      }
      else
        myText (0.56, 0.82, kBlack, "#it{pp}, 5.02 TeV", 0.04);

      c->SaveAs (Form ("%s/LeptonTrackDists/%sTrackDist_iCent%i.pdf", plotPath.Data (), iSpc == 0 ? "Electron" : (iSpc == 1 ? "Muon" : "Lepton"), iCent));
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot dR between leptons and tracks
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotLeptonTrackDRProjX () {
  
  const char* canvasName = Form ("c_trk_dr");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1400, 600);
    gDirectory->Add (c);
    c->Divide (2, 1);
  }

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    c->cd (iSpc+1);
    //gPad->SetLogy ();

    TH1D** harr = new TH1D*[numCentBins];
    for (short iCent = 0; iCent < numCentBins; iCent++) {

      TH2D* h2 = h_lepton_trk_dr[iCent][iSpc];

      TH1D* h = h2->ProjectionX ();
      harr[iCent] = h;

      int nEvt = 0;
      for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        nEvt += h_z_counts[iSpc][iPtZ][iCent]->GetBinContent (1);
      }
      h->Scale (1./nEvt);
      

      //TH1D* h = h2->ProjectionX ("temp0", h2->GetYaxis ()->FindBin (2), h2->GetYaxis ()->FindBin (3)-1);
      //h->Rebin (2);
      h->GetXaxis ()->SetTitle (Form ("Minimum #DeltaR (track, %s)", iSpc == 0 ? "electrons" : (iSpc == 1 ? "muons" : "leptons")));
      h->GetYaxis ()->SetTitle ("Tracks / Event");
      h->SetLineColor (finalColors[iCent]);
      h->SetMarkerColor (finalColors[iCent]);
      h->GetXaxis ()->SetRangeUser (0, 0.16);
      //h->GetYaxis ()->SetRangeUser (0.5, 0.5e4);
      h->GetYaxis ()->SetRangeUser (0, 1);
      h->GetYaxis ()->SetTitleOffset (1.1);

      h->DrawCopy (iCent == 0 ? "e1" : "e1 same");
      //delete h;

      if (iCent != 0) {
        myText (0.22, 0.80-iCent*0.06, finalColors[iCent], Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);
      }
      else {
        myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
        myText (0.22, 0.80, kBlack, "#it{pp}, 5.02 TeV", 0.04);
      }
    }

    TLine* l1, *l2;
    l1 = new TLine (0.05, 0, 0.05, 0.15);
    l2 = new TLine (0, 0.15, 0.05, 0.15);
    l1->SetLineStyle (2);
    l2->SetLineStyle (2);
    l1->Draw ("same");
    l2->Draw ("same");

    TPad* subpad = new TPad (Form ("p_%s", iSpc == 0 ? "ee":"mumu"), "", 0.5, 0.5, 1.00, 1.00);
    subpad->SetBorderSize (1);
    subpad->Draw ();
    subpad->cd ();

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      TH1D* h = harr[iCent];
      h->GetXaxis ()->SetTitle ("");
      h->GetYaxis ()->SetTitle ("");
      h->GetXaxis ()->SetRangeUser (0, 0.05);
      h->GetYaxis ()->SetRangeUser (0, 0.15);
      h->GetXaxis ()->SetNdivisions (5, 8, false);
      h->Draw (iCent == 0 ? "e1" : "e1 same");
    }
  }
  c->SaveAs (Form ("%s/LeptonTrackDists/LeptonTrackDist_ProjX.pdf", plotPath.Data ()));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z Pt spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZPtSpectra (FullAnalysis* a) {

  if (a)
    a->PlotZPtSpectra ();

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      const char* canvasName = Form ("c_z_pt_iCent%i_%s", iCent, spc);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      TPad* uPad = nullptr, *dPad = nullptr;
      if (canvasExists) {
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
        uPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_uPad", canvasName)));
        dPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_dPad", canvasName)));
      }
      else {
        c = new TCanvas (canvasName, "", 800, 800);
        c->cd ();
        uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 0.4, 1.0, 1.0);
        dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, 0.4);
        uPad->SetBottomMargin (0);
        dPad->SetTopMargin (0);
        dPad->SetBottomMargin (0.25);
        uPad->Draw ();
        dPad->Draw ();
        gDirectory->Add (c);
        gDirectory->Add (uPad);
        gDirectory->Add (dPad);
      }
      c->cd ();

      uPad->cd ();
      gPad->SetLogx ();
      gPad->SetLogy ();

      TH1D* h = h_z_pt[iCent][iSpc];

      if (!canvasExists) {
        TH1D* htemp = (TH1D*) h->Clone ("htemp");
        htemp->Reset ();
        htemp->GetYaxis ()->SetRangeUser (1.2e-6, 0.90);
        htemp->GetXaxis ()->SetTitle ("#it{p}_{T}^{ Z} [GeV]");
        htemp->GetYaxis ()->SetTitle ("(1/N_{Z}) (dN/d#it{p}_{T}) [GeV^{-1}]");
        htemp->GetXaxis ()->SetTitleSize (0);
        htemp->GetYaxis ()->SetTitleSize (0.04/0.6);
        htemp->GetXaxis ()->SetLabelSize (0);
        htemp->GetYaxis ()->SetLabelSize (0.04/0.6);
        htemp->GetYaxis ()->SetTitleOffset (1.1);
        htemp->DrawCopy ("hist");
        SaferDelete (&htemp);
      }

      if (plotFill) {
        //h->SetFillColorAlpha (fillColors[iCent], fillAlpha);
        h->SetFillColorAlpha (kAzure+10, fillAlpha);
        h->SetLineColor (kBlack);
        h->SetMarkerSize (0);
        h->SetLineWidth (0);

        h->DrawCopy ("bar same");
        h->SetLineWidth (1);
        h->Draw ("hist same");

        gPad->RedrawAxis ();
      } else {
        TGraphAsymmErrors* g = GetTGAE (h);
        ResetXErrors (g);

        const int markerStyle = kFullCircle;
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (1.25);
        g->SetLineWidth (2);
        g->SetLineColor (kBlack);
        g->SetMarkerColor (kBlack);
        //g->SetLineColor (colors[iCent]);
        //g->SetMarkerColor (colors[iCent]);
        g->Draw ("P");
      }

      if (!a) continue;

      dPad->cd ();
      gPad->SetLogx ();

      h = (TH1D*) h_z_pt[iCent][iSpc]->Clone (Form ("h_z_pt_ratio_%s_iCent%i_%s", spc, iCent, name.c_str ()));
      h->Divide (a->h_z_pt[iCent][iSpc]);
      h_z_pt_ratio[iCent][iSpc] = h;

      TH1D* htemp = (TH1D*) h->Clone ("htemp");
      htemp->Reset ();
      htemp->GetXaxis ()->SetMoreLogLabels ();
      htemp->GetYaxis ()->SetRangeUser (0.1, 1.9);
      htemp->GetXaxis ()->SetTitle ("#it{p}_{T}^{Z} [GeV]");
      htemp->GetYaxis ()->SetTitle ("Data / MC");
      htemp->GetXaxis ()->SetTitleSize (0.04/0.4);
      htemp->GetYaxis ()->SetTitleSize (0.04/0.4);
      htemp->GetXaxis ()->SetLabelSize (0.04/0.4);
      htemp->GetYaxis ()->SetLabelSize (0.04/0.4);
      htemp->GetXaxis ()->SetTitleOffset (1.2);
      htemp->GetYaxis ()->SetTitleOffset (1.1);
      htemp->GetYaxis ()->CenterTitle ();
      htemp->DrawCopy ("hist");
      SaferDelete (&htemp);

      TLine* l = new TLine (0, 1, 300, 1);
      l->SetLineColor (46);
      l->SetLineWidth (2);
      l->SetLineStyle (5);
      l->Draw ("same");

      if (h) {
        TGraphAsymmErrors* g = GetTGAE (h);
        ResetXErrors (g);

        const int markerStyle = kFullCircle;
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (1);
        g->SetLineWidth (2);
        g->SetLineColor (kBlack);
        g->SetMarkerColor (kBlack);
        g->Draw ("P");
      }
      else {
        cout << "Warning in FullAnalysis :: PlotZPtSpectra: Z pT spectra ratio not stored, needs to be calculated!" << endl;
      }

      uPad->cd ();

      myText (0.62, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.045/0.6);
      myText (0.26, 0.85, kBlack, Form ("#it{Z} #rightarrow %s", iSpc == 0 ? "#it{e}#it{e}" : (iSpc == 1 ? "#it{#mu}#it{#mu}" : "#it{l}#it{l}")), 0.04/0.6);
      //myMarkerText (0.753, 0.65, kBlack, kFullCircle, "Data", 1.25, 0.04/0.6);
      myMarkerText (0.81, 0.65, kBlack, kFullCircle, "Data", 1.25, 0.04/0.6);
      //myOnlyBoxText (0.76, 0.55, 1.2, kAzure+10, kBlack, 1, "MC", 0.04/0.6, 1001, 1);
      myOnlyBoxText (0.81, 0.58, 1.2, kAzure+10, kBlack, 1, "MC", 0.04/0.6, 1001, 1);

      if (iCent == 0) myText (0.62, 0.75, kBlack, Form ("#it{pp}, 5.02 TeV"), 0.04/0.6);
      else            myText (0.62, 0.75, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04/0.6);

      c->SaveAs (Form ("%s/ZPtSpectra/z_pt_spectrum_iCent%i_%s.pdf", plotPath.Data (), iCent, spc));
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z Pt spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZPtSpectraChannelCompare () {
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    const char* canvasName = Form ("c_z_pt_iCent%i_channelCompare", iCent);
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    TPad* uPad = nullptr, *dPad = nullptr;
    if (canvasExists) {
      c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      uPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_uPad", canvasName)));
      dPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_dPad", canvasName)));
    }
    else {
      c = new TCanvas (canvasName, "", 800, 800);
      c->cd ();
      uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 0.4, 1.0, 1.0);
      dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, 0.4);
      uPad->SetBottomMargin (0);
      dPad->SetTopMargin (0);
      dPad->SetBottomMargin (0.25);
      uPad->Draw ();
      dPad->Draw ();
      gDirectory->Add (c);
      gDirectory->Add (uPad);
      gDirectory->Add (dPad);
    }
    c->cd ();

    uPad->cd ();
    gPad->SetLogx ();
    gPad->SetLogy ();

    if (!canvasExists) {
      TH1D* htemp = new TH1D ("htemp", "", 1, 1, 250);
      htemp->Reset ();
      htemp->GetYaxis ()->SetRangeUser (1.2e-6, 0.90);
      htemp->GetXaxis ()->SetTitle ("#it{p}_{T}^{ Z} [GeV]");
      htemp->GetYaxis ()->SetTitle ("(1/N_{Z}) (dN/d#it{p}_{T}) [GeV^{-1}]");
      htemp->GetXaxis ()->SetTitleSize (0);
      htemp->GetYaxis ()->SetTitleSize (0.04/0.6);
      htemp->GetXaxis ()->SetLabelSize (0);
      htemp->GetYaxis ()->SetLabelSize (0.04/0.6);
      htemp->GetYaxis ()->SetTitleOffset (1.1);
      htemp->DrawCopy ("hist");
      SaferDelete (&htemp);
    }

    for (short iSpc : {0, 1}) {
      TH1D* h = h_z_pt[iCent][iSpc];

      TGraphAsymmErrors* g = GetTGAE (h);
      ResetXErrors (g);

      Style_t markerStyle = kFullCircle;
      const float markerSize = 1.25;
      Color_t markerColor = (iSpc == 0 ? finalColors[3] : finalColors[1]);

      g->SetMarkerStyle (markerStyle);
      g->SetMarkerColor (markerColor);
      g->SetMarkerSize (markerSize);
      g->SetLineWidth (2);
      g->SetLineColor (markerColor);
      ((TGAE*) g->Clone ())->Draw ("P");

      markerStyle = FullToOpenMarker (markerStyle);

      g->SetMarkerStyle (markerStyle);
      g->SetMarkerSize (markerSize); 
      if (iCent > 0) g->SetLineWidth (0);
      g->SetMarkerColor (kBlack);
      
      ((TGAE*) g->Clone ())->Draw ("P");

      SaferDelete (&g);
    }

    dPad->cd ();
    gPad->SetLogx ();

    {
      TH1D* htemp = new TH1D ("htemp", "", 1, 1, 250);;
      htemp->GetXaxis ()->SetMoreLogLabels ();
      htemp->GetYaxis ()->SetRangeUser (iCent == 0 ? 0.8: 0.1, iCent == 0 ? 1.2 : 1.9);
      htemp->GetXaxis ()->SetTitle ("#it{p}_{T}^{Z} [GeV]");
      htemp->GetYaxis ()->SetTitle ("#it{ee} / #it{#mu#mu}");
      htemp->GetXaxis ()->SetTitleSize (0.04/0.4);
      htemp->GetYaxis ()->SetTitleSize (0.04/0.4);
      htemp->GetXaxis ()->SetLabelSize (0.04/0.4);
      htemp->GetYaxis ()->SetLabelSize (0.04/0.4);
      htemp->GetXaxis ()->SetTitleOffset (1.2);
      htemp->GetYaxis ()->SetTitleOffset (1.1*0.4/0.6);
      htemp->GetYaxis ()->CenterTitle ();
      htemp->DrawCopy ("hist");
      SaferDelete (&htemp);

      TLine* l = new TLine (0, 1, 300, 1);
      l->SetLineColor (46);
      l->SetLineWidth (2);
      l->SetLineStyle (5);
      l->Draw ("same");
    }

    TH1D* h = (TH1D*) h_z_pt[iCent][0]->Clone (Form ("h_z_pt_ratio_%s_iCent%i_ee_over_mumu", name.c_str (), iCent));
    h->Divide (h_z_pt[iCent][1]);

    TGraphAsymmErrors* g = GetTGAE (h);
    ResetXErrors (g);

    const int markerStyle = kFullCircle;
    g->SetMarkerStyle (markerStyle);
    g->SetMarkerStyle (markerStyle);
    g->SetMarkerSize (1.25);
    g->SetLineWidth (2);
    g->SetLineColor (kBlack);
    g->SetMarkerColor (kBlack);
    g->Draw ("P");

    uPad->cd ();

    myText (0.62, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.045/0.6);
    myMarkerText (0.81, 0.65, finalColors[3], kFullCircle, "#it{Z} #rightarrow #it{ee}", 1.25, 0.04/0.6, true);
    myMarkerText (0.81, 0.58, finalColors[1], kFullCircle, "#it{Z} #rightarrow #it{#mu#mu}", 1.25, 0.04/0.6, true);

    if (iCent == 0) myText (0.62, 0.75, kBlack, Form ("#it{pp}, 5.02 TeV"), 0.04/0.6);
    else            myText (0.62, 0.75, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04/0.6);

    c->SaveAs (Form ("%s/ZPtSpectra/z_pt_spectrum_iCent%i_channelCompare.pdf", plotPath.Data (), iCent));
  }
  
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the per-bin Z pT distribution means.
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: CalculateZPtMeans () {

  double zpt_norm, mean_zpt, mean_zpt_err;
  double*** mean_zpts = Get3DArray <double> (3, numCentBins, nPtZBins);
  double*** mean_zpt_errs = Get3DArray <double> (3, numCentBins, nPtZBins);

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : "mumu");

    for (short iCent = 0; iCent < numCentBins; iCent++) {

      TH1D* h_zpt = h_z_pt[iCent][iSpc];
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        zpt_norm = 0;
        mean_zpt = 0;
        mean_zpt_err = 0;
        const int firstBin = h_zpt->FindBin (zPtBins[iPtZ]);
        const int lastBin = h_zpt->FindBin(zPtBins[iPtZ+1])-1;
        // calculate mean pT^Z
        for (int iX = firstBin; iX <= lastBin; iX++) {
          zpt_norm += h_zpt->GetBinContent (iX) * h_zpt->GetBinWidth (iX);
          mean_zpt += h_zpt->GetBinContent (iX) * h_zpt->GetBinCenter (iX) * h_zpt->GetBinWidth (iX);
        }
        mean_zpt = mean_zpt / zpt_norm;
        // now calculate the derivative vector and assume the covariance matrix is diagonal (events are independent)
        mean_zpt_err = 0;
        for (int iX = firstBin; iX <= lastBin; iX++) {
          mean_zpt_err += pow ((h_zpt->GetBinCenter (iX) - mean_zpt) * h_zpt->GetBinWidth (iX), 2) * pow (h_zpt->GetBinError (iX), 2) / pow (zpt_norm, 2);
        }
        mean_zpt_err = sqrt (mean_zpt_err);

        mean_zpts[iSpc][iCent][iPtZ] = mean_zpt;
        mean_zpt_errs[iSpc][iCent][iPtZ] = mean_zpt_err;

        cout << "<pT^Z> in iCent = " << iCent << ", iPtZ = " << iPtZ << ", iSpc = " << iSpc << " is " << mean_zpts[iSpc][iCent][iPtZ] << " +/- " << mean_zpt_errs[iSpc][iCent][iPtZ] << " GeV" << endl;
      } // end loop over iPtZ
    } // end loop over iCent
  } // end loop over iSpc
}




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Prints yield of Z's that meet the event selection criteria in each centrality
//////////////////////////////////////////////////////////////////////////////////////////////////
//void FullAnalysis :: PrintZYields () {
//  for (short iSpc = 0; iSpc < 3; iSpc++) {
//    const char* spc = (iSpc == 0 ? "Z->ee" : (iSpc == 1 ? "Z->mumu" : "Z->ll"));
//    for (short iCent = 0; iCent < numCentBins; iCent++) {
//      float yield = h_z_counts[iSpc][2][iCent]->GetBinContent (1);
//      if (iCent == 0) 
//        cout << "pp " << spc << " # Z's > 25 GeV  =  " << yield << endl;
//      else
//        cout << Form ("Pb+Pb %i-%i%% ", (int)centCuts[iCent], (int)centCuts[iCent-1]) << spc << " # Z's > 25 GeV  =  " << yield << endl;
//    }
//  }
//}





////////////////////////////////////////////////////////////////////////////////////////////////
// Plots the RAPIDITY-phi distribution of reconstructed Z bosons
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZYPhiMap () {
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      const char* canvasName = Form ("c_z_y_phi_%s_iCent%i", spc, iCent);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 800, 600);
        FormatTH2Canvas (c, true);
        gDirectory->Add (c);
      }
      c->cd ();

      TH2D* h = h_z_y_phi[iCent][iSpc][nPtZBins];

      h->GetXaxis ()->SetTitle ("y_{Z}");
      h->GetYaxis ()->SetTitle ("#phi_{Z}");
      h->GetZaxis ()->SetTitle ("Counts");

      h->GetXaxis ()->SetTitleOffset (1.4);
      h->GetYaxis ()->SetTitleOffset (1.1);

      h->Draw ("colz");

      myText (0.18, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);

      const char* spcLabel = iSpc == 0 ? "#it{Z} #rightarrow #it{e}#it{e}" : (iSpc == 1 ? "#it{Z} #rightarrow #it{#mu}#it{#mu}" : "#it{Z} #rightarrow #it{l}#it{l}");
      myText (0.62, 0.90, kBlack, spcLabel, 0.04);
      myText (0.62, 0.84, kBlack, "#it{p}_{T}^{Z} > 15 GeV", 0.04);

      if (iCent == 0) myText (0.18, 0.84, kBlack, Form ("#it{pp}, 5.02 TeV"), 0.04);
      else            myText (0.18, 0.84, kBlack, Form ("Pb+Pb %i-%i%%, 5.02 TeV", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);

      c->SaveAs (Form ("%s/ZYPhiDists/z%s_y_phi_iCent%i.pdf", plotPath.Data (), spc, iCent));
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots the PSEUDORAPIDITY distribution of reconstructed Z bosons
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZEtaMap (FullAnalysis* a) {

  if (a)
    a->PlotZEtaMap ();

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      const char* canvasName = Form ("c_z_eta_%s_iCent%i", spc, iCent);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      TPad* uPad = nullptr, *dPad = nullptr;
      if (canvasExists) {
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName)); 
        uPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_uPad", canvasName)));
        dPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_dPad", canvasName)));
      }
      else {
        c = new TCanvas (canvasName, "", 800, 800);
        c->cd ();
        uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 0.4, 1.0, 1.0);
        dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, 0.4);
        uPad->SetBottomMargin (0);
        dPad->SetTopMargin (0);
        dPad->SetBottomMargin (0.25);
        uPad->Draw ();
        dPad->Draw ();
        gDirectory->Add (c);
        gDirectory->Add (uPad);
        gDirectory->Add (dPad);
      }

      uPad->cd ();
      TH1D* h = h_z_eta[iCent][iSpc][nPtZBins];

      if (plotFill) {
        //h->SetFillColorAlpha (fillColors[iCent], fillAlpha);
        h->SetFillColorAlpha (kYellow-4, fillAlpha);
        h->SetLineColor (kBlack);
        h->SetMarkerSize (0);
        h->SetLineWidth (0);
        h->GetYaxis ()->SetRangeUser (0, 0.30);

        h->GetXaxis ()->SetTitle ("#eta_{Z}");
        h->GetYaxis ()->SetTitle ("(1/N_{Z}) (dN_{Z}/d#eta)");
        h->GetXaxis ()->SetTitleSize (0.04/0.6);
        h->GetYaxis ()->SetTitleSize (0.04/0.6);
        h->GetXaxis ()->SetLabelSize (0.04/0.6);
        h->GetYaxis ()->SetLabelSize (0.04/0.6);
        h->GetXaxis ()->SetTitleOffset (1.5*0.6);
        h->GetYaxis ()->SetTitleOffset (1.5*0.6);

        h->DrawCopy (!canvasExists ? "bar" : "bar same");
        h->SetLineWidth (1);
        h->Draw ("hist same");

        gPad->RedrawAxis ();
      }
      else {
        TGraphAsymmErrors* g = GetTGAE (h);
        ResetXErrors (g);

        const int markerStyle = kFullCircle;
        g->SetMarkerStyle (markerStyle);
        //g->SetMarkerSize (1);
        g->SetLineWidth (2);
        //g->SetLineColor (colors[iCent]);
        //g->SetMarkerColor (colors[iCent]);
        g->GetYaxis ()->SetRangeUser (0, 0.30);

        g->GetXaxis ()->SetTitle ("#eta_{Z}");
        g->GetYaxis ()->SetTitle ("(1/N_{Z}) (dN_{Z}/d#eta)");
        g->GetXaxis ()->SetTitleSize (0.04/0.6);
        g->GetYaxis ()->SetTitleSize (0.04/0.6);
        g->GetXaxis ()->SetLabelSize (0.04/0.6);
        g->GetYaxis ()->SetLabelSize (0.04/0.6);
        g->GetXaxis ()->SetTitleOffset (1.5*0.6);
        g->GetYaxis ()->SetTitleOffset (1.5*0.6);
        g->Draw (!canvasExists ? "AP" : "P");
      }

      if (!a)
        continue;

      myText (0.22, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.045/0.6);

      const char* spcLabel = iSpc == 0 ? "#it{Z} #rightarrow #it{e}#it{e}" : (iSpc == 1 ? "#it{Z} #rightarrow #it{#mu}#it{#mu}" : "#it{Z} #rightarrow #it{l}#it{l}");
      myText (0.66, 0.86, kBlack, spcLabel, 0.04/0.6);
      myText (0.66, 0.76, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[2]), 0.04/0.6);
      if (iCent == 0) myText (0.22, 0.76, kBlack, Form ("#it{pp}, 5.02 TeV"), 0.04/0.6);
      else            myText (0.22, 0.76, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04/0.6);
      //myMarkerText (0.273, 0.69, kBlack, kFullCircle, "Data", 1.25, 0.04/0.6);
      myMarkerText (0.30, 0.69, kBlack, kFullCircle, "Data", 1.25, 0.04/0.6);
      //myOnlyBoxText (0.28, 0.60, 1.2, kYellow-4, kBlack, 1, "MC", 0.04/0.6, 1001, 1);
      myOnlyBoxText (0.30, 0.60, 1.2, kYellow-4, kBlack, 1, "MC", 0.04/0.6, 1001, 1);


      dPad->cd ();

      h = (TH1D*) h_z_eta[iCent][iSpc][nPtZBins]->Clone (Form ("h_z_eta_ratio_%s_iCent%i_iPtZ%i_%s", spc, iCent, nPtZBins, name.c_str ()));
      h->Divide (a->h_z_eta[iCent][iSpc][nPtZBins]);
      h_z_eta_ratio[iCent][iSpc][nPtZBins] = h;

      if (h) {
        TGraphAsymmErrors* g = GetTGAE (h);
        ResetXErrors (g);

        const int markerStyle = kFullCircle;
        g->SetMarkerStyle (markerStyle);
        //g->SetMarkerSize (1);
        g->SetLineWidth (2);
        //g->SetLineColor (colors[iCent]);
        //g->SetMarkerColor (colors[iCent]);
        g->GetYaxis ()->SetRangeUser (0.66, 1.34);

        g->GetXaxis ()->SetTitle ("#eta_{Z}");
        //g->GetYaxis ()->SetTitle ("Pb+Pb / #it{pp}");
        g->GetYaxis ()->SetTitle ("Data / MC");
        g->GetXaxis ()->SetTitleSize (0.04/0.4);
        g->GetYaxis ()->SetTitleSize (0.04/0.4);
        g->GetXaxis ()->SetLabelSize (0.04/0.4);
        g->GetYaxis ()->SetLabelSize (0.04/0.4);
        g->GetXaxis ()->SetTitleOffset (2.5*0.4);
        g->GetYaxis ()->SetTitleOffset (1.5*0.4);

        g->GetYaxis ()->CenterTitle ();
        g->Draw ("AP");

        TLine* l = new TLine (-5.0, 1, 5.0, 1);
        l->SetLineColor (46);
        l->SetLineWidth (2);
        l->SetLineStyle (5);
        l->Draw ("same");
      }
      else {
        cout << "Warning in FullAnalysis :: PlotZEtaDst: Z Eta spectra ratio not stored, needs to be calculated!" << endl;
      }

      c->SaveAs (Form ("%s/ZEtaDists/z%s_eta_iCent%i.pdf", plotPath.Data (), spc, iCent));
    }
  }
  
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots the PSEUDORAPIDITY distribution of reconstructed Z bosons
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZYMap (FullAnalysis* a) {

  if (a)
    a->PlotZYMap ();

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iCent = 0; iCent < numCentBins; iCent++) {

      const char* canvasName = Form ("c_z_y_%s_iCent%i", spc, iCent);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      TPad* uPad = nullptr, *dPad = nullptr;
      if (canvasExists) {
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName)); 
        uPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_uPad", canvasName)));
        dPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_dPad", canvasName)));
      }
      else {
        c = new TCanvas (canvasName, "", 800, 800);
        c->cd ();
        uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 0.4, 1.0, 1.0);
        dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, 0.4);
        uPad->SetBottomMargin (0);
        dPad->SetTopMargin (0);
        dPad->SetBottomMargin (0.25);
        uPad->Draw ();
        dPad->Draw ();
        gDirectory->Add (c);
        gDirectory->Add (uPad);
        gDirectory->Add (dPad);
      }

      uPad->cd ();
      TH1D* h = h_z_y[iCent][iSpc][nPtZBins];

      if (plotFill) {
        //h->SetFillColorAlpha (fillColors[iCent], fillAlpha);
        h->SetFillColorAlpha (kYellow-4, fillAlpha);
        h->SetLineColor (kBlack);
        h->SetMarkerSize (0);
        h->SetLineWidth (0);
        h->GetYaxis ()->SetRangeUser (0, 0.6);

        h->GetXaxis ()->SetTitle ("#eta");
        h->GetYaxis ()->SetTitle ("1/N_{Z} dN_{Z}/dy_{Z}");
        h->GetXaxis ()->SetTitleSize (0.04/0.6);
        h->GetYaxis ()->SetTitleSize (0.04/0.6);
        h->GetXaxis ()->SetLabelSize (0.04/0.6);
        h->GetYaxis ()->SetLabelSize (0.04/0.6);
        h->GetXaxis ()->SetTitleOffset (1.5*0.6);
        h->GetYaxis ()->SetTitleOffset (1.5*0.6);

        h->DrawCopy (!canvasExists ? "bar" : "bar same");
        h->SetLineWidth (1);
        h->Draw ("hist same");

        gPad->RedrawAxis ();
      }
      else {
        TGraphAsymmErrors* g = GetTGAE (h);
        ResetXErrors (g);
        //deltaize (g, 0.1*(-1.5+iCent));

        const int markerStyle = kFullCircle;
        g->SetMarkerStyle (markerStyle);
        //g->SetMarkerSize (1);
        g->SetLineWidth (2);
        //g->SetLineColor (colors[iCent]);
        //g->SetMarkerColor (colors[iCent]);
        g->GetYaxis ()->SetRangeUser (0, 0.6);

        g->GetXaxis ()->SetTitle ("y_{Z}");
        g->GetYaxis ()->SetTitle ("(1/N_{Z}) (dN_{Z}/dy_{Z})");
        g->GetXaxis ()->SetTitleSize (0.04/0.6);
        g->GetYaxis ()->SetTitleSize (0.04/0.6);
        g->GetXaxis ()->SetLabelSize (0.04/0.6);
        g->GetYaxis ()->SetLabelSize (0.04/0.6);
        g->GetXaxis ()->SetTitleOffset (1.5*0.6);
        g->GetYaxis ()->SetTitleOffset (1.5*0.6);
        g->Draw (!canvasExists ? "AP" : "P");
      }

      if (!a)
        continue;

      myText (0.22, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.045/0.6);

      const char* spcLabel = iSpc == 0 ? "#it{Z} #rightarrow #it{e}#it{e}" : (iSpc == 1 ? "#it{Z} #rightarrow #it{#mu}#it{#mu}" : "#it{Z} #rightarrow #it{l}#it{l}");
      myText (0.66, 0.86, kBlack, spcLabel, 0.04/0.6);
      myText (0.66, 0.76, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV" , zPtBins[2]), 0.04/0.6);
      if (iCent == 0) myText (0.22, 0.76, kBlack, Form ("#it{pp}, 5.02 TeV"), 0.04/0.6);
      else            myText (0.22, 0.76, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04/0.6);
      //myMarkerText (0.273, 0.69, kBlack, kFullCircle, "Data", 1.25, 0.04/0.6);
      myMarkerText (0.30, 0.69, kBlack, kFullCircle, "Data", 1.25, 0.04/0.6);
      //myOnlyBoxText (0.28, 0.60, 1.2, kYellow-4, kBlack, 1, "MC", 0.04/0.6, 1001, 1);
      myOnlyBoxText (0.30, 0.60, 1.2, kYellow-4, kBlack, 1, "MC", 0.04/0.6, 1001, 1);
      

      dPad->cd ();

      h = (TH1D*) h_z_y[iCent][iSpc][nPtZBins]->Clone (Form ("h_z_y_ratio_%s_iCent%i_iPtZ%i_%s", spc, iCent, nPtZBins, name.c_str ()));
      h->Divide (a->h_z_y[iCent][iSpc][nPtZBins]);
      h_z_y_ratio[iCent][iSpc][nPtZBins] = h;

      if (h) {
        TGraphAsymmErrors* g = GetTGAE (h);
        ResetXErrors (g);

        const int markerStyle = kFullCircle;
        g->SetMarkerStyle (markerStyle);
        //g->SetMarkerSize (1);
        g->SetLineWidth (2);
        //g->SetLineColor (colors[iCent]);
        //g->SetMarkerColor (colors[iCent]);
        g->GetYaxis ()->SetRangeUser (0.66, 1.34);
        //g->GetYaxis ()->SetRangeUser (0, 2);

        g->GetXaxis ()->SetTitle ("y_{Z}");
        g->GetYaxis ()->SetTitle ("Data / MC");
        g->GetXaxis ()->SetTitleSize (0.04/0.4);
        g->GetYaxis ()->SetTitleSize (0.04/0.4);
        g->GetXaxis ()->SetLabelSize (0.04/0.4);
        g->GetYaxis ()->SetLabelSize (0.04/0.4);
        g->GetXaxis ()->SetTitleOffset (2.5*0.4);
        g->GetYaxis ()->SetTitleOffset (1.5*0.4);

        g->GetYaxis ()->CenterTitle ();
        g->Draw ("AP");

        TLine* l = new TLine (-2.5, 1, 2.5, 1);
        l->SetLineColor (46);
        l->SetLineWidth (2);
        l->SetLineStyle (5);
        l->Draw ("same");
      }
      else {
        cout << "Warning in FullAnalysis :: PlotZEtaDst: Z Eta spectra ratio not stored, needs to be calculated!" << endl;
      }

      c->SaveAs (Form ("%s/ZYDists/z%s_y_iCent%i.pdf", plotPath.Data (), spc, iCent));
    }
  }
  
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots the rapidity distribution of reconstructed Z->ee vs Z->mumu bosons
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZYMapSpcComp (const short pPtZ, FullAnalysis* a) {
  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    if (pPtZ != -1 && pPtZ != iPtZ)
      continue;

    const char* canvasName = Form ("c_z_y_SpcComp_iPtZ%i", iPtZ);
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    TPad* uPad = nullptr, *dPad = nullptr;
    if (canvasExists) {
      c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName)); 
      uPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_uPad", canvasName)));
      dPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_dPad", canvasName)));
    }
    else {
      c = new TCanvas (canvasName, "", 800, 800);
      c->cd ();
      uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 0.4, 1.0, 1.0);
      dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, 0.4);
      uPad->SetBottomMargin (0);
      dPad->SetTopMargin (0);
      dPad->SetBottomMargin (0.25);
      uPad->Draw ();
      dPad->Draw ();
      gDirectory->Add (c);
      gDirectory->Add (uPad);
      gDirectory->Add (dPad);
      c->cd ();
    }

    const short iCent = 0;

    for (short iSpc = 0; iSpc < 2; iSpc++) {
      //const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      TH1D* h = (TH1D*) h_z_y[iCent][iSpc][iPtZ]->Clone ();
      //h->Rebin (2);
      //h->Scale (0.5);

      uPad->cd ();
      if (plotFill) {
        h->SetFillColorAlpha (fillColors[iSpc+1], 0.2);
        h->SetLineColor (kBlack);
        h->SetMarkerSize (0);
        h->SetLineWidth (0);
        //h->GetYaxis ()->SetRangeUser (0, 1.3);
        h->GetYaxis ()->SetRangeUser (0, 0.5);

        h->GetXaxis ()->SetTitle ("y_{Z}");
        h->GetYaxis ()->SetTitle ("(1/N_{Z}) (dN_{Z}/dy_{Z})");
        h->GetXaxis ()->SetTitleSize (0.04/0.6);
        h->GetYaxis ()->SetTitleSize (0.04/0.6);
        h->GetXaxis ()->SetLabelSize (0.04/0.6);
        h->GetYaxis ()->SetLabelSize (0.04/0.6);
        h->GetXaxis ()->SetTitleOffset (1.5*0.6);
        h->GetYaxis ()->SetTitleOffset (1.5*0.6);

        if (!useAltMarker) {
          h->DrawCopy (iSpc == 0 && !canvasExists ? "bar" : "bar same");
          h->SetLineWidth (1);
        }
        else {
          h->SetLineStyle (2);
          h->SetLineColor (colors[iSpc+1]);
          h->SetFillColorAlpha (fillColors[iSpc+1], 0);
          h->SetLineWidth (2);
        }
        h->Draw ("hist same");

        gPad->RedrawAxis ();
      }
      else {
        TGraphAsymmErrors* g = GetTGAE (h);
        ResetXErrors (g);
        //deltaize (g, 0.1*(-1.5+iCent));

        const int markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (1);
        g->SetLineWidth (1);
        g->SetLineColor (colors[iSpc+1]);
        g->SetMarkerColor (colors[iSpc+1]);
        g->GetYaxis ()->SetRangeUser (0, 0.5);

        g->GetXaxis ()->SetTitle ("y_{Z}");
        g->GetYaxis ()->SetTitle ("1/N_{Z} dN_{Z}/dy_{Z}");
        g->GetXaxis ()->SetTitleSize (0.04/0.6);
        g->GetYaxis ()->SetTitleSize (0.04/0.6);
        g->GetXaxis ()->SetLabelSize (0.04/0.6);
        g->GetYaxis ()->SetLabelSize (0.04/0.6);
        g->GetXaxis ()->SetTitleOffset (1.5*0.6);
        g->GetYaxis ()->SetTitleOffset (1.5*0.6);
        g->Draw (iSpc == 0 && !canvasExists ? "AP" : "P");
      }

      if (a) {
        dPad->cd ();
        //while (a->h_z_y[iCent][iSpc][iPtZ]->GetNbinsX () >= 2*h->GetNbinsX ()) {
        //  a->h_z_y[iCent][iSpc][iPtZ]->Rebin (2);
        //  a->h_z_y[iCent][iSpc][iPtZ]->Scale (0.5);
        //}
        h->Divide (a->h_z_y[iCent][iSpc][iPtZ]);

        dPad->cd ();
        if (h) {
          TGraphAsymmErrors* g = GetTGAE (h);
          ResetXErrors (g);

          const int markerStyle = (useAltMarker ? kOpenCircle : kFullCircle);
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (1);
          g->SetLineWidth (1);
          g->SetLineColor (colors[iSpc+1]);
          g->SetMarkerColor (colors[iSpc+1]);
          g->GetYaxis ()->SetRangeUser (0.76, 1.24);
          //g->GetYaxis ()->SetRangeUser (0, 2);

          g->GetXaxis ()->SetTitle ("y_{Z}");
          g->GetYaxis ()->SetTitle ("Data / MC Truth");
          g->GetXaxis ()->SetTitleSize (0.04/0.4);
          g->GetYaxis ()->SetTitleSize (0.04/0.4);
          g->GetXaxis ()->SetLabelSize (0.04/0.4);
          g->GetYaxis ()->SetLabelSize (0.04/0.4);
          g->GetXaxis ()->SetTitleOffset (2.5*0.4);
          g->GetYaxis ()->SetTitleOffset (1.5*0.4);

          g->GetYaxis ()->CenterTitle ();
          g->Draw (iSpc == 0 && !canvasExists ? "AP" : "P");

          if (iSpc == 0 && !canvasExists) {
            TLine* l = new TLine (-2.5, 1, 2.5, 1);
            l->SetLineColor (46);
            l->SetLineWidth (2);
            l->SetLineStyle (5);
            l->Draw ("same");
          }
        }
      }
      else {
        cout << "Warning in FullAnalysis :: PlotZYMapSpcComp: Z y spectra ratio not stored, needs to be calculated!" << endl;
      }
    }
      
    uPad->cd ();
    myText (0.22, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.045/0.6);
    if (iCent == 0) myText (0.22, 0.82, kBlack, Form ("#it{pp}, 5.02 TeV"), 0.04/0.6);
    else            myText (0.22, 0.82, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04/0.6);
    if (iPtZ == nPtZBins-1) myText (0.22, 0.74, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.04/0.6);
    else                    myText (0.22, 0.74, kBlack, Form ("%g < #it{p}_{T} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04/0.6);

    myText (0.58, 0.90, kBlack, "Data", 0.032/0.6);
    myText (0.66, 0.90, kBlack, "MC Reco.", 0.032/0.6);
    for (int iSpc = 0; iSpc < 2; iSpc++) {
      const char* spcLabel = iSpc == 0 ? "#it{Z} #rightarrow #it{e}#it{e}" : (iSpc == 1 ? "#it{Z} #rightarrow #it{#mu}#it{#mu}" : "#it{Z} #rightarrow #it{l}#it{l}");
      myMarkerTextNoLine (0.64, 0.823-iSpc*0.06, colors[iSpc+1], kFullCircle, "", 1.8, 0.04/0.6);
      //myBoxText (0.72, 0.823-iSpc*0.06, colors[iSpc+1], kOpenCircle, "", 1.5, 0.04/0.6);
      myOnlyBoxText (0.75, 0.823-iSpc*0.06, 1.2, fillColors[iSpc+1], kBlack, 1, "", 0.06, 1001, 1);
      myText (0.79, 0.82-iSpc*0.06, colors[iSpc+1], spcLabel, 0.04/0.6);
    }

    c->SaveAs (Form ("%s/ZYDists/z_y_iPtZ%i.pdf", plotPath.Data (), iPtZ));
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z mass spectra
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZMassSpectra (FullAnalysis* a) {

  if (a)
    a->PlotZMassSpectra ();

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    for (short iReg = 0; iReg < 2; iReg++) {
      for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
        const char* canvasName = Form ("c_z_m_%s_iCent%i_iReg%i", spc, iCent, iReg);
        const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
        TCanvas* c = nullptr;
        TPad* uPad = nullptr, *dPad = nullptr;
        if (canvasExists) {
          c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName)); 
          uPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_uPad", canvasName)));
          dPad = dynamic_cast<TPad*>(gDirectory->Get (Form ("%s_dPad", canvasName)));
        }
        else {
          c = new TCanvas (canvasName, "", 800, 800);
          c->cd ();
          uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 0.4, 1.0, 1.0);
          dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, 0.4);
          uPad->SetBottomMargin (0);
          dPad->SetTopMargin (0);
          dPad->SetBottomMargin (0.25);
          uPad->Draw ();
          dPad->Draw ();
          gDirectory->Add (c);
          gDirectory->Add (uPad);
          gDirectory->Add (dPad);
        }

        uPad->cd ();
        TH1D* h = h_z_m[iCent][iSpc][iReg];

        //TF1* fit = new TF1 (Form ("f_z_m_%s_iCent%i_iReg%i", spc, iCent, iReg), "gaus(0)", 76, 106);
        //h->Fit (fit, "RN0Q");
        //{
        //  float mu = fit->GetParameter (1);
        //  float sigma = fit->GetParameter (2);
        //  delete fit;
        //  fit = new TF1 (Form ("f_z_m_%s_iCent%i_iReg%i", spc, iCent, iReg), "gaus(0)", mu-1.0*sigma, mu+1.2*sigma);
        //  h->Fit (fit, "RN0Q");
        //}

        if (plotFill) {
          //h->SetFillColorAlpha (fillColors[iCent], fillAlpha);
          h->SetFillColorAlpha (kYellow-4, fillAlpha);
          h->SetLineColor (kBlack);
          h->SetMarkerSize (0);
          h->SetLineWidth (0);
          //h->GetYaxis ()->SetRangeUser (0, 1.3);
          h->GetYaxis ()->SetRangeUser (0, 0.12);

          h->GetXaxis ()->SetTitle (Form ("m_{#it{%s}} [GeV]", (iSpc == 0 ? "ee" : (iSpc == 1 ? "#mu#mu" : "ll"))));
          //h->GetYaxis ()->SetTitle ("Arb. Units");
          h->GetYaxis ()->SetTitle ("Counts / Total");
          h->GetXaxis ()->SetTitleSize (0.04/0.6);
          h->GetYaxis ()->SetTitleSize (0.04/0.6);
          h->GetXaxis ()->SetLabelSize (0.04/0.6);
          h->GetYaxis ()->SetLabelSize (0.04/0.6);
          h->GetXaxis ()->SetTitleOffset (1.5*0.6);
          h->GetYaxis ()->SetTitleOffset (1.5*0.6);

          h->DrawCopy (!canvasExists ? "bar" : "bar same");
          h->SetLineWidth (1);
          h->Draw ("hist same");

          gPad->RedrawAxis ();

          //fit->SetLineColor (kRed+1);
          //fit->Draw ("same");
        }
        else {
          TGraphAsymmErrors* g = GetTGAE (h);
          ResetXErrors (g);
          //deltaize (g, 0.1*(-1.5+iCent));

          const int markerStyle = kFullCircle;
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (1);
          g->SetLineWidth (2);
          g->SetLineColor (kBlack);
          g->SetMarkerColor (kBlack);
          //g->GetYaxis ()->SetRangeUser (0, 1.3);
          g->GetYaxis ()->SetRangeUser (0, 0.12);

          g->GetXaxis ()->SetTitle (Form ("m_{#it{%s}} [GeV]", (iSpc == 0 ? "ee" : (iSpc == 1 ? "#mu#mu" : "ll"))));
          //g->GetYaxis ()->SetTitle ("Arb. Units");
          g->GetYaxis ()->SetTitle ("Counts / Total");
          g->GetXaxis ()->SetTitleSize (0.04/0.6);
          g->GetYaxis ()->SetTitleSize (0.04/0.6);
          g->GetXaxis ()->SetLabelSize (0.04/0.6);
          g->GetYaxis ()->SetLabelSize (0.04/0.6);
          g->GetXaxis ()->SetTitleOffset (1.5*0.6);
          g->GetYaxis ()->SetTitleOffset (1.5*0.6);
          g->Draw (!canvasExists ? "AP" : "P");

          //fit->SetLineColor (kBlack);
          //fit->Draw ("same");
        }

        if (!a)
          continue;

        LabelZMassSpectra (iSpc, iCent, iReg);
        

        dPad->cd ();

        h = (TH1D*) h_z_m[iCent][iSpc][iReg]->Clone (Form ("h_z_m_ratio_%s_iCent%i_iReg%i_%s", spc, iCent, iReg, name.c_str ()));
        h->Divide (a->h_z_m[iCent][iSpc][iReg]);
        h_z_m_ratio[iCent][iSpc][iReg] = h;

        if (h) {
          TGraphAsymmErrors* g = GetTGAE (h);
          ResetXErrors (g);

          const int markerStyle = kFullCircle;
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (1);
          g->SetLineWidth (2);
          g->SetLineColor (kBlack);
          g->SetMarkerColor (kBlack);
          g->GetYaxis ()->SetRangeUser (0.3, 1.7);

          g->GetXaxis ()->SetTitle (Form ("m_{#it{%s}} [GeV]", (iSpc == 0 ? "ee" : (iSpc == 1 ? "#mu#mu" : "ll"))));
          g->GetYaxis ()->SetTitle ("Data / MC");
          g->GetXaxis ()->SetTitleSize (0.04/0.4);
          g->GetYaxis ()->SetTitleSize (0.04/0.4);
          g->GetXaxis ()->SetLabelSize (0.04/0.4);
          g->GetYaxis ()->SetLabelSize (0.04/0.4);
          g->GetXaxis ()->SetTitleOffset (2.5*0.4);
          g->GetYaxis ()->SetTitleOffset (1.5*0.4);
          g->GetYaxis ()->CenterTitle ();
          g->Draw ("AP");

          TLine* l = new TLine (76, 1, 106, 1);
          l->SetLineColor (46);
          l->SetLineWidth (2);
          l->SetLineStyle (5);
          l->Draw ("same");
        }
        else {
          cout << "Warning in FullAnalysis :: PlotZMassSpectra: Z mass spectra ratio not stored, needs to be calculated!" << endl;
        }

        if (iReg == 2) 
          c->SaveAs (Form ("%s/ZMassSpectra/z%s_mass_spectrum_iCent%i.pdf", plotPath.Data (), spc, iCent));
        else
          c->SaveAs (Form ("%s/ZMassSpectra/z%s_mass_spectrum_iCent%i_iReg%i.pdf", plotPath.Data (), spc, iCent, iReg));
      }
    }
  }
}



void FullAnalysis :: LabelZMassSpectra (const short iSpc, const short iCent, const short iReg) {
  myText (0.22, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.045/0.6);
  const char* spc = iSpc == 0 ? "#it{Z} #rightarrow #it{e}#it{e}" : (iSpc == 1 ? "#it{Z} #rightarrow #it{#mu}#it{#mu}" : "#it{Z} #rightarrow #it{l}#it{l}");
  myText (0.71, 0.85, kBlack, spc, 0.04/0.6);
  if (iCent == 0) myText (0.22, 0.76, kBlack, Form ("#it{pp}, 5.02 TeV"), 0.04/0.6);
  else            myText (0.22, 0.76, kBlack, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04/0.6);

  //myMarkerText (0.753, 0.76, kBlack, kFullCircle, "Data", 1.25, 0.04/0.6);
  myMarkerText (0.78, 0.76, kBlack, kFullCircle, "Data", 1.25, 0.04/0.6);
  //myOnlyBoxText (0.76, 0.67, 1.2, kYellow-4, kBlack, 1, "MC", 0.04/0.6, 1001, 1);
  myOnlyBoxText (0.78, 0.67, 1.2, kYellow-4, kBlack, 1, "MC", 0.04/0.6, 1001, 1);

  if (iReg == 0)      myText (0.22, 0.67, kBlack, Form ("#left|y^{#it{%s}}#right| < 1", (iSpc == 0 ? "ee" : (iSpc == 1 ? "#mu#mu" : "ll"))), 0.04/0.6);
  else if (iReg == 1) myText (0.22, 0.67, kBlack, Form ("#left|y^{#it{%s}}#right| > 1", (iSpc == 0 ? "ee" : (iSpc == 1 ? "#mu#mu" : "ll"))), 0.04/0.6);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z yield with respect to the event plane angle
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZPhiYield (const short pSpc) {
  const char* canvasName = Form ("c_z_phi_%s", pSpc == 0 ? "ee" : (pSpc == 1 ? "#mu#mu" : "comb"));
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 600, 600);
    gDirectory->Add (c);
  }
  c->cd ();

  for (short iCent = 1; iCent < numCentBins; iCent++) {
    TH1D* h = h_z_phi[iCent][pSpc];
    //h->Rebin (8);
    //h->Scale (1. / h->Integral (), "width");
    float v2 = 0, v2err = 0;
    float c = 0, cerr = 0;

    for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
      c += h->GetBinContent (ix) * h->GetBinWidth (ix);
      cerr += pow (h->GetBinError (ix) * h->GetBinWidth (ix), 2);

      v2 += h->GetBinContent (ix) * cos (h->GetBinCenter (ix)) * h->GetBinWidth (ix);
      v2err += pow (h->GetBinError (ix) * cos (h->GetBinCenter (ix)) * h->GetBinWidth (ix), 2);
    }
    cerr = sqrt (cerr);
    v2err = sqrt (v2err);
    v2 = v2 / (c);
    v2err = sqrt (pow (v2err / c, 2) + pow (v2 * cerr / (c*c), 2));
    c = c / pi;
    cerr = cerr / pi;

    TGraphAsymmErrors* g = GetTGAE (h);
    deltaize (g, (1.5-iCent)*0.02, false);

    g->GetXaxis ()->SetTitle ("2#left|#phi_{Z} - #Psi_{2}#right|");
    g->GetYaxis ()->SetTitle ("(1/N_{Z}) (dN/d#Delta#phi)");

    g->GetYaxis ()->SetRangeUser (0.16, 0.5);

    g->SetLineColor (colors[iCent]);
    g->SetMarkerColor (colors[iCent]);
    g->SetMarkerSize (1);
    if (iCent == 1)
      g->Draw ("AP");
    else
      g->Draw ("P");

    TF1* fit = new TF1 (Form ("fit_iCent%i", iCent), "[0]*(1+2*[1]*cos(x))", 0, pi);
    fit->FixParameter (0, c);
    fit->FixParameter (1, v2);
    fit->SetLineColor (colors[iCent]);
    fit->SetLineStyle (2);
    fit->SetLineWidth (2);
    fit->Draw ("same");

    myText (0.22, 0.38-0.055*iCent, colors[iCent], Form ("v_{2} = %s", FormatMeasurement (v2, v2err, 2)), 0.04);
    myText (0.66, 0.88-0.055*iCent, colors[iCent], Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.04);

  } // end loop over cents

  //myText (0.66, 0.88, colors[0], "#it{pp}", 0.04);
  myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.66, 0.88, kBlack, "Pb+Pb, 5.02 TeV", 0.04);
  myText (0.22, 0.81, kBlack, "#it{p}_{T}^{Z} > 5 GeV", 0.04);
  if (pSpc != 2)
    myText (0.25, 0.72, kBlack, Form ("#it{Z} #rightarrow %s", pSpc == 0 ? "#it{e}#it{e}" : "#it{#mu}#it{#mu}"), 0.04);

  c->SaveAs (Form ("%s/q2_Mixing/ZPhiYields_%s.pdf", plotPath.Data (), pSpc == 0 ? "ee" : (pSpc == 1 ? "mumu" : "comb")));
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Z yield with respect to the event plane angle
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotZLeptonDPhi () {
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    const char* canvasName = Form ("c_z_lepton_dphi_iCent%i", iCent);
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 600, 600);
      gDirectory->Add (c);
      c->cd ();
    }
    c->cd ();

    for (short iSpc = 0; iSpc < 3; iSpc++) {
      //const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      TH1D* h = h_z_lepton_dphi[iCent][iSpc];

      h->SetLineColor (colors[iSpc]);
      h->SetMarkerColor (colors[iSpc]);

      h->GetXaxis ()->SetTitle ("#Delta#phi");
      h->GetYaxis ()->SetTitle ("Counts");

      h->Draw (canvasExists || iSpc != 0 ? "e1 same" : "e1");
    }
    c->SaveAs (Form ("%s/ZLeptonDPhi_iCent%i.pdf", plotPath.Data (), iCent));
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot distribution of Y(pT or xZh) for each event, binned in Z Pt
////////////////////////////////////////////////////////////////////////////////////////////////
void FullAnalysis :: PlotAllYields_Scatter_dPtZ (const bool useTrkPt, const short pSpc) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

    const char* canvasName = Form ("c_TrkYield_zpt_%s", useTrkPt ? "pTch" : "xhZ");
    const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
    TCanvas* c = nullptr;
    if (canvasExists)
      c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
    else {
      c = new TCanvas (canvasName, "", 400*numCentBins, 400);
      gDirectory->Add (c);
      c->Divide (4, 1);
    }

    for (short iCent = 0; iCent < numCentBins; iCent++) {
      c->cd (iCent+1);

      gPad->SetLogx ();
      gPad->SetLogy ();

      TH1D* h = new TH1D ("", "", useTrkPt ? nPtchBins[nPtZBins-1] : nXhZBins[nPtZBins-1], useTrkPt ? pTchBins[nPtZBins-1] : xhZBins[nPtZBins-1]);
      h->GetXaxis ()->SetRangeUser (useTrkPt ? pTchBins[nPtZBins-1][0] : xhZBins[nPtZBins-1][0], useTrkPt ? pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]] : xhZBins[nPtZBins-1][nXhZBins[nPtZBins-1]]);
      useTrkPt ? h->GetYaxis ()->SetRangeUser (8e-1, 1e3) : h->GetYaxis ()->SetRangeUser (8e-1, 1e3);

      //h->GetXaxis ()->SetMoreLogLabels ();

      h->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
      h->GetYaxis ()->SetTitle (useTrkPt ? "Y (#it{p}_{T}^{ ch})" : "Y (#it{x}_{hZ})");

      h->Draw ();

      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        TGraph* g = (useTrkPt ? g_trk_pt_ptz : g_trk_xhz_ptz)[iSpc][iPtZ][iCent];

        double x, y;
        for (int i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          g->SetPoint (i, (1.+0.07*(iPtZ-3))*x, y);
        }

        g->SetMarkerStyle (kOpenCircle);
        g->SetMarkerColor (colors[iPtZ-1]);
        g->SetMarkerSize (0.3);
        g->Draw ("P");
        
      } // end loop over iPtZ

      if (iCent == 0) {
        myText (0.22, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
        myText (0.22, 0.80, kBlack, "#it{pp}, 5.02 TeV", 0.06);
        myText (0.22, 0.74, kBlack, Form ("#it{Z} #rightarrow %s", iSpc == 0 ? "#it{e}#it{e}" : (iSpc == 1 ? "#it{#mu}#it{#mu}" : "#it{l}#it{l}")), 0.06);
      }
      else
        myText (0.55, 0.87, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);
    } // end loop over iCent

    c->SaveAs (Form ("%s/TrkYields/ScatterPlots/%s_dists_%s.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ", spc));
  } // end loop over iSpc
}



#endif
