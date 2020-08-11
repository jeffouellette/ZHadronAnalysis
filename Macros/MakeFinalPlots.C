#ifndef __MakeFinalPlots_C__
#define __MakeFinalPlots_C__

#include "Params.h"
//#include "PlotHybridModel.h"

#include <Utilities.h>
#include <ArrayTemplates.h>

#include <TFile.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TColor.h>
#include <TLatex.h>
#include <TPave.h>
#include <TMarker.h>

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

typedef TGraphErrors TGE;
typedef TGraphAsymmErrors TGAE;

const Color_t starColor     = (Color_t) tcolor->GetColor (50, 175, 50); //(Color_t) tcolor->GetColor (240, 197,  40);
const Style_t starMarker    = kFullCrossX;
//const Color_t cmsColor      = (Color_t) tcolor->GetColor (153,  55,  55);
const Color_t cmsColor      = (Color_t) tcolor->GetColor ( 71,  97, 130);
const Style_t cmsMarker     = kFullCircle;
const Color_t phenixColor   = kMagenta-6;
const Style_t phenixMarker  = kFullSquare;
const Color_t atlasColor    = finalColors[3];
const Color_t atlasFillColor = finalFillColors[3];

const Color_t hybridColor = (Color_t) tcolor->GetColor (82, 82, 173);
const double  hybridAlpha = 0.5;
const Color_t scetgColor  = (Color_t) tcolor->GetColor (39, 180, 66);
const double  scetgAlpha  = 0.6;
const Color_t jewelColor  = tcolor->GetColor (255, 170, 50);
const double  jewelAlpha  = 0.5;
const Color_t colbtColor  = kAzure+1;
const double  colbtAlpha  = 0.5;

const double minModelUnc = 0.06;


void SetMinErrors (TGAE* g, const double minerr, const bool logy) {
  double x, y;
  for (int i = 0; i < g->GetN (); i++) {
    g->GetPoint (i, x, y);
    if (y > 0 && g->GetErrorYhigh (i) / pow (y, logy ? 1 : 0) < minerr) g->SetPointEYhigh (i, minerr * pow (y, logy ? 1 : 0));
    if (y > 0 && g->GetErrorYlow (i) / pow (y, logy ? 1 : 0) < minerr) g->SetPointEYlow (i, minerr * pow (y, logy ? 1 : 0));
  }
}


//void MakeTheoryBox (const double x, const double y, const Color_t color, const double colorAlpha, const double boxMultiplier = 1.) {
//  const double ytsize = 0.07;
//  const double xtsize = 0.18;
//  const double y1 = y - 0.25*ytsize;
//  const double y2 = y + 0.25*ytsize;
//  const double x2 = x - 0.15*xtsize;
//  const double x1 = x - 0.55*xtsize*boxMultiplier;
//  TPave *mbox = new TPave (x1, y1, x2, y2, 0, "NDC");
//  mbox->SetFillColorAlpha (color, colorAlpha);
//  mbox->SetFillStyle (1001);
//  mbox->Draw ();
//
//  TLine mline;
//  mline.SetLineWidth (1);
//  mline.SetLineColor (color);
//  //mline.SetLineStyle (lstyle);
//  mline.SetLineStyle (0);
//  Double_t y_new = (y1+y2)/2.;
//  //mline.DrawLineNDC (x1, y_new, x2, y_new);
//  mline.DrawLineNDC (x1, y1, x2, y1);
//  mline.DrawLineNDC (x1, y2, x2, y2);
//  mline.DrawLineNDC (x1, y1, x1, y2);
//  mline.DrawLineNDC (x2, y1, x2, y2);
//  return;
//}
//
//
//void MakeDataBox (const double x, const double y, const Color_t color, const double colorAlpha, const Style_t mstyle, const double msize) {
//  MakeTheoryBox (x, y, color, colorAlpha);
//
//  const double ytsize = 0.07;
//  const double xtsize = 0.18;
//
//  const double y1 = y - 0.25*ytsize;
//  const double y2 = y + 0.25*ytsize;
//  const double x2 = x - 0.15*xtsize;
//  const double x1 = x - 0.55*xtsize;
//
//  TLine* ml = new TLine ();
//  ml->SetNDC();
//  ml->SetLineColor (color);
//  ml->SetLineStyle (1);
//  ml->SetLineWidth (2);
//
//  ml->DrawLineNDC (0.9*x1+0.1*x2, 0.5*(y1+y2), 0.1*x1+0.9*x2, 0.5*(y1+y2));
//  ml->DrawLineNDC (0.5*(x1+x2), 0.9*y1+0.1*y2, 0.5*(x1+x2), 0.1*y1+0.9*y2);
//
//  TMarker* marker = new TMarker (x-0.35*0.18, y, 0);
//  marker->SetNDC();
//  marker->SetMarkerColor (color);
//  marker->SetMarkerStyle (mstyle);
//  marker->SetMarkerSize (msize);
//  marker->Draw ();
//
//  if (IsFullMarker (mstyle)) {
//    TMarker* marker2 = new TMarker (x-0.35*0.18, y, 0);
//    marker2->SetNDC();
//    marker2->SetMarkerColor (kBlack);
//    marker2->SetMarkerStyle (FullToOpenMarker (mstyle));
//    marker2->SetMarkerSize (msize);
//    marker2->Draw();
//  }
//  return;
//}

TH1D*** h_trk_pt_ptz_iaa_stat_conf = Get2DArray <TH1D*> (nPtZBins, numCentBins);
TH1D*** h_trk_xhz_ptz_iaa_stat_conf = Get2DArray <TH1D*> (nPtZBins, numCentBins);
//TH1D*** h_trk_pt_ptz_sub_stat_conf = Get2DArray <TH1D*> (nPtZBins, numCentBins);
//TH1D*** h_trk_xhz_ptz_sub_stat_conf = Get2DArray <TH1D*> (nPtZBins, numCentBins);


TGAE*** g_trk_pt_ptz_iaa_syst_conf = Get2DArray <TGAE*> (nPtZBins, numCentBins);
TGAE*** g_trk_xhz_ptz_iaa_syst_conf = Get2DArray <TGAE*> (nPtZBins, numCentBins);
//TGAE*** g_trk_pt_ptz_sub_syst_conf = Get2DArray <TGAE*> (nPtZBins, numCentBins);
//TGAE*** g_trk_xhz_ptz_sub_syst_conf = Get2DArray <TGAE*> (nPtZBins, numCentBins);


TH1D*** h_trk_dphi_ptz_gt4_sub_stat = Get2DArray <TH1D*> (nPtZBins, numCentBins+1);
TH1D*** h_trk_dphi_ptz_lt4_sub_stat = Get2DArray <TH1D*> (nPtZBins, numCentBins+1);
TH1D*** h_trk_pt_ptz_sub_stat = Get2DArray <TH1D*> (nPtZBins, numCentBins);
TH1D*** h_trk_xhz_ptz_sub_stat = Get2DArray <TH1D*> (nPtZBins, numCentBins);
TH1D*** h_trk_pt_ptz_iaa_stat = Get2DArray <TH1D*> (nPtZBins, numCentBins);
//TH1D*** h_trk_pt_ptz_iaa_stat_2015proj = Get2DArray <TH1D*> (nPtZBins, numCentBins);
TH1D*** h_trk_xhz_ptz_iaa_stat = Get2DArray <TH1D*> (nPtZBins, numCentBins);
//TH1D*** h_trk_xhz_ptz_iaa_stat_2015proj = Get2DArray <TH1D*> (nPtZBins, numCentBins);
TGAE** g_trk_avg_pt_ptz_stat = Get1DArray <TGAE*> (numCentBins);
TGAE** g_trk_avg_xhz_ptz_stat = Get1DArray <TGAE*> (numCentBins);
TGAE** g_avg_ntrk_ptz_stat = Get1DArray <TGAE*> (numCentBins);


TGAE*** g_trk_dphi_ptz_lt4_sub_syst = Get2DArray <TGAE*> (nPtZBins, numCentBins+1);
TGAE*** g_trk_dphi_ptz_gt4_sub_syst = Get2DArray <TGAE*> (nPtZBins, numCentBins+1);
TGAE*** g_trk_pt_ptz_sub_syst = Get2DArray <TGAE*> (nPtZBins, numCentBins);
TGAE*** g_trk_xhz_ptz_sub_syst = Get2DArray <TGAE*> (nPtZBins, numCentBins);
TGAE*** g_trk_pt_ptz_iaa_syst = Get2DArray <TGAE*> (nPtZBins, numCentBins);
TGAE*** g_trk_xhz_ptz_iaa_syst = Get2DArray <TGAE*> (nPtZBins, numCentBins);
TGAE** g_trk_avg_pt_ptz_syst = Get1DArray <TGAE*> (numCentBins);
TGAE** g_trk_avg_xhz_ptz_syst = Get1DArray <TGAE*> (numCentBins);
TGAE** g_avg_ntrk_ptz_syst = Get1DArray <TGAE*> (numCentBins);


TGAE* tg_CMS_Zh_IAA_0_30_syst = nullptr;
TGAE* tg_CMS_Zh_IAA_0_30_stat = nullptr;


TGAE** g_pythia_pth_ptz = Get1DArray <TGAE*> (nPtZBins);
TGAE** g_pythia_finepth_ptz = Get1DArray <TGAE*> (nPtZBins);
TGAE** g_pythia_xhz_ptz = Get1DArray <TGAE*> (nPtZBins);
TGAE** g_pythia_finexhz_ptz = Get1DArray <TGAE*> (nPtZBins);


TGAE** g_colbt_pth = Get1DArray <TGAE*> (nPtZBins);
TGAE** g_colbt_xhz = Get1DArray <TGAE*> (nPtZBins);


TGAE** g_hybrid_xhz = Get1DArray <TGAE*> (nPtZBins);
TGAE** g_hybrid_pth  = Get1DArray <TGAE*> (nPtZBins);


TGAE** g_scetg_pth = Get1DArray <TGAE*> (nPtZBins);
TGAE** g_scetg_xhz = Get1DArray <TGAE*> (nPtZBins);


TGAE** g_jewel_xhz = Get1DArray <TGAE*> (nPtZBins);
TGAE** g_jewel_pth = Get1DArray <TGAE*> (nPtZBins);


TGE* tg_PHENIX_IAA_stat = nullptr;
TGE* tg_PHENIX_IAA_syst = nullptr;
TGAE* tg_STAR_IAA_stat = nullptr;
TGE* tg_STAR_IAA_syst = nullptr;
TGE* tg_CMS_IAA_stat = nullptr;
TGE* tg_CMS_IAA_syst = nullptr;



void MakeFinalPlots () {

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Load our histograms and systematics graphs from CONF note
  ////////////////////////////////////////////////////////////////////////////////////////////////
  TFile* CONFResultsFile = new TFile ("../rootFiles/Results/CONFNote/finalResults.root", "read");
  TFile* CONFSysFile = new TFile ("../rootFiles/Results/CONFNote/CombinedSys.root", "read");

  for (short iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
    for (short iCent = 1; iCent < numCentBins; iCent++) {
      h_trk_pt_ptz_iaa_stat_conf[iPtZ][iCent] = (TH1D*) CONFResultsFile->Get (Form ("h_ztrk_pt_iaa_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
      h_trk_xhz_ptz_iaa_stat_conf[iPtZ][iCent] = (TH1D*) CONFResultsFile->Get (Form ("h_ztrk_xhz_iaa_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
    }
    //for (short iCent = 0; iCent < numCentBins; iCent++) {
    //  h_trk_pt_ptz_sub_stat_conf[iPtZ][iCent] = (TH1D*) CONFResultsFile->Get (Form ("h_ztrk_pt_sub_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
    //  h_trk_xhz_ptz_sub_stat_conf[iPtZ][iCent] = (TH1D*) CONFResultsFile->Get (Form ("h_ztrk_xhz_sub_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
    //}
  }

  for (short iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
    for (short iCent = 1; iCent < numCentBins; iCent++) {
      g_trk_pt_ptz_iaa_syst_conf[iPtZ][iCent] = (TGAE*) CONFSysFile->Get (Form ("g_z_trk_zpt_iaa_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
      g_trk_xhz_ptz_iaa_syst_conf[iPtZ][iCent] = (TGAE*) CONFSysFile->Get (Form ("g_z_trk_zxzh_iaa_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
    }
    //for (short iCent = 0; iCent < numCentBins; iCent++) {
    //  g_trk_pt_ptz_sub_syst_conf[iPtZ][iCent] = (TGAE*) CONFSysFile->Get (Form ("g_z_trk_zpt_sub_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
    //  g_trk_xhz_ptz_sub_syst_conf[iPtZ][iCent] = (TGAE*) CONFSysFile->Get (Form ("g_z_trk_zxzh_sub_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
    //}
  }



  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Load our histograms and systematics graphs
  ////////////////////////////////////////////////////////////////////////////////////////////////
  TFile* resultsFile = new TFile ("../rootFiles/Results/finalResults.root", "read");
  TFile* sysFile = new TFile ("../rootFiles/Systematics/CombinedSys.root", "read");

  const double abbrev_xhz_bins[4] = {1./8., 1./4., 1./2., 1.};
  const short num_abbrev_xhz_bins = 3;

  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    for (short iCent : {0, numCentBins}) {
      h_trk_dphi_ptz_lt4_sub_stat[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_trk_dphi_sub_lt4_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
      h_trk_dphi_ptz_gt4_sub_stat[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_trk_dphi_sub_gt4_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
    } // end loop over iCent
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      h_trk_pt_ptz_sub_stat[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_trk_pt_ptz_sub_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
      h_trk_xhz_ptz_sub_stat[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_trk_xhz_ptz_sub_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
      if (iPtZ == 2) {
        TH1D* h = new TH1D (Form ("h_trk_xhz_ptz_sub_comb_iPtZ%i_iCent%i_data18_abbrev", iPtZ, iCent), "", num_abbrev_xhz_bins, abbrev_xhz_bins);
        for (int iX = 1; iX <= num_abbrev_xhz_bins; iX++) {
          h->SetBinContent (iX, h_trk_xhz_ptz_sub_stat[iPtZ][iCent]->GetBinContent (h_trk_xhz_ptz_sub_stat[iPtZ][iCent]->FindFixBin (h->GetBinCenter (iX))));
          h->SetBinError (iX, h_trk_xhz_ptz_sub_stat[iPtZ][iCent]->GetBinError (h_trk_xhz_ptz_sub_stat[iPtZ][iCent]->FindFixBin (h->GetBinCenter (iX))));
        }
        h_trk_xhz_ptz_sub_stat[iPtZ][iCent] = h;
      }
    } // end loop over iCent
    for (short iCent = 1; iCent < numCentBins; iCent++) {
      h_trk_pt_ptz_iaa_stat[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_trk_pt_ptz_iaa_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
      //h_trk_pt_ptz_iaa_stat_2015proj[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_trk_pt_ptz_iaa_2015proj_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
      h_trk_xhz_ptz_iaa_stat[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_trk_xhz_ptz_iaa_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));
      //h_trk_xhz_ptz_iaa_stat_2015proj[iPtZ][iCent] = (TH1D*) resultsFile->Get (Form ("h_trk_xhz_ptz_iaa_2015proj_comb_iPtZ%i_iCent%i_data18", iPtZ, iCent));

      if (iPtZ == 2) {
        TH1D* h = new TH1D (Form ("h_trk_xhz_ptz_iaa_comb_iPtZ%i_iCent%i_data18_abbrev", iPtZ, iCent), "", num_abbrev_xhz_bins, abbrev_xhz_bins);
        for (int iX = 1; iX <= num_abbrev_xhz_bins; iX++) {
          h->SetBinContent (iX, h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]->GetBinContent (h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]->FindFixBin (h->GetBinCenter (iX))));
          h->SetBinError (iX, h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]->GetBinError (h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]->FindFixBin (h->GetBinCenter (iX))));
        }
        h_trk_xhz_ptz_iaa_stat[iPtZ][iCent] = h;
      }
    } // end loop over iCent
  } // end loop over iPtZ
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    g_trk_avg_pt_ptz_stat[iCent] = (TGAE*) resultsFile->Get (Form ("g_trk_avg_pt_ptz_comb_iCent%i_data18", iCent));
    g_trk_avg_xhz_ptz_stat[iCent] = (TGAE*) resultsFile->Get (Form ("g_trk_avg_xhz_ptz_comb_iCent%i_data18", iCent));
    g_avg_ntrk_ptz_stat[iCent] = (TGAE*) resultsFile->Get (Form ("g_avg_ntrk_ptz_comb_iCent%i_data18", iCent));
  } // end loop over iCent

  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    for (short iCent : {0, numCentBins}) {
      g_trk_dphi_ptz_lt4_sub_syst[iPtZ][iCent] = (TGAE*) sysFile->Get (Form ("g_trk_dphi_sub_lt4_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
      g_trk_dphi_ptz_gt4_sub_syst[iPtZ][iCent] = (TGAE*) sysFile->Get (Form ("g_trk_dphi_sub_gt4_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
    } // end iCent = numCentBins scope
    for (short iCent = 1; iCent < numCentBins; iCent++) {
      g_trk_pt_ptz_iaa_syst[iPtZ][iCent] = (TGAE*) sysFile->Get (Form ("g_trk_pt_ptz_iaa_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
      g_trk_xhz_ptz_iaa_syst[iPtZ][iCent] = (TGAE*) sysFile->Get (Form ("g_trk_xhz_ptz_iaa_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
      if (iPtZ == 2) {
        g_trk_xhz_ptz_iaa_syst[iPtZ][iCent]->RemovePoint (0);
      }
    } // end loop over iCent
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      g_trk_pt_ptz_sub_syst[iPtZ][iCent] = (TGAE*) sysFile->Get (Form ("g_trk_pt_ptz_sub_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
      g_trk_xhz_ptz_sub_syst[iPtZ][iCent] = (TGAE*) sysFile->Get (Form ("g_trk_xhz_ptz_sub_comb_iPtZ%i_iCent%i_combSys", iPtZ, iCent));
      if (iPtZ == 2) {
        g_trk_xhz_ptz_sub_syst[iPtZ][iCent]->RemovePoint (0);
      }
    } // end loop over iCent
  } // end loop over iPtZ
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    g_trk_avg_pt_ptz_syst[iCent] = (TGAE*) sysFile->Get (Form ("g_trk_avg_pt_ptz_comb_iCent%i_combSys", iCent));
    g_trk_avg_xhz_ptz_syst[iCent] = (TGAE*) sysFile->Get (Form ("g_trk_avg_xhz_ptz_comb_iCent%i_combSys", iCent));
    g_avg_ntrk_ptz_syst[iCent] = (TGAE*) sysFile->Get (Form ("g_avg_ntrk_ptz_comb_iCent%i_combSys", iCent));
  } // end loop over iCent


  tg_CMS_Zh_IAA_0_30_syst = new TGAE ();
  tg_CMS_Zh_IAA_0_30_stat = new TGAE ();

  {
    const string fileName = "../DataThieving/CMS_Zh_IAA_HP2020/0_30centrality.dat";
    ifstream f;
    f.open (fileName.c_str ());

    double xlo, x, xhi, y, ystathi, ystatlo, ysysthi, ysystlo;
    f >> xlo >> x >> xhi >> y >> ystathi >> ystatlo >> ysysthi >> ysystlo;
    do {
      tg_CMS_Zh_IAA_0_30_syst->SetPoint (tg_CMS_Zh_IAA_0_30_syst->GetN (), x, y);
      tg_CMS_Zh_IAA_0_30_stat->SetPoint (tg_CMS_Zh_IAA_0_30_stat->GetN (), x, y);

      tg_CMS_Zh_IAA_0_30_syst->SetPointEXlow (tg_CMS_Zh_IAA_0_30_syst->GetN () - 1, fabs (x-xlo));
      tg_CMS_Zh_IAA_0_30_syst->SetPointEXhigh (tg_CMS_Zh_IAA_0_30_syst->GetN () - 1, fabs (xhi-x));
      tg_CMS_Zh_IAA_0_30_syst->SetPointEYlow (tg_CMS_Zh_IAA_0_30_syst->GetN () - 1, fabs (y-ysystlo));
      tg_CMS_Zh_IAA_0_30_syst->SetPointEYhigh (tg_CMS_Zh_IAA_0_30_syst->GetN () - 1, fabs (ysysthi-y));
      tg_CMS_Zh_IAA_0_30_stat->SetPointEYlow (tg_CMS_Zh_IAA_0_30_stat->GetN () - 1, fabs (y-ystatlo));
      tg_CMS_Zh_IAA_0_30_stat->SetPointEYhigh (tg_CMS_Zh_IAA_0_30_stat->GetN () - 1, fabs (ystathi-y));
      f >> xlo >> x >> xhi >> y >> ystathi >> ystatlo >> ysysthi >> ysystlo;
    } while (f);
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Get Powheg+Pythia hadron yield predictions
  ////////////////////////////////////////////////////////////////////////////////////////////////
  {
    TFile* pythiaFile = new TFile ("../rootFiles/TruthAnalysis/Nominal/pythiaCompare.root", "read");
    TH1D* h_z_counts = (TH1D*) pythiaFile->Get ("h_z_counts");
    TH1D* h = nullptr, *hSig = nullptr;
    TH2D* h2 = nullptr;
    TGAE* g = nullptr;
    double x, y;

    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      const float nTotalEvents = h_z_counts->GetBinContent (iPtZ+1);

      h = (TH1D*) pythiaFile->Get (Form ("h_trk_pTch_yield_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_pTch_cov_iPtZ%i", iPtZ));

      h->Scale (1./nTotalEvents, "width");
      h2->Scale (1., "width");
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        for (short iY = 1; iY <= h2->GetNbinsY (); iY++)
          h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (nTotalEvents)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
      h2->Scale (1./(nTotalEvents*nTotalEvents - nTotalEvents));
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        h->SetBinError (iX, sqrt (h2->GetBinContent (iX)));
      h->Scale (4./pi);

      hSig = h;

      h = (TH1D*) pythiaFile->Get (Form ("h_trk_pTch_yield_bkg_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_pTch_cov_bkg_iPtZ%i", iPtZ));

      h->Scale (1./nTotalEvents, "width");
      h2->Scale (1., "width");
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        for (short iY = 1; iY <= h2->GetNbinsY (); iY++)
          h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (nTotalEvents)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
      h2->Scale (1./(nTotalEvents*nTotalEvents - nTotalEvents));
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        h->SetBinError (iX, sqrt (h2->GetBinContent (iX)));
      h->Scale (8./pi);

      hSig->Add (h, -1);
      g_pythia_pth_ptz[iPtZ] = make_graph (hSig);



      h = (TH1D*) pythiaFile->Get (Form ("h_trk_FinePtch_yield_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_FinePtch_cov_iPtZ%i", iPtZ));

      h->Scale (1./nTotalEvents, "width");
      h2->Scale (1., "width");
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        for (short iY = 1; iY <= h2->GetNbinsY (); iY++)
          h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (nTotalEvents)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
      h2->Scale (1./(nTotalEvents*nTotalEvents - nTotalEvents));
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        h->SetBinError (iX, sqrt (h2->GetBinContent (iX)));
      h->Scale (4./pi);

      hSig = h;

      h = (TH1D*) pythiaFile->Get (Form ("h_trk_FinePtch_yield_bkg_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_FinePtch_cov_bkg_iPtZ%i", iPtZ));

      h->Scale (1./nTotalEvents, "width");
      h2->Scale (1., "width");
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        for (short iY = 1; iY <= h2->GetNbinsY (); iY++)
          h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (nTotalEvents)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
      h2->Scale (1./(nTotalEvents*nTotalEvents - nTotalEvents));
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        h->SetBinError (iX, sqrt (h2->GetBinContent (iX)));
      h->Scale (8./pi);

      hSig->Add (h, -1);
      g_pythia_finepth_ptz[iPtZ] = make_graph (hSig);

      g = g_pythia_finepth_ptz[iPtZ];

      
      h = (TH1D*) pythiaFile->Get (Form ("h_trk_xhZ_yield_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_xhZ_cov_iPtZ%i", iPtZ));

      h->Scale (1./nTotalEvents, "width");
      h2->Scale (1., "width");
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        for (short iY = 1; iY <= h2->GetNbinsY (); iY++)
          h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (nTotalEvents)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
      h2->Scale (1./(nTotalEvents*nTotalEvents - nTotalEvents));
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        h->SetBinError (iX, sqrt (h2->GetBinContent (iX)));
      h->Scale (4./pi);

      hSig = h;

      h = (TH1D*) pythiaFile->Get (Form ("h_trk_xhZ_yield_bkg_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_xhZ_cov_bkg_iPtZ%i", iPtZ));

      h->Scale (1./nTotalEvents, "width");
      h2->Scale (1., "width");
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        for (short iY = 1; iY <= h2->GetNbinsY (); iY++)
          h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (nTotalEvents)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
      h2->Scale (1./(nTotalEvents*nTotalEvents - nTotalEvents));
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        h->SetBinError (iX, sqrt (h2->GetBinContent (iX)));
      h->Scale (8./pi);

      hSig->Add (h, -1);
      g_pythia_xhz_ptz[iPtZ] = make_graph (hSig);



      h = (TH1D*) pythiaFile->Get (Form ("h_trk_FineXhZ_yield_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_FineXhZ_cov_iPtZ%i", iPtZ));

      h->Scale (1./nTotalEvents, "width");
      h2->Scale (1., "width");
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        for (short iY = 1; iY <= h2->GetNbinsY (); iY++)
          h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (nTotalEvents)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
      h2->Scale (1./(nTotalEvents*nTotalEvents - nTotalEvents));
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        h->SetBinError (iX, sqrt (h2->GetBinContent (iX)));
      h->Scale (4./pi);

      hSig = h;

      h = (TH1D*) pythiaFile->Get (Form ("h_trk_FineXhZ_yield_bkg_iPtZ%i", iPtZ));
      h2 = (TH2D*) pythiaFile->Get (Form ("h2_trk_FineXhZ_cov_bkg_iPtZ%i", iPtZ));

      h->Scale (1./nTotalEvents, "width");
      h2->Scale (1., "width");
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        for (short iY = 1; iY <= h2->GetNbinsY (); iY++)
          h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - (nTotalEvents)*(h->GetBinContent (iX))*(h->GetBinContent (iY)));
      h2->Scale (1./(nTotalEvents*nTotalEvents - nTotalEvents));
      for (short iX = 1; iX <= h2->GetNbinsX (); iX++)
        h->SetBinError (iX, sqrt (h2->GetBinContent (iX)));
      h->Scale (8./pi);

      hSig->Add (h, -1);
      g_pythia_finexhz_ptz[iPtZ] = make_graph (hSig);
    }
  }



  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Get CoLBT predictions
  ////////////////////////////////////////////////////////////////////////////////////////////////
  {
    for (int iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
      g_colbt_pth[iPtZ] = new TGAE ();
      g_colbt_xhz[iPtZ] = new TGAE ();
    }

    string modelFileName, dummyLine;
    short iPtZ;
    ifstream f;
    float x = 0, withResp = 0, statUnc = 0;

    modelFileName = "../CoLBT/Iaa_pt_zpt30_cent0_10.dat";
    iPtZ = 3;
    f.open (modelFileName.c_str ());
    getline (f, dummyLine);
    for (int i = 0; i < 5; i++) {
      f >> x >> withResp >> statUnc;
      g_colbt_pth[iPtZ]->SetPoint (g_colbt_pth[iPtZ]->GetN (), x, withResp);
      const double xUnc = g_trk_pt_ptz_iaa_syst[iPtZ][3]->GetErrorX (i);
      g_colbt_pth[iPtZ]->SetPointError (g_colbt_pth[iPtZ]->GetN ()-1, xUnc, xUnc, statUnc, statUnc);
    }
    f.close ();

    modelFileName = "../CoLBT/Iaa_pt_zpt60_cent0_10.dat";
    iPtZ = 4;
    f.open (modelFileName.c_str ());
    getline (f, dummyLine);
    for (int i = 0; i < 6; i++) {
      f >> x >> withResp >> statUnc;
      g_colbt_pth[iPtZ]->SetPoint (g_colbt_pth[iPtZ]->GetN (), x, withResp);
      const double xUnc = g_trk_pt_ptz_iaa_syst[iPtZ][3]->GetErrorX (i);
      g_colbt_pth[iPtZ]->SetPointError (g_colbt_pth[iPtZ]->GetN ()-1, xUnc, xUnc, statUnc, statUnc);
    }
    f.close ();

    modelFileName = "../CoLBT/Iaa_zt_zpt30_cent0_10.dat";
    iPtZ = 3;
    f.open (modelFileName.c_str ());
    getline (f, dummyLine);
    for (int i = 0; i < 5; i++) {
      f >> x >> withResp >> statUnc;
      g_colbt_xhz[iPtZ]->SetPoint (g_colbt_xhz[iPtZ]->GetN (), x, withResp);
      const double xUnc = g_trk_xhz_ptz_iaa_syst[iPtZ][3]->GetErrorX (i);
      g_colbt_xhz[iPtZ]->SetPointError (g_colbt_xhz[iPtZ]->GetN ()-1, xUnc, xUnc, statUnc, statUnc);
    }
    f.close ();

    modelFileName = "../CoLBT/Iaa_zt_zpt60_cent0_10.dat";
    iPtZ = 4;
    f.open (modelFileName.c_str ());
    getline (f, dummyLine);
    for (int i = 0; i < 6; i++) {
      f >> x >> withResp >> statUnc;
      g_colbt_xhz[iPtZ]->SetPoint (g_colbt_xhz[iPtZ]->GetN (), x, withResp);
      const double xUnc = g_trk_xhz_ptz_iaa_syst[iPtZ][3]->GetErrorX (i);
      g_colbt_xhz[iPtZ]->SetPointError (g_colbt_xhz[iPtZ]->GetN ()-1, xUnc, xUnc, statUnc, statUnc);
    }
    f.close ();
  }



  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Get Hybrid model predictions
  ////////////////////////////////////////////////////////////////////////////////////////////////
  {
    for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      g_hybrid_xhz[iPtZ] = new TGAE ();
      g_hybrid_pth[iPtZ] = new TGAE ();
    }

    const bool useTrkPt = false;
    //string modelFileName = "../HybridModel/IAAs/010_IAA_z_wake_1_ignore_neg_0.dat"; // means medium response including only the positive contribution from the wake
    string modelFileName = "../HybridModel/IAAs/010_IAA_z_wake_1_ignore_neg_1.dat"; // means full medium response, including also the negative contribution from the wake
    ifstream f;
    f.open (modelFileName.c_str ());
    float dummy = 0, x = 0, y1 = 0, y2 = 0, y = 0, yerr = 0;
    short ix = 0;
    while (f) {
      f >> x;
      for (int iPtZ = 2; iPtZ <= 4; iPtZ++) {
        f >> dummy >> dummy >> dummy >> dummy >> y1 >> y2;
        if (x*zPtBins[iPtZ+1] < 1)
          continue;
        y = 0.5*(y1+y2);
        yerr = 0.5*(y1-y2); 
        g_hybrid_xhz[iPtZ]->SetPoint (ix, x, y);
        g_hybrid_xhz[iPtZ]->SetPointEYhigh (ix, yerr);
        g_hybrid_xhz[iPtZ]->SetPointEYlow (ix, yerr);
      }
      ix++;
    }
    f.close ();

    //modelFileName = "../HybridModel/IAAs/010_IAA_pt_wake_1_ignore_neg_0.dat";
    modelFileName = "../HybridModel/IAAs/010_IAA_pt_wake_1_ignore_neg_1.dat";
    f.open (modelFileName.c_str ());
    ix = 0;
    while (f) {
      f >> x;
      for (int iPtZ = 2; iPtZ <= 4; iPtZ++) {
        f >> dummy >> dummy >> dummy >> dummy >> y1 >> y2;
        if (x > zPtBins[iPtZ])
          continue;
        y = 0.5*(y1+y2);
        yerr = 0.5*(y1-y2); 
        g_hybrid_pth[iPtZ]->SetPoint (ix, x, y);
        g_hybrid_pth[iPtZ]->SetPointEYhigh (ix, yerr);
        g_hybrid_pth[iPtZ]->SetPointEYlow (ix, yerr);
      }
      ix++;
    }
    f.close ();
  }



  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Get SCET_G predictions
  ////////////////////////////////////////////////////////////////////////////////////////////////
  {
    for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      g_scetg_pth[iPtZ] = new TGAE ();
      g_scetg_xhz[iPtZ] = new TGAE ();
    }

    ifstream f;
    double dummy = 0, x = 0, y = 0;

    TGAE* g;
    vector<double> xarr (0);
    vector<double> yarr_g1_8 (0), yarr_g2_0 (0), yarr_g2_2 (0); 
    string modelFileName;

    g = g_scetg_pth[2];
    xarr.clear ();
    yarr_g1_8.clear ();
    yarr_g2_0.clear ();
    yarr_g2_2.clear ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.ATLASptlowG1.8Z15-30";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      xarr.push_back (x);
      yarr_g1_8.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.ATLASptlowG2.0Z15-30";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_0.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.ATLASptlowG2.2Z15-30";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_2.push_back (y);
    }
    f.close ();

    for (int ix = 0; ix < xarr.size (); ix++) {
      double yarr[] = {yarr_g1_8[ix], yarr_g2_0[ix], yarr_g2_2[ix]};
      sort (yarr, yarr+sizeof(yarr)/sizeof(yarr[0]));
      g->SetPoint (ix, xarr[ix], yarr[1]);
      g->SetPointEYhigh (ix, yarr[2]-yarr[1]);
      g->SetPointEYlow  (ix, yarr[1]-yarr[0]);
    }



    g = g_scetg_xhz[2];
    xarr.clear ();
    yarr_g1_8.clear ();
    yarr_g2_0.clear ();
    yarr_g2_2.clear ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.ATLASlowG1.8Z15-30";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      xarr.push_back (x);
      yarr_g1_8.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.ATLASlowG2.0Z15-30";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_0.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.ATLASlowG2.2Z15-30";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_2.push_back (y);
    }
    f.close ();

    for (int ix = 0; ix < xarr.size (); ix++) {
      double yarr[] = {yarr_g1_8[ix], yarr_g2_0[ix], yarr_g2_2[ix]};
      sort (yarr, yarr+sizeof(yarr)/sizeof(yarr[0]));
      g->SetPoint (ix, xarr[ix], yarr[1]);
      g->SetPointEYhigh (ix, yarr[2]-yarr[1]);
      g->SetPointEYlow  (ix, yarr[1]-yarr[0]);
    }



    g = g_scetg_pth[3];
    xarr.clear ();
    yarr_g1_8.clear ();
    yarr_g2_0.clear ();
    yarr_g2_2.clear ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g1.8ATLASpt.Z30-60";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      xarr.push_back (x);
      yarr_g1_8.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g2.0ATLASpt.Z30-60";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_0.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g2.2ATLASpt.Z30-60";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_2.push_back (y);
    }
    f.close ();

    for (int ix = 0; ix < xarr.size (); ix++) {
      double yarr[] = {yarr_g1_8[ix], yarr_g2_0[ix], yarr_g2_2[ix]};
      sort (yarr, yarr+sizeof(yarr)/sizeof(yarr[0]));
      g->SetPoint (ix, xarr[ix], yarr[1]);
      g->SetPointEYhigh (ix, yarr[2]-yarr[1]);
      g->SetPointEYlow  (ix, yarr[1]-yarr[0]);
    }



    g = g_scetg_xhz[3];
    xarr.clear ();
    yarr_g1_8.clear ();
    yarr_g2_0.clear ();
    yarr_g2_2.clear ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g1.8ATLASnew.Z30-60";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      xarr.push_back (x);
      yarr_g1_8.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g2.0ATLASnew.Z30-60";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_0.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g2.2ATLASnew.Z30-60";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_2.push_back (y);
    }
    f.close ();

    for (int ix = 0; ix < xarr.size (); ix++) {
      double yarr[] = {yarr_g1_8[ix], yarr_g2_0[ix], yarr_g2_2[ix]};
      sort (yarr, yarr+sizeof(yarr)/sizeof(yarr[0]));
      g->SetPoint (ix, xarr[ix], yarr[1]);
      g->SetPointEYhigh (ix, yarr[2]-yarr[1]);
      g->SetPointEYlow  (ix, yarr[1]-yarr[0]);
    }



    g = g_scetg_pth[4];
    xarr.clear ();
    yarr_g1_8.clear ();
    yarr_g2_0.clear ();
    yarr_g2_2.clear ();
    
    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g1.8ATLASpt.Z60-100";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      xarr.push_back (x);
      yarr_g1_8.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g2.0ATLASpt.Z60-100";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_0.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g2.2ATLASpt.Z60-100";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_2.push_back (y);
    }
    f.close ();
    for (int ix = 0; ix < xarr.size (); ix++) {
      double yarr[] = {yarr_g1_8[ix], yarr_g2_0[ix], yarr_g2_2[ix]};
      sort (yarr, yarr+sizeof(yarr)/sizeof(yarr[0]));
      g->SetPoint (ix, xarr[ix], yarr[1]);
      g->SetPointEYhigh (ix, yarr[2]-yarr[1]);
      g->SetPointEYlow  (ix, yarr[1]-yarr[0]);
    }



    g = g_scetg_xhz[4];
    xarr.clear ();
    yarr_g1_8.clear ();
    yarr_g2_0.clear ();
    yarr_g2_2.clear ();
    
    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g1.8ATLASnew.Z60-100";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      xarr.push_back (x);
      yarr_g1_8.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g2.0ATLASnew.Z60-100";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_0.push_back (y);
    }
    f.close ();

    modelFileName = "../SCETg_HP2020/R_SigZ0Had5020_010.g2.2ATLASnew.Z60-100";
    f.open (modelFileName.c_str ());
    while (f) {
      f >> x >> y;
      yarr_g2_2.push_back (y);
    }
    f.close ();
    for (int ix = 0; ix < xarr.size (); ix++) {
      double yarr[] = {yarr_g1_8[ix], yarr_g2_0[ix], yarr_g2_2[ix]};
      sort (yarr, yarr+sizeof(yarr)/sizeof(yarr[0]));
      g->SetPoint (ix, xarr[ix], yarr[1]);
      g->SetPointEYhigh (ix, yarr[2]-yarr[1]);
      g->SetPointEYlow  (ix, yarr[1]-yarr[0]);
    }

  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Get JEWEL predictions
  ////////////////////////////////////////////////////////////////////////////////////////////////
  {
    TFile* jewelFile = new TFile ("../rootFiles/Jewel/hists.root", "read");
    for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      TH1D* h_num = (TH1D*) jewelFile->Get (Form ("h_z_trk_pth_medium_iSpc2_iPtZ%i", iPtZ));
      TH1D* h_den = (TH1D*) jewelFile->Get (Form ("h_z_trk_pth_vacuum_iSpc2_iPtZ%i", iPtZ));
      h_num->Divide (h_den);
      g_jewel_pth[iPtZ] = make_graph (h_num);

      h_num = (TH1D*) jewelFile->Get (Form ("h_z_trk_xhz_medium_iSpc2_iPtZ%i", iPtZ));
      h_den = (TH1D*) jewelFile->Get (Form ("h_z_trk_xhz_vacuum_iSpc2_iPtZ%i", iPtZ));
      h_num->Divide (h_den);
      g_jewel_xhz[iPtZ] = make_graph (h_num);
    }
    jewelFile->Close ();
  }



  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Get data comparison plots
  ////////////////////////////////////////////////////////////////////////////////////////////////
  TFile* dataCompFile = new TFile ("../DataComparisons/Zhadron_IAA_data_comparisons.root", "read");
  tg_PHENIX_IAA_stat = (TGE*) dataCompFile->Get ("tg_PHENIX_IAA_stat");
  tg_PHENIX_IAA_syst = (TGE*) dataCompFile->Get ("tg_PHENIX_IAA_syst");
  tg_STAR_IAA_stat = (TGAE*) dataCompFile->Get ("tg_STAR_IAA_stat");
  tg_STAR_IAA_syst = (TGE*) dataCompFile->Get ("tg_STAR_IAA_syst");
  tg_CMS_IAA_stat = (TGE*) dataCompFile->Get ("tg_CMS_IAA_stat");
  tg_CMS_IAA_syst = (TGE*) dataCompFile->Get ("tg_CMS_IAA_syst");

  TLatex* tl = new TLatex ();
  TLine* l = new TLine ();




  //{
  //  const char* canvasName = "c1";
  //  TCanvas* c1 = new TCanvas (canvasName, "", 1600, 1200);

  //  const double llMargin = 0.17;
  //  const double lrMargin = 0.032;
  //  const double clMargin = 0.032;
  //  const double crMargin = 0.032;
  //  const double rlMargin = 0.032;
  //  const double rrMargin = 0.040;
  //  const double bMargin = 0.15;
  //  const double tMargin = 0.04;

  //  const double deltaL = (1. - llMargin - lrMargin);
  //  const double deltaC = (1. - clMargin - crMargin);
  //  const double deltaR = (1. - rlMargin - rrMargin);

  //  const double a = (double) (deltaR * deltaC / (deltaL*deltaR + deltaC*deltaR + deltaL*deltaC));
  //  const double b = (double) (deltaR * deltaL / (deltaL*deltaR + deltaC*deltaR + deltaL*deltaC));

  //  const double xPadLCMiddle = a;
  //  const double xPadCRMiddle = a+b;

  //  const double yPadMiddle = 0.5;

  //  TPad* luPad = nullptr;
  //  TPad* cuPad = nullptr;
  //  TPad* ruPad = nullptr;
  //  TPad* ldPad = nullptr;
  //  TPad* cdPad = nullptr;
  //  TPad* rdPad = nullptr;

  //  luPad = new TPad (Form ("%s_luPad", canvasName), "", 0, yPadMiddle, xPadLCMiddle, 1);
  //  cuPad = new TPad (Form ("%s_cuPad", canvasName), "", xPadLCMiddle, yPadMiddle, xPadCRMiddle, 1);
  //  ruPad = new TPad (Form ("%s_ruPad", canvasName), "", xPadCRMiddle, yPadMiddle, 1, 1);
  //  ldPad = new TPad (Form ("%s_ldPad", canvasName), "", 0, 0, xPadLCMiddle, yPadMiddle);
  //  cdPad = new TPad (Form ("%s_cdPad", canvasName), "", xPadLCMiddle, 0, xPadCRMiddle, yPadMiddle);
  //  rdPad = new TPad (Form ("%s_rdPad", canvasName), "", xPadCRMiddle, 0, 1, yPadMiddle);

  //  luPad->SetLeftMargin (llMargin);
  //  luPad->SetRightMargin (lrMargin);
  //  cuPad->SetLeftMargin (clMargin);
  //  cuPad->SetRightMargin (crMargin);
  //  ruPad->SetLeftMargin (rlMargin);
  //  ruPad->SetRightMargin (rrMargin);
  //  luPad->SetBottomMargin (bMargin);
  //  luPad->SetTopMargin (tMargin);
  //  cuPad->SetBottomMargin (bMargin);
  //  cuPad->SetTopMargin (tMargin);
  //  ruPad->SetBottomMargin (bMargin);
  //  ruPad->SetTopMargin (tMargin);
  //  ldPad->SetLeftMargin (llMargin);
  //  ldPad->SetRightMargin (lrMargin);
  //  cdPad->SetLeftMargin (clMargin);
  //  cdPad->SetRightMargin (crMargin);
  //  rdPad->SetLeftMargin (rlMargin);
  //  rdPad->SetRightMargin (rrMargin);
  //  ldPad->SetBottomMargin (bMargin);
  //  ldPad->SetTopMargin (tMargin);
  //  cdPad->SetBottomMargin (bMargin);
  //  cdPad->SetTopMargin (tMargin);
  //  rdPad->SetBottomMargin (bMargin);
  //  rdPad->SetTopMargin (tMargin);

  //  TPad* uPads[3] = {luPad, cuPad, ruPad};
  //  TPad* dPads[3] = {ldPad, cdPad, rdPad};

  //  for (int i = 0; i < 3; i++) {
  //    uPads[i]->Draw ();
  //    dPads[i]->Draw ();
  //  }
  //  for (int i = 0; i < 3; i++) {
  //    uPads[i]->SetLogx ();
  //    uPads[i]->SetLogy ();
  //    dPads[i]->SetLogx ();
  //    dPads[i]->SetLogy ();
  //  }


  //  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
  //    uPads[iPtZ-2]->cd ();

  //    const double xmin = pTchBins[iPtZ][0];
  //    const double xmax = pTchBins[iPtZ][nPtchBins[iPtZ]];

  //    TH1D* h = new TH1D ("", "", 1, xmin, xmax);

  //    TAxis* xax = h->GetXaxis ();
  //    TAxis* yax = h->GetYaxis ();

  //    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
  //    xax->SetRangeUser (xmin, xmax);
  //    xax->SetTitleFont (43);
  //    xax->SetTitleSize (30);
  //    xax->SetTitleOffset (2.5);
  //    xax->SetLabelSize (0);

  //    yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
  //    const double ymin = 0.10;
  //    const double ymax = 10;
  //    yax->SetRangeUser (ymin, ymax);
  //    yax->SetTitleFont (43);
  //    yax->SetTitleSize (30);
  //    yax->SetTitleOffset (3.2);
  //    yax->SetLabelSize (0);

  //    h->SetLineWidth (0);

  //    h->DrawCopy ("");
  //    SaferDelete (&h);

  //    tl->SetTextFont (43);
  //    tl->SetTextSize (28);
  //    tl->SetTextAlign (21);
  //    const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
  //    tl->DrawLatex (1,  yoff, "1");
  //    tl->DrawLatex (2,  yoff, "2");
  //    tl->DrawLatex (3,  yoff, "3");
  //    tl->DrawLatex (4,  yoff, "4");
  //    tl->DrawLatex (5,  yoff, "5");
  //    tl->DrawLatex (6,  yoff, "6");
  //    tl->DrawLatex (7,  yoff, "7");
  //    tl->DrawLatex (10, yoff, "10");
  //    tl->DrawLatex (20, yoff, "20");
  //    tl->DrawLatex (30, yoff, "30");
  //    tl->DrawLatex (40, yoff, "40");
  //    tl->DrawLatex (60, yoff, "60");

  //    if (iPtZ == 2) {
  //      tl->SetTextAlign (32);

  //      const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-llMargin-lrMargin));
  //      tl->DrawLatex (xoff, 0.1, "0.1");
  //      tl->DrawLatex (xoff, 0.2, "0.2");
  //      tl->DrawLatex (xoff, 0.3, "0.3");
  //      tl->DrawLatex (xoff, 0.5, "0.5");
  //      tl->DrawLatex (xoff, 0.7, "0.7");
  //      tl->DrawLatex (xoff, 1, "1");
  //      tl->DrawLatex (xoff, 2, "2");
  //      tl->DrawLatex (xoff, 3, "3");
  //      tl->DrawLatex (xoff, 5, "5");
  //      tl->DrawLatex (xoff, 7, "7");
  //      tl->DrawLatex (xoff, 10, "10");
  //    }

  //    l->SetLineStyle (2);
  //    l->SetLineWidth (2);
  //    l->SetLineColor (kBlack);
  //    l->DrawLine (xmin, 1, xmax, 1);
  //  } // end loop over iPtZ

  //  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
  //    uPads[iPtZ-2]->cd ();

  //    for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
  //      TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

  //      RecenterGraph (g_syst);
  //      ResetXErrors (g_syst);
  //      deltaize (g_syst, 0.09*(iCent-1 - 0.5*(numCentBins-2)), true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
  //      ResetXErrors (g_syst);
  //      SetConstantXErrors (g_syst, 0.040, true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);

  //      g_syst->SetMarkerSize (0);
  //      g_syst->SetLineWidth (1);
  //      g_syst->SetMarkerColor (finalColors[iCent]);
  //      g_syst->SetLineColor (finalColors[iCent]);
  //      g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

  //      g_syst->Draw ( "5P");

  //    } // end loop over iCent

  //    for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
  //      TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

  //      RecenterGraph (g_stat);
  //      ResetXErrors (g_stat);
  //      deltaize (g_stat, 0.09*(iCent-1 - 0.5*(numCentBins-2)), true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
  //      ResetXErrors (g_stat);

  //      Style_t markerStyle = markerStyles[iCent-1];
  //      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
  //      g_stat->SetMarkerStyle (markerStyle);
  //      g_stat->SetMarkerSize (markerSize);
  //      g_stat->SetLineWidth (3);
  //      g_stat->SetMarkerColor (finalColors[iCent]);
  //      g_stat->SetLineColor (finalColors[iCent]);

  //      ((TGAE*) g_stat->Clone ())->Draw ("P");

  //      markerStyle = FullToOpenMarker (markerStyle);
  //      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
  //      
  //      g_stat->SetMarkerStyle (markerStyle);
  //      g_stat->SetMarkerSize (markerSize);
  //      g_stat->SetLineWidth (0);
  //      g_stat->SetMarkerColor (kBlack);

  //      ((TGAE*) g_stat->Clone ())->Draw ("P");

  //      SaferDelete (&g_stat);
  //    } // end loop over iCent
  //  } // end loop over iPtZ


  //  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
  //    dPads[iPtZ-2]->cd ();

  //    const double xmin = (iPtZ == 2 ? 1./8. : xhZBins[iPtZ][0]);
  //    const double xmax = xhZBins[iPtZ][nXhZBins[iPtZ]];

  //    TH1D* h = new TH1D ("", "", 1, xmin, xmax);

  //    TAxis* xax = h->GetXaxis ();
  //    TAxis* yax = h->GetYaxis ();

  //    xax->SetTitle ("#it{x}_{h#it{Z}} = #it{p}_{T}^{ch} #/#it{p}_{T}^{#it{Z}}");

  //    xax->SetRangeUser (xmin, xmax);
  //    xax->SetTitleFont (43);
  //    xax->SetTitleSize (30);
  //    xax->SetTitleOffset (2.5);
  //    xax->SetLabelSize (0);

  //    yax->SetTitle ("#it{I}_{AA} (#it{x}_{h#it{Z}})");
  //    const double ymin = 0.05;
  //    const double ymax = 10;
  //    yax->SetRangeUser (ymin, ymax);
  //    yax->SetTitleFont (43);
  //    yax->SetTitleSize (30);
  //    yax->SetTitleOffset (3.2);
  //    yax->SetLabelSize (0);

  //    h->SetLineWidth (0);

  //    h->DrawCopy ("");
  //    SaferDelete (&h);

  //    tl->SetTextFont (43);
  //    tl->SetTextSize (28);
  //    tl->SetTextAlign (21);
  //    const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
  //    if (iPtZ > 2) {
  //      if (iPtZ > 3) tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
  //      tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
  //    }
  //    if (iPtZ >= 3) tl->DrawLatex (1e-1,  yoff, "10^{-1}");
  //    tl->DrawLatex (2e-1,  yoff, "2#times10^{-1}");
  //    tl->DrawLatex (5e-1,  yoff, "5#times10^{-1}");
  //    tl->DrawLatex (1,     yoff, "1");

  //    if (iPtZ == 2) {
  //      tl->SetTextAlign (32);
  //
  //      const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-llMargin-lrMargin));
  //      tl->DrawLatex (xoff, 0.05, "0.05");
  //      tl->DrawLatex (xoff, 0.1, "0.1");
  //      tl->DrawLatex (xoff, 0.2, "0.2");
  //      tl->DrawLatex (xoff, 0.3, "0.3");
  //      tl->DrawLatex (xoff, 0.5, "0.5");
  //      tl->DrawLatex (xoff, 1, "1");
  //      tl->DrawLatex (xoff, 2, "2");
  //      tl->DrawLatex (xoff, 3, "3");
  //      tl->DrawLatex (xoff, 5, "5");
  //      tl->DrawLatex (xoff, 10, "10");
  //    }

  //    l->SetLineStyle (2);
  //    l->SetLineWidth (2);
  //    l->DrawLine (xmin, 1, xmax, 1);
  //  } // end loop over iPtZ

  //  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
  //    dPads[iPtZ-2]->cd ();

  //    const double xmin = (iPtZ == 2 ? 1./8. : xhZBins[iPtZ][0]);
  //    const double xmax = xhZBins[iPtZ][nXhZBins[iPtZ]];

  //    for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
  //      TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_iaa_syst[iPtZ][iCent]->Clone ();

  //      RecenterGraph (g_syst);
  //      ResetXErrors (g_syst);
  //      deltaize (g_syst, 0.09*(iCent-1 - 0.5*(numCentBins-2)), true, xmin, xmax);
  //      ResetXErrors (g_syst);
  //      SetConstantXErrors (g_syst, 0.040, true, xmin, xmax);

  //      g_syst->SetMarkerSize (0);
  //      g_syst->SetLineWidth (1);
  //      g_syst->SetMarkerColor (finalColors[iCent]);
  //      g_syst->SetLineColor (finalColors[iCent]);
  //      g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

  //      g_syst->Draw ("5P");
  //    } // end loop over iCent

  //    for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
  //      TGAE* g_stat = make_graph (h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]);

  //      Style_t markerStyle = markerStyles[iCent-1];
  //      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

  //      RecenterGraph (g_stat);
  //      ResetXErrors (g_stat);
  //      deltaize (g_stat, 0.09*(iCent-1 - 0.5*(numCentBins-2)), true, xmin, xmax);
  //      ResetXErrors (g_stat);

  //      g_stat->SetMarkerStyle (markerStyle);
  //      g_stat->SetMarkerSize (markerSize);
  //      g_stat->SetLineWidth (3);
  //      g_stat->SetMarkerColor (finalColors[iCent]);
  //      g_stat->SetLineColor (finalColors[iCent]);

  //      ((TGAE*) g_stat->Clone ())->Draw ("P");

  //      markerStyle = FullToOpenMarker (markerStyle);
  //      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
  //      
  //      g_stat->SetMarkerStyle (markerStyle);
  //      g_stat->SetMarkerSize (markerSize);
  //      g_stat->SetLineWidth (0);
  //      g_stat->SetMarkerColor (kBlack);

  //      ((TGAE*) g_stat->Clone ())->Draw ("P");

  //      SaferDelete (&g_stat);
  //    } // end loop over iCent
  //  } // end loop over iPtZ
  //  

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextAlign (11);

  //  luPad->cd ();
  //  tl->SetTextSize (32);
  //  tl->DrawLatexNDC (0.21, 0.86, "#bf{#it{ATLAS}} Internal");
  //  tl->SetTextSize (28);
  //  tl->DrawLatexNDC (0.21, 0.79, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.21, 0.73, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

  //  ldPad->cd ();
  //  tl->SetTextSize (32);
  //  tl->DrawLatexNDC (0.21, 0.86, "#bf{#it{ATLAS}} Internal");
  //  tl->SetTextSize (28);
  //  tl->DrawLatexNDC (0.21, 0.79, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.21, 0.73, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

  //  tl->SetTextSize (30);
  //  luPad->cd ();
  //  tl->DrawLatexNDC (0.56, 0.86, "15 < #it{p}_{T}^{#it{Z}} < 30 GeV");
  //  cuPad->cd ();
  //  tl->DrawLatexNDC (0.50, 0.86, "30 < #it{p}_{T}^{#it{Z}} < 60 GeV");
  //  ruPad->cd ();
  //  tl->DrawLatexNDC (0.60, 0.86, "#it{p}_{T}^{#it{Z}} > 60 GeV");

  //  ldPad->cd ();
  //  tl->DrawLatexNDC (0.56, 0.86, "15 < #it{p}_{T}^{#it{Z}} < 30 GeV");
  //  cdPad->cd ();
  //  tl->DrawLatexNDC (0.50, 0.86, "30 < #it{p}_{T}^{#it{Z}} < 60 GeV");
  //  rdPad->cd ();
  //  tl->DrawLatexNDC (0.60, 0.86, "#it{p}_{T}^{#it{Z}} > 60 GeV");

  //  ruPad->cd ();
  //  myMarkerAndBoxAndLineText (0.52, 0.780, 3.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "30-80\% #/#it{pp}", 0.02 / (gPad->GetWNDC ()));
  //  myMarkerAndBoxAndLineText (0.52, 0.715, 3.0, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.8, "10-30\% #/#it{pp}", 0.02 / (gPad->GetWNDC ()));
  //  myMarkerAndBoxAndLineText (0.52, 0.650, 3.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "0-10\% #/#it{pp}", 0.02 / (gPad->GetWNDC ()));

  //  rdPad->cd ();
  //  myMarkerAndBoxAndLineText (0.52, 0.780, 3.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "30-80\% #/#it{pp}", 0.02 / (gPad->GetWNDC ()));
  //  myMarkerAndBoxAndLineText (0.52, 0.715, 3.0, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.8, "10-30\% #/#it{pp}", 0.02 / (gPad->GetWNDC ()));
  //  myMarkerAndBoxAndLineText (0.52, 0.650, 3.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "0-10\% #/#it{pp}", 0.02 / (gPad->GetWNDC ()));
  //  

  //  c1->SaveAs (Form ("%s/FinalPlots/iaa_allptz.pdf", plotPath.Data ()));
  //}




  //{
  //  TCanvas* c2 = new TCanvas ("c2", "", 800, 800);

  //  const double lMargin = 0.15;
  //  const double rMargin = 0.04;
  //  const double bMargin = 0.15;
  //  const double tMargin = 0.04;

  //  c2->SetLeftMargin (lMargin);
  //  c2->SetRightMargin (rMargin);
  //  c2->SetBottomMargin (bMargin);
  //  c2->SetTopMargin (tMargin);

  //  c2->SetLogx ();
  //  c2->SetLogy ();

  //  {
  //    TH1D* h = new TH1D ("", "", nPtchBins[nPtZBins-1], pTchBins[nPtZBins-1]);

  //    TAxis* xax = h->GetXaxis ();
  //    TAxis* yax = h->GetYaxis ();

  //    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
  //    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
  //    xax->SetRangeUser (pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
  //    //xax->SetTitleFont (43);
  //    //xax->SetTitleSize (30);
  //    //xax->SetTitleOffset (1.25);
  //    xax->SetLabelSize (0);

  //    yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
  //    const double ymin = 0.08;
  //    const double ymax = 1200;
  //    yax->SetRangeUser (ymin, ymax);
  //    //yax->SetTitleFont (43);
  //    //yax->SetTitleSize (30);
  //    //yax->SetTitleOffset (1.30);
  //    //yax->SetLabelFont (43);
  //    //yax->SetLabelSize (28);

  //    h->SetLineWidth (0);

  //    h->DrawCopy ("");
  //    SaferDelete (&h);

  //    tl->SetTextFont (43);
  //    tl->SetTextSize (36);
  //    tl->SetTextAlign (21);
  //    const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));

  //    tl->DrawLatex (1,  yoff, "1");
  //    tl->DrawLatex (2,  yoff, "2");
  //    tl->DrawLatex (3,  yoff, "3");
  //    tl->DrawLatex (4,  yoff, "4");
  //    tl->DrawLatex (5,  yoff, "5");
  //    tl->DrawLatex (6,  yoff, "6");
  //    tl->DrawLatex (7,  yoff, "7");
  //    tl->DrawLatex (10, yoff, "10");
  //    tl->DrawLatex (20, yoff, "20");
  //    tl->DrawLatex (30, yoff, "30");
  //    tl->DrawLatex (40, yoff, "40");
  //    tl->DrawLatex (60, yoff, "60");

  //    l->SetLineStyle (2);
  //    l->SetLineWidth (2);
  //    l->SetLineColor (kBlack);
  //    l->DrawLine (1, 1, 14, 1);
  //    l->DrawLine (1, 10, 28, 10);
  //    l->DrawLine (1, 100, 60, 100);
  //  }

  //  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
  //    for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
  //      TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

  //      OffsetYAxis (g_syst, pow (10, iPtZ-2), true);
  //      RecenterGraph (g_syst);
  //      ResetXErrors (g_syst);
  //      deltaize (g_syst, 0.06*(iCent-1 - 0.5*(numCentBins-2)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
  //      ResetXErrors (g_syst);
  //      SetConstantXErrors (g_syst, 0.040, true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);

  //      g_syst->SetMarkerSize (0);
  //      g_syst->SetLineWidth (1);
  //      g_syst->SetMarkerColor (finalColors[iCent]);
  //      g_syst->SetLineColor (finalColors[iCent]);
  //      g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

  //      g_syst->Draw ( "5P");

  //    } // end loop over iCent

  //    for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
  //      TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

  //      OffsetYAxis (g_stat, pow (10, iPtZ-2), true);
  //      RecenterGraph (g_stat);
  //      ResetXErrors (g_stat);
  //      deltaize (g_stat, 0.06*(iCent-1 - 0.5*(numCentBins-2)), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
  //      ResetXErrors (g_stat);

  //      Style_t markerStyle = markerStyles[iCent-1];
  //      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);
  //      g_stat->SetMarkerStyle (markerStyle);
  //      g_stat->SetMarkerSize (markerSize);
  //      g_stat->SetLineWidth (3);
  //      g_stat->SetMarkerColor (finalColors[iCent]);
  //      g_stat->SetLineColor (finalColors[iCent]);

  //      ((TGAE*) g_stat->Clone ())->Draw ("P");

  //      markerStyle = FullToOpenMarker (markerStyle);
  //      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);
  //      
  //      g_stat->SetMarkerStyle (markerStyle);
  //      g_stat->SetMarkerSize (markerSize);
  //      g_stat->SetLineWidth (0);
  //      g_stat->SetMarkerColor (kBlack);

  //      ((TGAE*) g_stat->Clone ())->Draw ("P");

  //      SaferDelete (&g_stat);
  //    } // end loop over iCent
  //  } // end loop over iPtZ

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextAlign (11);

  //  tl->SetTextSize (28);
  //  tl->DrawLatexNDC (0.19, 0.890, "#bf{#it{ATLAS}} Internal");
  //  tl->SetTextSize (26);
  //  tl->DrawLatexNDC (0.43, 0.890, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.43, 0.845, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

  //  tl->SetTextSize (24);
  //  tl->DrawLatex (3.40, 1.28, "15 < #it{p}_{T}^{#it{Z]} < 30 GeV (#times 1)");
  //  tl->DrawLatex (6.35, 12.8, "30 < #it{p}_{T}^{#it{Z]} < 60 GeV (#times 10)");
  //  tl->DrawLatex (13, 128, "#it{p}_{T}^{#it{Z}} > 60 GeV (#times 10^{2})");

  //  tl->SetTextSize (26);
  //  tl->SetTextAlign (11);
  //  tl->DrawLatexNDC (0.765, 0.320-0.011, "30-80\% #/#it{pp}");
  //  tl->DrawLatexNDC (0.765, 0.270-0.011, "10-30\% #/#it{pp}");
  //  tl->DrawLatexNDC (0.765, 0.220-0.011, "0-10\% #/#it{pp}");
  //  MakeDataBox   (0.78, 0.320, finalFillColors[1], 0.30, markerStyles[0], 1.6);
  //  MakeDataBox   (0.78, 0.270, finalFillColors[2], 0.30, markerStyles[1], 1.6);
  //  MakeDataBox   (0.78, 0.220, finalFillColors[3], 0.30, markerStyles[2], 2.2);
  //  //myMarkerAndBoxAndLineText (0.76, 0.280, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.6, "30-80\% #/#it{pp}", 0.032);
  //  //myMarkerAndBoxAndLineText (0.76, 0.230, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.6, "10-30\% #/#it{pp}", 0.032);
  //  //myMarkerAndBoxAndLineText (0.76, 0.180, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.2, "0-10\% #/#it{pp}", 0.032);

  //  c2->SaveAs (Form ("%s/FinalPlots/iaa_allptz_pTchOnly_onePlot.pdf", plotPath.Data ()));
  //}




  //{
  //  const char* canvasName = "c3";
  //  TCanvas* c3 = new TCanvas (canvasName, "", 800, 960);

  //  const double lMargin = 0.12;
  //  const double rMargin = 0.04;
  //  const double bMargin = 0.10;
  //  const double tMargin = 0.04;

  //  const double uPad_tMargin = 0.0659094;
  //  const double dPad_bMargin = 0.254384;
  //  const double y1 = 0.393106;
  //  const double y2 = 0.686212;

  //  TPad* uPad = new TPad (Form ("%s_uPad", canvasName), "", 0, y2, 1, 1);
  //  TPad* cPad = new TPad (Form ("%s_cPad", canvasName), "", 0, y1, 1, y2);
  //  TPad* dPad = new TPad (Form ("%s_dPad", canvasName), "", 0, 0, 1, y1);

  //  uPad->SetLeftMargin (lMargin);
  //  uPad->SetRightMargin (rMargin);
  //  cPad->SetLeftMargin (lMargin);
  //  cPad->SetRightMargin (rMargin);
  //  dPad->SetLeftMargin (lMargin);
  //  dPad->SetRightMargin (rMargin);

  //  uPad->SetTopMargin (uPad_tMargin);
  //  uPad->SetBottomMargin (0);
  //  cPad->SetTopMargin (0);
  //  cPad->SetBottomMargin (0);
  //  dPad->SetTopMargin (0);
  //  dPad->SetBottomMargin (dPad_bMargin);

  //  uPad->SetLogx ();
  //  uPad->SetLogy ();
  //  cPad->SetLogx ();
  //  cPad->SetLogy ();
  //  dPad->SetLogx ();
  //  dPad->SetLogy ();

  //  c3->cd ();
  //  uPad->Draw ();
  //  cPad->Draw ();
  //  dPad->Draw ();

  //  std::vector<TPad*> pads = {dPad, cPad, uPad};
  //  std::map<TPad*, double> ymaxvals;
  //  std::map<TPad*, double> yminvals;

  //  ymaxvals[uPad] = 4;
  //  yminvals[uPad] = 0.3;
  //  ymaxvals[cPad] = 3;
  //  yminvals[cPad] = 0.2;
  //  ymaxvals[dPad] = 2;
  //  yminvals[dPad] = 0.12;

  //  const double xmin = pTchBins[nPtZBins-1][0];
  //  const double xmax = pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]];

  //  for (TPad* pad : pads) {
  //    pad->cd ();

  //    TH1D* h = new TH1D ("", "", 1, xmin, xmax);

  //    TAxis* xax = h->GetXaxis ();
  //    TAxis* yax = h->GetYaxis ();

  //    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
  //    xax->SetRangeUser (xmin, xmax);
  //    xax->SetTitleFont (43);
  //    xax->SetTitleSize (30);
  //    xax->SetTitleOffset (3.4);
  //    xax->SetLabelSize (0);

  //    xax->SetTickLength (0.04 * (1.-gPad->GetTopMargin ()-gPad->GetBottomMargin ()) / gPad->GetHNDC ());
  //    xax->SetTickLength (0.04);

  //    yax->SetRangeUser (yminvals[pad], ymaxvals[pad]);
  //    if (pad == uPad) {
  //      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
  //      yax->SetTitleFont (43);
  //      yax->SetTitleSize (30);
  //      yax->SetTitleOffset (1.80);
  //    }
  //    else {
  //      yax->SetTitleSize (0);
  //    }
  //    yax->SetLabelSize (0);

  //    h->SetLineWidth (0);

  //    h->DrawCopy ("");
  //    SaferDelete (&h);
  //  }

  //  {
  //    dPad->cd ();

  //    tl->SetTextFont (43);
  //    tl->SetTextSize (32);
  //    tl->SetTextAlign (21);

  //    const double yoff = yminvals[dPad] / exp (0.086 * (log (ymaxvals[dPad]) - log (yminvals[dPad])) / (1.-0-dPad_bMargin));
  //    tl->DrawLatex (1,  yoff, "1");
  //    tl->DrawLatex (2,  yoff, "2");
  //    tl->DrawLatex (3,  yoff, "3");
  //    tl->DrawLatex (4,  yoff, "4");
  //    tl->DrawLatex (5,  yoff, "5");
  //    tl->DrawLatex (6,  yoff, "6");
  //    tl->DrawLatex (7,  yoff, "7");
  //    tl->DrawLatex (10, yoff, "10");
  //    tl->DrawLatex (20, yoff, "20");
  //    tl->DrawLatex (30, yoff, "30");
  //    tl->DrawLatex (40, yoff, "40");
  //    tl->DrawLatex (60, yoff, "60");
  //  }

  //  for (TPad* pad : pads) {
  //    pad->cd ();

  //    tl->SetTextAlign (32);

  //    const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-lMargin-rMargin));
  //    if (yminvals[pad] == 0.12) tl->DrawLatex (xoff, 0.2, "0.2");
  //    if (yminvals[pad] <= 0.2) tl->DrawLatex (xoff, 0.3, "0.3");
  //    if (yminvals[pad] == 0.3) tl->DrawLatex (xoff, 0.4, "0.4");
  //    tl->DrawLatex (xoff, 0.6, "0.6");
  //    tl->DrawLatex (xoff, 1, "1");
  //    if (ymaxvals[pad] > 2) tl->DrawLatex (xoff, 2, "2");
  //    if (ymaxvals[pad] > 3) tl->DrawLatex (xoff, 3, "3");

  //    l->SetLineStyle (2);
  //    l->SetLineWidth (2);
  //    l->SetLineColor (kBlack);
  //    l->DrawLine (xmin, 1, xmax, 1);
  //  }

  //  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
  //    pads[iPtZ-2]->cd ();

  //    const double xmin = pTchBins[nPtZBins-1][0];
  //    const double xmax = pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]];

  //    for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
  //      TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

  //      RecenterGraph (g_syst);
  //      ResetXErrors (g_syst);
  //      deltaize (g_syst, 0.06*(iCent-1 - 0.5*(numCentBins-2)), true, xmin, xmax);
  //      ResetXErrors (g_syst);
  //      SetConstantXErrors (g_syst, 0.040, true, xmin, xmax);

  //      g_syst->SetMarkerSize (0);
  //      g_syst->SetLineWidth (1);
  //      g_syst->SetMarkerColor (finalColors[iCent]);
  //      g_syst->SetLineColor (finalColors[iCent]);
  //      g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

  //      g_syst->Draw ( "5P");

  //    } // end loop over iCent

  //    for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
  //      TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

  //      RecenterGraph (g_stat);
  //      ResetXErrors (g_stat);
  //      deltaize (g_stat, 0.06*(iCent-1 - 0.5*(numCentBins-2)), true, xmin, xmax);
  //      ResetXErrors (g_stat);

  //      Style_t markerStyle = markerStyles[iCent-1];
  //      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);
  //      g_stat->SetMarkerStyle (markerStyle);
  //      g_stat->SetMarkerSize (markerSize);
  //      g_stat->SetLineWidth (3);
  //      g_stat->SetMarkerColor (finalColors[iCent]);
  //      g_stat->SetLineColor (finalColors[iCent]);

  //      ((TGAE*) g_stat->Clone ())->Draw ("P");

  //      markerStyle = FullToOpenMarker (markerStyle);
  //      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);
  //      
  //      g_stat->SetMarkerStyle (markerStyle);
  //      g_stat->SetMarkerSize (markerSize);
  //      g_stat->SetLineWidth (0);
  //      g_stat->SetMarkerColor (kBlack);

  //      ((TGAE*) g_stat->Clone ())->Draw ("P");

  //      SaferDelete (&g_stat);
  //    } // end loop over iCent
  //  } // end loop over iPtZ

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextAlign (11);

  //  c3->cd ();
  //  tl->SetTextSize (30);
  //  tl->DrawLatexNDC (0.39, 0.930, "#bf{#it{ATLAS}} Internal");
  //  tl->SetTextSize (28);
  //  tl->DrawLatexNDC (0.41, 0.885, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.41, 0.845, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

  //  tl->SetTextSize (28);
  //  tl->DrawLatexNDC (0.745, 0.938, "#it{p}_{T}^{#it{Z}} > 60 GeV");
  //  tl->DrawLatexNDC (0.673, 0.646, "30 < #it{p}_{T}^{#it{Z}} < 60 GeV");
  //  tl->DrawLatexNDC (0.671, 0.350, "15 < #it{p}_{T}^{#it{Z}} < 30 GeV");

  //  c3->cd ();

  //  tl->SetTextSize (28);
  //  tl->SetTextAlign (11);
  //  tl->DrawLatexNDC (0.775, 0.260-0.009, "30-80\%^{}/^{}#it{pp}");
  //  tl->DrawLatexNDC (0.775, 0.210-0.009, "10-30\%^{}/^{}#it{pp}");
  //  tl->DrawLatexNDC (0.775, 0.160-0.009, "0-10\%^{}/^{}#it{pp}");
  //  MakeDataBox   (0.79, 0.260, finalFillColors[1], 0.30, markerStyles[0], 1.6);
  //  MakeDataBox   (0.79, 0.210, finalFillColors[2], 0.30, markerStyles[1], 1.6);
  //  MakeDataBox   (0.79, 0.160, finalFillColors[3], 0.30, markerStyles[2], 2.2);

  //  c3->SaveAs (Form ("%s/FinalPlots/iaa_allptz_pTchOnly_threePlots.pdf", plotPath.Data ()));
  //}




  //{
  //  TCanvas* c4 = new TCanvas ("c4", "", 800, 800);
  //  const double lMargin = 0.12;
  //  const double rMargin = 0.04;
  //  const double bMargin = 0.15;
  //  const double tMargin = 0.04;

  //  c4->SetLeftMargin (lMargin);
  //  c4->SetRightMargin (rMargin);
  //  c4->SetBottomMargin (bMargin);
  //  c4->SetTopMargin (tMargin);

  //  c4->SetLogx ();
  //  c4->SetLogy ();

  //  short iPtZ = nPtZBins-1;
  //  const short iCent = 3;

  //  const double xmin = (iPtZ == 2 ? 1./8. : xhZBins[iPtZ][0]);
  //  const double xmax = xhZBins[iPtZ][nXhZBins[iPtZ]];

  //  {
  //    TH1D* h = new TH1D ("", "", 1, xmin, xmax);

  //    TAxis* xax = h->GetXaxis ();
  //    TAxis* yax = h->GetYaxis ();

  //    xax->SetTitle ("#it{x}_{h,#gamma/#it{Z}} = #it{p}_{T}^{ch} #/#it{p}_{T}^{#gamma/#it{Z}}");
  //    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
  //    xax->SetRangeUser (xmin, xmax);
  //    xax->SetLabelSize (0);

  //    yax->SetTitle ("#it{I}_{AA} (#it{x}_{h#it{Z}})");
  //    yax->SetTitleOffset (0.8 * yax->GetTitleOffset ());
  //    yax->SetLabelSize (0);
  //    const double ymin = 0.05;
  //    const double ymax = 7;
  //    yax->SetRangeUser (ymin, ymax);

  //    h->SetLineWidth (0);

  //    h->DrawCopy ("");
  //    SaferDelete (&h);

  //    tl->SetTextFont (43);
  //    tl->SetTextSize (36);
  //    tl->SetTextAlign (21);

  //    const double yoff = ymin / exp (0.054 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
  //    if (iPtZ > 2) {
  //      if (iPtZ > 3) tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
  //      tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
  //    }
  //    tl->DrawLatex (1e-1,  yoff, "10^{-1}");
  //    tl->DrawLatex (2e-1,  yoff, "2#times10^{-1}");
  //    tl->DrawLatex (5e-1,  yoff, "5#times10^{-1}");
  //    tl->DrawLatex (1,     yoff, "1");

  //    tl->SetTextAlign (32);
  //
  //    const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-lMargin-lMargin));
  //    tl->DrawLatex (xoff, 0.1, "0.1");
  //    tl->DrawLatex (xoff, 0.2, "0.2");
  //    tl->DrawLatex (xoff, 0.3, "0.3");
  //    tl->DrawLatex (xoff, 0.5, "0.5");
  //    tl->DrawLatex (xoff, 1, "1");
  //    tl->DrawLatex (xoff, 2, "2");
  //    tl->DrawLatex (xoff, 3, "3");
  //    tl->DrawLatex (xoff, 5, "5");

  //    l->SetLineStyle (2);
  //    l->SetLineWidth (2);
  //    l->DrawLine (xmin, 1, xmax, 1);
  //  }

  //  TGE* g = nullptr;

  //  g = (TGE*) tg_CMS_IAA_syst->Clone ();
  //  g->SetMarkerSize (0);
  //  g->SetLineColor (cmsColor);
  //  g->SetLineWidth (1);
  //  g->SetFillColorAlpha (cmsColor, 0.3);
  //  ((TGE*) g->Clone ())->Draw ("5P");

  //  SaferDelete (&g);

  //  g = (TGE*) tg_CMS_IAA_stat->Clone ();
  //  g->SetMarkerStyle (cmsMarker);
  //  g->SetMarkerSize (1.8);
  //  g->SetMarkerColor (cmsColor);
  //  g->SetLineColor (cmsColor);
  //  g->SetLineWidth (3);
  //  ((TGE*) g->Clone ())->Draw ("P");

  //  if (IsFullMarker (cmsMarker)) {
  //    g->SetMarkerStyle (FullToOpenMarker (cmsMarker));
  //    g->SetMarkerSize (1.8);
  //    g->SetLineWidth (0);
  //    g->SetMarkerColor (kBlack);
  //    ((TGE*) g->Clone ())->Draw ("P");
  //  }

  //  SaferDelete (&g);

  //  TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_iaa_syst[iPtZ][iCent]->Clone ();
  //  RecenterGraph (g_syst);
  //  ResetXErrors (g_syst);
  //  SetConstantXErrors (g_syst, 0.060, true, xmin, xmax);

  //  g_syst->SetMarkerSize (0);
  //  g_syst->SetMarkerColor (atlasColor);
  //  g_syst->SetLineColor (atlasColor);
  //  g_syst->SetLineWidth (1);
  //  g_syst->SetFillColorAlpha (atlasFillColor, 0.3);

  //  g_syst->Draw ("5P");

  //  TGAE* g_stat = make_graph (h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]);

  //  RecenterGraph (g_stat);
  //  ResetXErrors (g_stat);

  //  g_stat->SetMarkerStyle (kFullDiamond);
  //  g_stat->SetMarkerSize (2.4);
  //  g_stat->SetLineWidth (3);
  //  g_stat->SetMarkerColor (atlasColor);
  //  g_stat->SetLineColor (atlasColor);

  //  ((TGAE*) g_stat->Clone ())->Draw ("P");

  //  g_stat->SetMarkerStyle (kOpenDiamond);
  //  g_stat->SetMarkerSize (2.4);
  //  g_stat->SetLineWidth (0);
  //  g_stat->SetMarkerColor (kBlack);

  //  ((TGAE*) g_stat->Clone ())->Draw ("P");

  //  SaferDelete (&g_stat);

  //  myMarkerAndBoxAndLineText (0.284, 0.888, 1.4, 1001, atlasFillColor, 0.30, atlasColor, kFullDiamond,  2.4, "", 0.032);
  //  myMarkerAndBoxAndLineText (0.284, 0.838, 1.4, 1001, cmsColor,       0.30, cmsColor,   cmsMarker, 1.8, "", 0.032);

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextSize (26);
  //  tl->SetTextAlign (11);
  //  tl->DrawLatexNDC (0.290, 0.886, "ATLAS #it{p}_{T}^{#it{Z}} > 60 GeV, 0-10% Pb+Pb");
  //  tl->DrawLatexNDC (0.290, 0.836, "CMS #it{p}_{T}^{#gamma} > 60 GeV, #it{p}_{T}^{jet} > 30 GeV, 0-10\% Pb+Pb");

  //  tl->SetTextSize (32);
  //  tl->DrawLatexNDC (0.24, 0.300, "#bf{#it{ATLAS}} Internal");
  //  tl->SetTextSize (30);
  //  tl->DrawLatexNDC (0.26, 0.255, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.26, 0.210, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

  //  c4->SaveAs (Form ("%s/FinalPlots/iaa_xhz_cmsComp.pdf", plotPath.Data ()));
  //}




  //{
  //  TCanvas* c5 = new TCanvas ("c5", "", 800, 800);
  //  const double lMargin = 0.12;
  //  const double rMargin = 0.04;
  //  const double bMargin = 0.15;
  //  const double tMargin = 0.04;

  //  c5->SetLeftMargin (lMargin);
  //  c5->SetRightMargin (rMargin);
  //  c5->SetBottomMargin (bMargin);
  //  c5->SetTopMargin (tMargin);

  //  c5->SetLogx ();
  //  c5->SetLogy ();

  //  short iPtZ = nPtZBins-1;
  //  const short iCent = 3;

  //  const double xmin = (iPtZ == 2 ? 1./8. : xhZBins[iPtZ][0]);
  //  const double xmax = xhZBins[iPtZ][nXhZBins[iPtZ]];

  //  {
  //    TH1D* h = new TH1D ("", "", 1, xmin, xmax);

  //    TAxis* xax = h->GetXaxis ();
  //    TAxis* yax = h->GetYaxis ();

  //    xax->SetTitle ("#it{x}_{h,#gamma/#it{Z}} = #it{p}_{T}^{ch} #/#it{p}_{T}^{#gamma/#it{Z}}");
  //    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
  //    xax->SetRangeUser (xmin, xmax);
  //    xax->SetLabelSize (0);

  //    yax->SetTitle ("#it{I}_{AA} (#it{x}_{h#it{Z}})");
  //    yax->SetTitleOffset (0.8 * yax->GetTitleOffset ());
  //    yax->SetLabelSize (0);
  //    const double ymin = 0.05;
  //    const double ymax = 7;
  //    yax->SetRangeUser (ymin, ymax);

  //    h->SetLineWidth (0);

  //    h->DrawCopy ("");
  //    SaferDelete (&h);

  //    tl->SetTextFont (43);
  //    tl->SetTextSize (36);
  //    tl->SetTextAlign (21);

  //    const double yoff = ymin / exp (0.054 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
  //    if (iPtZ > 2) {
  //      if (iPtZ > 3) tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
  //      tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
  //    }
  //    tl->DrawLatex (1e-1,  yoff, "10^{-1}");
  //    tl->DrawLatex (2e-1,  yoff, "2#times10^{-1}");
  //    tl->DrawLatex (5e-1,  yoff, "5#times10^{-1}");
  //    tl->DrawLatex (1,     yoff, "1");

  //    tl->SetTextAlign (32);
  //
  //    const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-lMargin-lMargin));
  //    tl->DrawLatex (xoff, 0.1, "0.1");
  //    tl->DrawLatex (xoff, 0.2, "0.2");
  //    tl->DrawLatex (xoff, 0.3, "0.3");
  //    tl->DrawLatex (xoff, 0.5, "0.5");
  //    tl->DrawLatex (xoff, 1, "1");
  //    tl->DrawLatex (xoff, 2, "2");
  //    tl->DrawLatex (xoff, 3, "3");
  //    tl->DrawLatex (xoff, 5, "5");

  //    l->SetLineStyle (2);
  //    l->SetLineWidth (2);
  //    l->DrawLine (xmin, 1, xmax, 1);
  //  }

  //  TGE* g = nullptr;

  //  g = (TGE*) tg_PHENIX_IAA_syst->Clone ();
  //  g->SetMarkerSize (0);
  //  g->SetMarkerColor (phenixColor);
  //  g->SetLineColor (phenixColor);
  //  g->SetLineWidth (1);
  //  g->SetFillColorAlpha (phenixColor, 0.3);
  //  ((TGE*) g->Clone ())->Draw ("5P");

  //  SaferDelete (&g);

  //  g = (TGE*) tg_PHENIX_IAA_stat->Clone ();
  //  g->SetMarkerStyle (phenixMarker);
  //  g->SetMarkerSize (1.8);
  //  g->SetMarkerColor (phenixColor);
  //  g->SetLineColor (phenixColor);
  //  g->SetLineWidth (3);
  //  ((TGE*) g->Clone ())->Draw ("P");

  //  if (IsFullMarker (phenixMarker)) {
  //    g->SetMarkerStyle (FullToOpenMarker (phenixMarker));
  //    g->SetMarkerSize (1.8);
  //    g->SetLineWidth (0);
  //    g->SetMarkerColor (kBlack);
  //    ((TGE*) g->Clone ())->Draw ("P");
  //  }

  //  SaferDelete (&g);

  //  g = (TGE*) tg_STAR_IAA_syst->Clone ();
  //  g->SetMarkerSize (0);
  //  g->SetLineColor (starColor);
  //  g->SetLineWidth (1);
  //  g->SetFillColorAlpha (starColor, 0.2);
  //  ((TGE*) g->Clone ())->Draw ("5P");

  //  SaferDelete (&g);

  //  g = (TGE*) tg_STAR_IAA_stat->Clone ();
  //  g->SetMarkerStyle (starMarker);
  //  g->SetMarkerSize (1.8);
  //  g->SetMarkerColor (starColor);
  //  g->SetLineColor (starColor);
  //  g->SetLineWidth (3);
  //  ((TGE*) g->Clone ())->Draw ("P");

  //  if (IsFullMarker (starMarker)) {
  //    g->SetMarkerStyle (FullToOpenMarker (starMarker));
  //    g->SetMarkerSize (1.8);
  //    g->SetLineWidth (0);
  //    g->SetMarkerColor (kBlack);
  //    ((TGE*) g->Clone ())->Draw ("P");
  //  }

  //  SaferDelete (&g);

  //  TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_iaa_syst[iPtZ][iCent]->Clone ();
  //  RecenterGraph (g_syst);
  //  ResetXErrors (g_syst);
  //  SetConstantXErrors (g_syst, 0.060, true, xmin, xmax);

  //  g_syst->SetMarkerSize (0);
  //  g_syst->SetMarkerColor (atlasColor);
  //  g_syst->SetLineColor (atlasColor);
  //  g_syst->SetLineWidth (1);
  //  g_syst->SetFillColorAlpha (atlasFillColor, 0.3);

  //  g_syst->Draw ("5P");

  //  TGAE* g_stat = make_graph (h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]);

  //  RecenterGraph (g_stat);
  //  ResetXErrors (g_stat);

  //  g_stat->SetMarkerStyle (kFullDiamond);
  //  g_stat->SetMarkerSize (2.4);
  //  g_stat->SetLineWidth (3);
  //  g_stat->SetMarkerColor (atlasColor);
  //  g_stat->SetLineColor (atlasColor);

  //  ((TGAE*) g_stat->Clone ())->Draw ("P");

  //  g_stat->SetMarkerStyle (kOpenDiamond);
  //  g_stat->SetMarkerSize (2.4);
  //  g_stat->SetLineWidth (0);
  //  g_stat->SetMarkerColor (kBlack);

  //  ((TGAE*) g_stat->Clone ())->Draw ("P");

  //  SaferDelete (&g_stat);

  //  //myMarkerAndBoxAndLineText (0.31, 0.90-0.012, 1.4, 1001, atlasFillColor, 0.30, atlasColor, kFullDiamond,  2.4, "ATLAS #it{p}_{T}^{#it{Z}} > 60 GeV, 0-10% Pb+Pb", 0.032);
  //  //myMarkerAndBoxAndLineText (0.31, 0.85-0.012, 1.4, 1001, phenixColor,    0.30, phenixColor,  phenixMarker, 1.8, "PHENIX 5 < #it{p}_{T}^{#gamma} < 9 GeV, 0-40\% Au+Au", 0.032);
  //  //myMarkerAndBoxAndLineText (0.31, 0.80-0.012, 1.4, 1001, starColor,      0.20, starColor,    starMarker, 1.8, "STAR 12 < #it{p}_{T}^{#gamma} < 20 GeV, 0-12\% Au+Au", 0.032);
  //  myMarkerAndBoxAndLineText (0.31, 0.888, 1.4, 1001, atlasFillColor, 0.30, atlasColor, kFullDiamond,  2.4, "", 0.032);
  //  myMarkerAndBoxAndLineText (0.31, 0.838, 1.4, 1001, phenixColor,    0.30, phenixColor,  phenixMarker, 1.8, "", 0.032);
  //  myMarkerAndBoxAndLineText (0.31, 0.788, 1.4, 1001, starColor,      0.20, starColor,    starMarker, 1.8, "", 0.032);

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextSize (26);
  //  tl->SetTextAlign (11);
  //  tl->DrawLatexNDC (0.316, 0.85+0.036, "ATLAS #it{p}_{T}^{#it{Z}} > 60 GeV, 0-10% Pb+Pb");
  //  tl->DrawLatexNDC (0.316, 0.80+0.036, "PHENIX 5 < #it{p}_{T}^{#gamma} < 9 GeV, 0-40\% Au+Au");
  //  tl->DrawLatexNDC (0.316, 0.75+0.036, "STAR 12 < #it{p}_{T}^{#gamma} < 20 GeV, 0-12\% Au+Au");

  //  tl->SetTextSize (32);
  //  tl->DrawLatexNDC (0.24, 0.310, "#bf{#it{ATLAS}} Internal");
  //  tl->SetTextSize (30);
  //  tl->DrawLatexNDC (0.26, 0.260, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.26, 0.210, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

  //  c5->SaveAs (Form ("%s/FinalPlots/iaa_xhz_rhicComp.pdf", plotPath.Data ()));
  //}




  //{
  //  const char* canvasName = "c6";
  //  TCanvas* c6 = new TCanvas (canvasName, "", 1600, 1200);

  //  const double llMargin = 0.17;
  //  const double lrMargin = 0.032;
  //  const double clMargin = 0.032;
  //  const double crMargin = 0.032;
  //  const double rlMargin = 0.032;
  //  const double rrMargin = 0.040;
  //  const double bMargin = 0.15;
  //  const double tMargin = 0.04;

  //  const double deltaL = (1. - llMargin - lrMargin);
  //  const double deltaC = (1. - clMargin - crMargin);
  //  const double deltaR = (1. - rlMargin - rrMargin);

  //  const double a = (double) (deltaR * deltaC / (deltaL*deltaR + deltaC*deltaR + deltaL*deltaC));
  //  const double b = (double) (deltaR * deltaL / (deltaL*deltaR + deltaC*deltaR + deltaL*deltaC));

  //  const double xPadLCMiddle = a;
  //  const double xPadCRMiddle = a+b;

  //  const double yPadMiddle = 0.5;

  //  TPad* luPad = nullptr;
  //  TPad* cuPad = nullptr;
  //  TPad* ruPad = nullptr;
  //  TPad* ldPad = nullptr;
  //  TPad* cdPad = nullptr;
  //  TPad* rdPad = nullptr;

  //  luPad = new TPad (Form ("%s_luPad", canvasName), "", 0, yPadMiddle, xPadLCMiddle, 1);
  //  cuPad = new TPad (Form ("%s_cuPad", canvasName), "", xPadLCMiddle, yPadMiddle, xPadCRMiddle, 1);
  //  ruPad = new TPad (Form ("%s_ruPad", canvasName), "", xPadCRMiddle, yPadMiddle, 1, 1);
  //  ldPad = new TPad (Form ("%s_ldPad", canvasName), "", 0, 0, xPadLCMiddle, yPadMiddle);
  //  cdPad = new TPad (Form ("%s_cdPad", canvasName), "", xPadLCMiddle, 0, xPadCRMiddle, yPadMiddle);
  //  rdPad = new TPad (Form ("%s_rdPad", canvasName), "", xPadCRMiddle, 0, 1, yPadMiddle);

  //  luPad->SetLeftMargin (llMargin);
  //  luPad->SetRightMargin (lrMargin);
  //  cuPad->SetLeftMargin (clMargin);
  //  cuPad->SetRightMargin (crMargin);
  //  ruPad->SetLeftMargin (rlMargin);
  //  ruPad->SetRightMargin (rrMargin);
  //  luPad->SetBottomMargin (bMargin);
  //  luPad->SetTopMargin (tMargin);
  //  cuPad->SetBottomMargin (bMargin);
  //  cuPad->SetTopMargin (tMargin);
  //  ruPad->SetBottomMargin (bMargin);
  //  ruPad->SetTopMargin (tMargin);
  //  ldPad->SetLeftMargin (llMargin);
  //  ldPad->SetRightMargin (lrMargin);
  //  cdPad->SetLeftMargin (clMargin);
  //  cdPad->SetRightMargin (crMargin);
  //  rdPad->SetLeftMargin (rlMargin);
  //  rdPad->SetRightMargin (rrMargin);
  //  ldPad->SetBottomMargin (bMargin);
  //  ldPad->SetTopMargin (tMargin);
  //  cdPad->SetBottomMargin (bMargin);
  //  cdPad->SetTopMargin (tMargin);
  //  rdPad->SetBottomMargin (bMargin);
  //  rdPad->SetTopMargin (tMargin);

  //  TPad* uPads[3] = {luPad, cuPad, ruPad};
  //  TPad* dPads[3] = {ldPad, cdPad, rdPad};

  //  for (int i = 0; i < 3; i++) {
  //    uPads[i]->Draw ();
  //    dPads[i]->Draw ();
  //  }
  //  for (int i = 0; i < 3; i++) {
  //    uPads[i]->SetLogx ();
  //    uPads[i]->SetLogy ();
  //    dPads[i]->SetLogx ();
  //    dPads[i]->SetLogy ();
  //  }


  //  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
  //    uPads[iPtZ-2]->cd ();

  //    const double xmin = pTchBins[iPtZ][0];
  //    const double xmax = pTchBins[iPtZ][nPtchBins[iPtZ]];

  //    TH1D* h = new TH1D ("", "", 1, xmin, xmax);

  //    TAxis* xax = h->GetXaxis ();
  //    TAxis* yax = h->GetYaxis ();

  //    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
  //    xax->SetRangeUser (xmin, xmax);
  //    xax->SetTitleFont (43);
  //    xax->SetTitleSize (30);
  //    xax->SetTitleOffset (2.5);
  //    xax->SetLabelSize (0);

  //    yax->SetTitle ("(1/N_{Z}) (d^{2}N_{ch} / d#it{p}_{T} d#Delta#phi) [GeV^{-1}]");
  //    yax->SetRangeUser (2e-3, 200);
  //    yax->SetTitleFont (43);
  //    yax->SetTitleSize (30);
  //    yax->SetTitleOffset (3.2);
  //    yax->SetLabelFont (43);
  //    if (iPtZ == 2) yax->SetLabelSize (28);
  //    else yax->SetLabelSize (0);

  //    h->SetLineWidth (0);

  //    h->DrawCopy ("");
  //    SaferDelete (&h);

  //    tl->SetTextFont (43);
  //    tl->SetTextSize (28);
  //    tl->SetTextAlign (21);
  //    const double yoff = 0.0009;
  //    tl->DrawLatex (1,  yoff, "1");
  //    tl->DrawLatex (2,  yoff, "2");
  //    tl->DrawLatex (3,  yoff, "3");
  //    tl->DrawLatex (4,  yoff, "4");
  //    tl->DrawLatex (5,  yoff, "5");
  //    tl->DrawLatex (6,  yoff, "6");
  //    tl->DrawLatex (7,  yoff, "7");
  //    tl->DrawLatex (10, yoff, "10");
  //    tl->DrawLatex (20, yoff, "20");
  //    tl->DrawLatex (30, yoff, "30");
  //    tl->DrawLatex (40, yoff, "40");
  //    tl->DrawLatex (60, yoff, "60");
  //  } // end loop over iPtZ

  //  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
  //    uPads[iPtZ-2]->cd ();

  //    const double xmin = pTchBins[iPtZ][0];
  //    const double xmax = pTchBins[iPtZ][nPtchBins[iPtZ]];

  //    for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
  //      TGAE* g_syst = (TGAE*) g_trk_pt_ptz_sub_syst[iPtZ][iCent]->Clone ();

  //      RecenterGraph (g_syst);
  //      ResetXErrors (g_syst);
  //      const short iCentDelta = (iCent == 0 ? numCentBins-1 : iCent-1);
  //      deltaize (g_syst, 0.05*(iCentDelta - 0.5*(numCentBins-1)), true, xmin, xmax);
  //      ResetXErrors (g_syst);
  //      SetConstantXErrors (g_syst, 0.04, true, xmin, xmax);

  //      g_syst->SetMarkerSize (0);
  //      g_syst->SetLineWidth (1);
  //      g_syst->SetMarkerColor (finalColors[iCent]);
  //      g_syst->SetLineColor (finalColors[iCent]);
  //      g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

  //      ((TGAE*) g_syst->Clone ())->Draw ("5P");

  //      SaferDelete (&g_syst);
  //    } // end loop over iCent

  //    for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
  //      TGAE* g_stat = make_graph (h_trk_pt_ptz_sub_stat[iPtZ][iCent]);

  //      Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
  //      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

  //      RecenterGraph (g_stat);
  //      ResetXErrors (g_stat);
  //      const short iCentDelta = (iCent == 0 ? numCentBins-1 : iCent-1);
  //      deltaize (g_stat, 0.05*(iCentDelta - 0.5*(numCentBins-1)), true, xmin, xmax);
  //      ResetXErrors (g_stat);

  //      g_stat->SetMarkerStyle (markerStyle);
  //      g_stat->SetMarkerSize (markerSize);
  //      g_stat->SetLineWidth (3);
  //      g_stat->SetMarkerColor (finalColors[iCent]);
  //      g_stat->SetLineColor (finalColors[iCent]);

  //      ((TGAE*) g_stat->Clone ())->Draw ("P");

  //      if (iCent != 0) {
  //        markerStyle = FullToOpenMarker (markerStyle);
  //        markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
  //        
  //        g_stat->SetMarkerStyle (markerStyle);
  //        g_stat->SetMarkerSize (markerSize);
  //        if (iCent > 0) g_stat->SetLineWidth (0);
  //        g_stat->SetMarkerColor (kBlack);

  //        ((TGAE*) g_stat->Clone ())->Draw ("P");
  //      }
  //      else {
  //        TLine* g_stat_line = new TLine ();
  //        g_stat_line->SetLineWidth (2);
  //        g_stat_line->SetLineColor (kBlack);
  //        for (int i = 0; i < g_stat->GetN (); i++) {
  //          double x, y, yhi, ylo;
  //          g_stat->GetPoint (i, x, y);
  //          yhi = y + g_stat->GetErrorYhigh (i);
  //          ylo = y - g_stat->GetErrorYlow (i);
  //          g_stat_line->DrawLine (x, ylo, x, yhi);
  //        }
  //        SaferDelete (&g_stat_line);
  //      }

  //      SaferDelete (&g_stat);
  //    } // end loop over iCent
  //  } // end loop over iPtZ


  //  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
  //    dPads[iPtZ-2]->cd ();

  //    const double xmin = (iPtZ == 2 ? 1./8. : xhZBins[iPtZ][0]);
  //    const double xmax = xhZBins[iPtZ][nXhZBins[iPtZ]];

  //    TH1D* h = new TH1D ("", "", 1, xmin, xmax);

  //    TAxis* xax = h->GetXaxis ();
  //    TAxis* yax = h->GetYaxis ();

  //    xax->SetTitle ("#it{x}_{h#it{Z}} = #it{p}_{T}^{ch} #/#it{p}_{T}^{#it{Z}}");
  //    xax->SetRangeUser (xmin, xmax);
  //    xax->SetTitleFont (43);
  //    xax->SetTitleSize (30);
  //    xax->SetTitleOffset (2.5);
  //    xax->SetLabelSize (0);

  //    yax->SetTitle ("#it{I}_{AA} (#it{x}_{h#it{Z}})");
  //    yax->SetTitle ("(1/N_{Z}) (d^{2}N_{ch} / d#it{x} d#Delta#phi)");
  //    yax->SetRangeUser (5e-3, 4e3);
  //    yax->SetTitleFont (43);
  //    yax->SetTitleSize (30);
  //    yax->SetTitleOffset (3.2);
  //    yax->SetLabelFont (43);
  //    if (iPtZ == 2) yax->SetLabelSize (28);
  //    else yax->SetLabelSize (0);

  //    h->SetLineWidth (0);

  //    h->DrawCopy ("");
  //    SaferDelete (&h);

  //    tl->SetTextFont (43);
  //    tl->SetTextSize (28);
  //    tl->SetTextAlign (21);
  //    const double yoff = 1.95e-3;
  //    if (iPtZ > 2) {
  //      if (iPtZ > 3) tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
  //      else tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
  //    }
  //    if (iPtZ >= 3) tl->DrawLatex (1e-1,  yoff, "10^{-1}");
  //    tl->DrawLatex (2e-1,  yoff, "2#times10^{-1}");
  //    tl->DrawLatex (5e-1,  yoff, "5#times10^{-1}");
  //    tl->DrawLatex (1,     yoff, "1");
  //  } // end loop over iPtZ

  //  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
  //    dPads[iPtZ-2]->cd ();

  //    const double xmin = (iPtZ == 2 ? 1./8. : xhZBins[iPtZ][0]);
  //    const double xmax = xhZBins[iPtZ][nXhZBins[iPtZ]];

  //    for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
  //      TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_sub_syst[iPtZ][iCent]->Clone ();

  //      RecenterGraph (g_syst);
  //      ResetXErrors (g_syst);
  //      const short iCentDelta = (iCent == 0 ? numCentBins-1 : iCent-1);
  //      deltaize (g_syst, 0.05*(iCentDelta - 0.5*(numCentBins-1)), true, xmin, xmax);
  //      ResetXErrors (g_syst);
  //      SetConstantXErrors (g_syst, 0.04, true, xmin, xmax);

  //      g_syst->SetMarkerSize (0);
  //      g_syst->SetLineWidth (1);
  //      g_syst->SetMarkerColor (finalColors[iCent]);
  //      g_syst->SetLineColor (finalColors[iCent]);
  //      g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

  //      ((TGAE*) g_syst->Clone ())->Draw ("5P");

  //      SaferDelete (&g_syst);
  //    } // end loop over iCent

  //    for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
  //      TGAE* g_stat = make_graph (h_trk_xhz_ptz_sub_stat[iPtZ][iCent]);

  //      Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
  //      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

  //      RecenterGraph (g_stat);
  //      ResetXErrors (g_stat);
  //      const short iCentDelta = (iCent == 0 ? numCentBins-1 : iCent-1);
  //      deltaize (g_stat, 0.05*(iCentDelta - 0.5*(numCentBins-1)), true, xmin, xmax);
  //      ResetXErrors (g_stat);

  //      g_stat->SetMarkerStyle (markerStyle);
  //      g_stat->SetMarkerSize (markerSize);
  //      g_stat->SetLineWidth (3);
  //      g_stat->SetMarkerColor (finalColors[iCent]);
  //      g_stat->SetLineColor (finalColors[iCent]);

  //      ((TGAE*) g_stat->Clone ())->Draw ("P");

  //      if (iCent != 0) {
  //        markerStyle = FullToOpenMarker (markerStyle);
  //        markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
  //        
  //        g_stat->SetMarkerStyle (markerStyle);
  //        g_stat->SetMarkerSize (markerSize);
  //        if (iCent > 0) g_stat->SetLineWidth (0);
  //        g_stat->SetMarkerColor (kBlack);

  //        ((TGAE*) g_stat->Clone ())->Draw ("P");
  //      }
  //      else {
  //        TLine* g_stat_line = new TLine ();
  //        g_stat_line->SetLineWidth (2);
  //        g_stat_line->SetLineColor (kBlack);
  //        for (int i = 0; i < g_stat->GetN (); i++) {
  //          double x, y, yhi, ylo;
  //          g_stat->GetPoint (i, x, y);
  //          yhi = y + g_stat->GetErrorYhigh (i);
  //          ylo = y - g_stat->GetErrorYlow (i);
  //          g_stat_line->DrawLine (x, ylo, x, yhi);
  //        }
  //        SaferDelete (&g_stat_line);
  //      }

  //      SaferDelete (&g_stat);
  //    } // end loop over iCent
  //  } // end loop over iPtZ

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextAlign (11);

  //  luPad->cd ();
  //  tl->SetTextSize (32);
  //  tl->DrawLatexNDC (0.21, 0.86, "#bf{#it{ATLAS}} Internal");
  //  tl->SetTextSize (28);
  //  tl->DrawLatexNDC (0.21, 0.79, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.21, 0.73, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

  //  ldPad->cd ();
  //  tl->SetTextSize (32);
  //  tl->DrawLatexNDC (0.21, 0.86, "#bf{#it{ATLAS}} Internal");
  //  tl->SetTextSize (28);
  //  tl->DrawLatexNDC (0.21, 0.79, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.21, 0.73, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

  //  tl->SetTextSize (30);
  //  luPad->cd ();
  //  tl->DrawLatexNDC (0.56, 0.86, "15 < #it{p}_{T}^{#it{Z}} < 30 GeV");
  //  cuPad->cd ();
  //  tl->DrawLatexNDC (0.50, 0.86, "30 < #it{p}_{T}^{#it{Z}} < 60 GeV");
  //  ruPad->cd ();
  //  tl->DrawLatexNDC (0.60, 0.86, "#it{p}_{T}^{#it{Z}} > 60 GeV");

  //  ldPad->cd ();
  //  tl->DrawLatexNDC (0.56, 0.86, "15 < #it{p}_{T}^{#it{Z}} < 30 GeV");
  //  cdPad->cd ();
  //  tl->DrawLatexNDC (0.50, 0.86, "30 < #it{p}_{T}^{#it{Z}} < 60 GeV");
  //  rdPad->cd ();
  //  tl->DrawLatexNDC (0.60, 0.86, "#it{p}_{T}^{#it{Z}} > 60 GeV");

  //  cuPad->cd ();
  //  myMarkerAndBoxAndLineText (0.22, 0.280, 3.0, 1001, finalFillColors[0], 0.30, finalColors[0], kOpenCircle,     1.8, "#it{pp}", 0.016 / (gPad->GetWNDC ()));
  //  myMarkerAndBoxAndLineText (0.22, 0.220, 3.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "Pb+Pb 30-80\%", 0.016 / (gPad->GetWNDC ()));
  //  ruPad->cd ();
  //  myMarkerAndBoxAndLineText (0.22, 0.280, 3.0, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.8, "Pb+Pb 10-30\%", 0.016 / (gPad->GetWNDC ()));
  //  myMarkerAndBoxAndLineText (0.22, 0.220, 3.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "Pb+Pb 0-10\%", 0.016 / (gPad->GetWNDC ()));

  //  cdPad->cd ();
  //  myMarkerAndBoxAndLineText (0.22, 0.280, 3.0, 1001, finalFillColors[0], 0.30, finalColors[0], kOpenCircle,     1.8, "#it{pp}", 0.016 / (gPad->GetWNDC ()));
  //  myMarkerAndBoxAndLineText (0.22, 0.220, 3.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "Pb+Pb 30-80\%", 0.016 / (gPad->GetWNDC ()));
  //  rdPad->cd ();
  //  myMarkerAndBoxAndLineText (0.22, 0.280, 3.0, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.8, "Pb+Pb 10-30\%", 0.016 / (gPad->GetWNDC ()));
  //  myMarkerAndBoxAndLineText (0.22, 0.220, 3.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "Pb+Pb 0-10\%", 0.016 / (gPad->GetWNDC ()));

  //  c6->SaveAs (Form ("%s/FinalPlots/yield_allptz.pdf", plotPath.Data ()));
  //}



  //{
  //  TCanvas* c7 = new TCanvas ("c7", "", 800, 960);

  //  const double lMargin = 0.12;
  //  const double rMargin = 0.04;
  //  const double bMargin = 0.10;
  //  const double tMargin = 0.04;

  //  c7->SetLeftMargin (lMargin);
  //  c7->SetRightMargin (rMargin);
  //  c7->SetBottomMargin (bMargin);
  //  c7->SetTopMargin (tMargin);

  //  c7->SetLogx ();
  //  c7->SetLogy ();

  //  const double xmin = pTchBins[nPtZBins-1][0];
  //  const double xmax = pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]];

  //  TH1D* h = new TH1D ("", "", 1, xmin, xmax);

  //  TAxis* xax = h->GetXaxis ();
  //  TAxis* yax = h->GetYaxis ();

  //  xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
  //  xax->SetRangeUser (xmin, xmax);
  //  xax->SetTitleFont (43);
  //  xax->SetTitleSize (30);
  //  xax->SetTitleOffset (1.3);
  //  xax->SetLabelSize (0);

  //  yax->SetTitle ("(1/N_{Z}) (d^{2}N_{ch} / d#it{p}_{T} d#Delta#phi) [GeV^{-1}]");
  //  double ymin = 2e-3;
  //  double ymax = 2e3;
  //  yax->SetRangeUser (ymin, ymax);
  //  yax->SetTitleFont (43);
  //  yax->SetTitleSize (30);
  //  yax->SetTitleOffset (1.8);
  //  yax->SetLabelFont (43);
  //  yax->SetLabelSize (32);

  //  h->SetLineWidth (0);

  //  h->DrawCopy ("");
  //  SaferDelete (&h);

  //  tl->SetTextFont (43);
  //  tl->SetTextSize (32);
  //  tl->SetTextAlign (21);
  //  const double yoff = ymin / exp (0.034 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
  //  tl->DrawLatex (1,  yoff, "1");
  //  tl->DrawLatex (2,  yoff, "2");
  //  tl->DrawLatex (3,  yoff, "3");
  //  tl->DrawLatex (4,  yoff, "4");
  //  tl->DrawLatex (5,  yoff, "5");
  //  tl->DrawLatex (6,  yoff, "6");
  //  tl->DrawLatex (7,  yoff, "7");
  //  tl->DrawLatex (10, yoff, "10");
  //  tl->DrawLatex (20, yoff, "20");
  //  tl->DrawLatex (30, yoff, "30");
  //  tl->DrawLatex (40, yoff, "40");
  //  tl->DrawLatex (60, yoff, "60");

  //  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
  //    for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
  //      TGAE* g_syst = (TGAE*) g_trk_pt_ptz_sub_syst[iPtZ][iCent]->Clone ();

  //      OffsetYAxis (g_syst, pow (10, iPtZ-2), true);
  //      RecenterGraph (g_syst);
  //      ResetXErrors (g_syst);
  //      const short iCentDelta = (iCent == 0 ? numCentBins-1 : iCent-1);
  //      deltaize (g_syst, 0.05*(iCentDelta - 0.5*(numCentBins-1)), true, xmin, xmax);
  //      ResetXErrors (g_syst);
  //      SetConstantXErrors (g_syst, 0.04, true, xmin, xmax);

  //      g_syst->SetMarkerSize (0);
  //      g_syst->SetLineWidth (1);
  //      g_syst->SetMarkerColor (finalColors[iCent]);
  //      g_syst->SetLineColor (finalColors[iCent]);
  //      g_syst->SetFillColorAlpha (finalFillColors[iCent], 0.3);

  //      ((TGAE*) g_syst->Clone ())->Draw ("5P");

  //      SaferDelete (&g_syst);
  //    } // end loop over iCent

  //    for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
  //      TGAE* g_stat = make_graph (h_trk_pt_ptz_sub_stat[iPtZ][iCent]);

  //      Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
  //      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);

  //      OffsetYAxis (g_stat, pow (10, iPtZ-2), true);
  //      RecenterGraph (g_stat);
  //      ResetXErrors (g_stat);
  //      const short iCentDelta = (iCent == 0 ? numCentBins-1 : iCent-1);
  //      deltaize (g_stat, 0.05*(iCentDelta - 0.5*(numCentBins-1)), true, xmin, xmax);
  //      ResetXErrors (g_stat);

  //      g_stat->SetMarkerStyle (markerStyle);
  //      g_stat->SetMarkerSize (markerSize);
  //      g_stat->SetLineWidth (3);
  //      g_stat->SetMarkerColor (finalColors[iCent]);
  //      g_stat->SetLineColor (finalColors[iCent]);

  //      ((TGAE*) g_stat->Clone ())->Draw ("P");

  //      if (iCent != 0) {
  //        markerStyle = FullToOpenMarker (markerStyle);
  //        markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);
  //        
  //        g_stat->SetMarkerStyle (markerStyle);
  //        g_stat->SetMarkerSize (markerSize);
  //        if (iCent > 0) g_stat->SetLineWidth (0);
  //        g_stat->SetMarkerColor (kBlack);

  //        ((TGAE*) g_stat->Clone ())->Draw ("P");
  //      }
  //      else {
  //        TLine* g_stat_line = new TLine ();
  //        g_stat_line->SetLineWidth (2);
  //        g_stat_line->SetLineColor (kBlack);
  //        for (int i = 0; i < g_stat->GetN (); i++) {
  //          double x, y, yhi, ylo;
  //          g_stat->GetPoint (i, x, y);
  //          yhi = y + g_stat->GetErrorYhigh (i);
  //          ylo = y - g_stat->GetErrorYlow (i);
  //          g_stat_line->DrawLine (x, ylo, x, yhi);
  //        }
  //        SaferDelete (&g_stat_line);
  //      }

  //      SaferDelete (&g_stat);
  //    } // end loop over iCent
  //  } // end loop over iPtZ

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextAlign (11);

  //  tl->SetTextSize (30);
  //  tl->DrawLatexNDC (0.37, 0.900, "#bf{#it{ATLAS}} Internal");
  //  tl->SetTextSize (28);
  //  tl->DrawLatexNDC (0.39, 0.8625, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.39, 0.825, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

  //  tl->SetTextSize (24);
  //  tl->SetTextAngle (322);
  //  tl->DrawLatex (2.5, 2.1, "15 < #it{p}_{T}^{#it{Z}} < 30 GeV (#times 1)");
  //  tl->DrawLatex (5.5, 8.6, "30 < #it{p}_{T}^{#it{Z}} < 60 GeV (#times 10)");
  //  tl->DrawLatex (14, 30, "#it{p}_{T}^{#it{Z}} > 60 GeV (#times 10^{2})");
  //  tl->SetTextAngle (0);

  //  tl->SetTextSize (28);
  //  tl->SetTextAlign (11);
  //  tl->DrawLatexNDC (0.345, 0.270-0.009, "#it{pp}");
  //  tl->DrawLatexNDC (0.345, 0.230-0.009, "30-80\%");
  //  tl->DrawLatexNDC (0.345, 0.190-0.009, "10-30\%");
  //  tl->DrawLatexNDC (0.345, 0.150-0.009, "0-10\%");
  //  MakeDataBox   (0.36, 0.270, finalFillColors[0], 0.30, kOpenCircle, 1.6);
  //  MakeDataBox   (0.36, 0.230, finalFillColors[1], 0.30, markerStyles[0], 1.6);
  //  MakeDataBox   (0.36, 0.190, finalFillColors[2], 0.30, markerStyles[1], 1.6);
  //  MakeDataBox   (0.36, 0.150, finalFillColors[3], 0.30, markerStyles[2], 2.2);

  //  c7->SaveAs (Form ("%s/FinalPlots/yield_allptz_pTchOnly_onePlot.pdf", plotPath.Data ()));
  //}




  //{
  //  TCanvas* c8 = new TCanvas ("c8", "", 800, 800);

  //  const double lMargin = 0.15;
  //  const double rMargin = 0.04;
  //  const double bMargin = 0.15;
  //  const double tMargin = 0.04;

  //  c8->SetLeftMargin (lMargin);
  //  c8->SetRightMargin (rMargin);
  //  c8->SetBottomMargin (bMargin);
  //  c8->SetTopMargin (tMargin);

  //  c8->SetLogx ();

  //  short iPtZ = 4;
  //  const short iCent = 3;

  //  const double xmin = pTchBins[iPtZ][0];
  //  const double xmax = pTchBins[iPtZ][nPtchBins[iPtZ]];
  //  const double ymin = 0.;
  //  const double ymax = 3.6;

  //  {
  //    TH1D* h = new TH1D ("", "", 1, xmin, xmax);

  //    TAxis* xax = h->GetXaxis ();
  //    TAxis* yax = h->GetYaxis ();

  //    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
  //    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
  //    xax->SetRangeUser (xmin, xmax);
  //    xax->SetLabelSize (0);

  //    yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
  //    yax->SetTitleOffset (1.0 * yax->GetTitleOffset ());
  //    yax->SetLabelSize (0);
  //    yax->SetRangeUser (ymin, ymax);

  //    h->SetLineWidth (0);

  //    h->DrawCopy ("");
  //    SaferDelete (&h);

  //    tl->SetTextFont (43);
  //    tl->SetTextSize (36);
  //    tl->SetTextAlign (21);

  //    const double yoff = ymin - 0.05 * (ymax-ymin) / (1.-tMargin-bMargin);
  //    tl->DrawLatex (1,  yoff, "1");
  //    tl->DrawLatex (2,  yoff, "2");
  //    tl->DrawLatex (3,  yoff, "3");
  //    tl->DrawLatex (4,  yoff, "4");
  //    tl->DrawLatex (5,  yoff, "5");
  //    tl->DrawLatex (6,  yoff, "6");
  //    tl->DrawLatex (7,  yoff, "7");
  //    tl->DrawLatex (10, yoff, "10");
  //    tl->DrawLatex (20, yoff, "20");
  //    tl->DrawLatex (30, yoff, "30");
  //    tl->DrawLatex (40, yoff, "40");
  //    tl->DrawLatex (60, yoff, "60");

  //    tl->SetTextAlign (32);

  //    const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-lMargin-rMargin));
  //    tl->DrawLatex (xoff, 0, "0");
  //    tl->DrawLatex (xoff, 0.5, "0.5");
  //    tl->DrawLatex (xoff, 1, "1");
  //    tl->DrawLatex (xoff, 1.5, "1.5");
  //    tl->DrawLatex (xoff, 2, "2");
  //    tl->DrawLatex (xoff, 2.5, "2.5");
  //    tl->DrawLatex (xoff, 3, "3");
  //    tl->DrawLatex (xoff, 3.5, "3.5");

  //    l->SetLineStyle (2);
  //    l->SetLineWidth (2);
  //    l->SetLineColor (kBlack);
  //    l->DrawLine (xmin, 1, xmax, 1);
  //  }

  //  for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
  //    TGAE* g = (TGAE*) g_hybrid_pth[iPtZ]->Clone ();
  //    SetMinErrors (g, minModelUnc, false);
  //    ResetXErrors (g);
  //    deltaize (g, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, xmin, xmax);
  //    ResetXErrors (g);

  //    g->SetFillColorAlpha (finalModelFillColors[iPtZ-1], 0.4);
  //    ((TGAE*) g->Clone ())->Draw ("3");
  //    ResetTGAEErrors (g);
  //    ResetXErrors (g);
  //    g->SetLineColor (finalModelFillColors[iPtZ-1]);
  //    g->SetLineWidth (4);
  //    ((TGAE*) g->Clone ())->Draw ("L");
  //    SaferDelete (&g);
  //  }

  //  for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
  //    TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

  //    RecenterGraph (g_syst);
  //    ResetXErrors (g_syst);
  //    deltaize (g_syst, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, xmin, xmax);
  //    ResetXErrors (g_syst);
  //    SetConstantXErrors (g_syst, 0.04, true, xmin, xmax);

  //    g_syst->SetMarkerSize (0);
  //    g_syst->SetLineWidth (1);
  //    g_syst->SetLineColor (finalColors[iPtZ-1]);
  //    g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

  //    g_syst->Draw ("5P");
  //  }

  //  for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
  //    TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

  //    RecenterGraph (g_stat);
  //    ResetXErrors (g_stat);
  //    deltaize (g_stat, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, xmin, xmax);
  //    ResetXErrors (g_stat);

  //    Style_t markerStyle = markerStyles[iPtZ-2];
  //    float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

  //    g_stat->SetMarkerStyle (markerStyle);
  //    g_stat->SetMarkerSize (markerSize);
  //    g_stat->SetLineWidth (3);
  //    g_stat->SetMarkerColor (finalColors[iPtZ-1]);
  //    g_stat->SetLineColor (finalColors[iPtZ-1]);

  //    ((TGAE*) g_stat->Clone ())->Draw ("P");

  //    markerStyle = FullToOpenMarker (markerStyle);
  //    markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

  //    g_stat->SetMarkerStyle (markerStyle);
  //    g_stat->SetMarkerSize (markerSize);
  //    g_stat->SetLineWidth (0);
  //    g_stat->SetMarkerColor (kBlack);

  //    ((TGAE*) g_stat->Clone ())->Draw ("P");

  //    SaferDelete (&g_stat);
  //  }

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextAlign (11);

  //  tl->SetTextSize (32);
  //  tl->DrawLatexNDC (0.31, 0.890, "#bf{#it{ATLAS}} Internal");
  //  tl->SetTextSize (30);
  //  tl->DrawLatexNDC (0.34, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.34, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

  //  tl->SetTextSize (27);
  //  tl->SetTextAlign (21);
  //  tl->DrawLatexNDC (0.418, 0.74, "15-30");
  //  tl->DrawLatexNDC (0.518, 0.74, "30-60");
  //  tl->DrawLatexNDC (0.610, 0.74, "60+");

  //  tl->SetTextAlign (11);
  //  tl->DrawLatexNDC (0.66, 0.74, "#it{p}_{T}^{#it{Z}} [GeV]");
  //  tl->DrawLatexNDC (0.66, 0.685+0.007, "ATLAS 0-10\%");
  //  tl->DrawLatexNDC (0.66, 0.630+0.014, "Hybrid Model");
  //  myMarkerAndBoxAndLineText (0.654, 0.685+0.007, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.560, 0.685+0.007, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.8, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.464, 0.685+0.007, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.654, 0.630+0.014, 1.4, 1001, finalModelFillColors[3], 0.40, finalModelFillColors[3], -1, 1, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.560, 0.630+0.014, 1.4, 1001, finalModelFillColors[2], 0.40, finalModelFillColors[2], -1, 1, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.464, 0.630+0.014, 1.4, 1001, finalModelFillColors[1], 0.40, finalModelFillColors[1], -1, 1, "", 0.045);

  //  l->SetLineStyle (1);
  //  l->SetLineWidth (4);
  //  l->SetLineColor (finalModelFillColors[3]);
  //  l->DrawLineNDC (0.654 - (0.8*0.045) + 0.02 - (0.04*1.4), 0.630+0.014+(0.25*0.045), 0.654 - (0.8*0.045) + 0.02, 0.630+0.014+(0.25*0.045));
  //  l->SetLineColor (finalModelFillColors[2]);
  //  l->DrawLineNDC (0.560 - (0.8*0.045) + 0.02 - (0.04*1.4), 0.630+0.014+(0.25*0.045), 0.560 - (0.8*0.045) + 0.02, 0.630+0.014+(0.25*0.045));
  //  l->SetLineColor (finalModelFillColors[1]);
  //  l->DrawLineNDC (0.464 - (0.8*0.045) + 0.02 - (0.04*1.4), 0.630+0.014+(0.25*0.045), 0.464 - (0.8*0.045) + 0.02, 0.630+0.014+(0.25*0.045));

  //  c8->SaveAs (Form ("%s/FinalPlots/iaa_ptch_hybridComp.pdf", plotPath.Data ()));
  //}





  //{
  //  TCanvas* c9 = new TCanvas ("c9", "", 800, 800);

  //  const double lMargin = 0.15;
  //  const double rMargin = 0.04;
  //  const double bMargin = 0.15;
  //  const double tMargin = 0.04;

  //  c9->SetLeftMargin (lMargin);
  //  c9->SetRightMargin (rMargin);
  //  c9->SetBottomMargin (bMargin);
  //  c9->SetTopMargin (tMargin);

  //  c9->SetLogx ();

  //  short iPtZ = 4;
  //  const short iCent = 3;

  //  const double xmin = pTchBins[iPtZ][0];
  //  const double xmax = pTchBins[iPtZ][nPtchBins[iPtZ]];
  //  const double ymin = 0.;
  //  const double ymax = 3.6;

  //  {
  //    TH1D* h = new TH1D ("", "", 1, xmin, xmax);

  //    TAxis* xax = h->GetXaxis ();
  //    TAxis* yax = h->GetYaxis ();

  //    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
  //    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
  //    xax->SetRangeUser (xmin, xmax);
  //    xax->SetLabelSize (0);

  //    yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
  //    yax->SetTitleOffset (1.0 * yax->GetTitleOffset ());
  //    yax->SetLabelSize (0);
  //    yax->SetRangeUser (ymin, ymax);

  //    h->SetLineWidth (0);

  //    h->DrawCopy ("");
  //    SaferDelete (&h);

  //    tl->SetTextFont (43);
  //    tl->SetTextSize (36);
  //    tl->SetTextAlign (21);

  //    const double yoff = ymin - 0.05 * (ymax-ymin) / (1.-tMargin-bMargin);
  //    tl->DrawLatex (1,  yoff, "1");
  //    tl->DrawLatex (2,  yoff, "2");
  //    tl->DrawLatex (3,  yoff, "3");
  //    tl->DrawLatex (4,  yoff, "4");
  //    tl->DrawLatex (5,  yoff, "5");
  //    tl->DrawLatex (6,  yoff, "6");
  //    tl->DrawLatex (7,  yoff, "7");
  //    tl->DrawLatex (10, yoff, "10");
  //    tl->DrawLatex (20, yoff, "20");
  //    tl->DrawLatex (30, yoff, "30");
  //    tl->DrawLatex (40, yoff, "40");
  //    tl->DrawLatex (60, yoff, "60");

  //    tl->SetTextAlign (32);

  //    const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-lMargin-rMargin));
  //    tl->DrawLatex (xoff, 0, "0");
  //    tl->DrawLatex (xoff, 0.5, "0.5");
  //    tl->DrawLatex (xoff, 1, "1");
  //    tl->DrawLatex (xoff, 1.5, "1.5");
  //    tl->DrawLatex (xoff, 2, "2");
  //    tl->DrawLatex (xoff, 2.5, "2.5");
  //    tl->DrawLatex (xoff, 3, "3");
  //    tl->DrawLatex (xoff, 3.5, "3.5");

  //    l->SetLineStyle (2);
  //    l->SetLineWidth (2);
  //    l->SetLineColor (kBlack);
  //    l->DrawLine (xmin, 1, xmax, 1);
  //  }

  //  for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
  //    TGAE* g = (TGAE*) g_jewel_pth[iPtZ]->Clone ();
  //    SetMinErrors (g, minModelUnc, false);
  //    ResetXErrors (g);
  //    deltaize (g, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, xmin, xmax);
  //    ResetXErrors (g);

  //    g->SetFillColorAlpha (finalModelFillColors[iPtZ-1], 0.4);
  //    ((TGAE*) g->Clone ())->Draw ("3");
  //    ResetTGAEErrors (g);
  //    ResetXErrors (g);
  //    g->SetLineColor (finalModelFillColors[iPtZ-1]);
  //    g->SetLineWidth (4);
  //    ((TGAE*) g->Clone ())->Draw ("L");
  //    SaferDelete (&g);
  //  }

  //  for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
  //    TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

  //    RecenterGraph (g_syst);
  //    ResetXErrors (g_syst);
  //    deltaize (g_syst, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, xmin, xmax);
  //    ResetXErrors (g_syst);
  //    SetConstantXErrors (g_syst, 0.04, true, pTchBins[4][0], pTchBins[4][nPtchBins[4]]);

  //    g_syst->SetMarkerSize (0);
  //    g_syst->SetLineWidth (1);
  //    g_syst->SetLineColor (finalColors[iPtZ-1]);
  //    g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

  //    g_syst->Draw ("5P");
  //  }

  //  for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
  //    TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

  //    RecenterGraph (g_stat);
  //    ResetXErrors (g_stat);
  //    deltaize (g_stat, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, xmin, xmax);
  //    ResetXErrors (g_stat);

  //    Style_t markerStyle = markerStyles[iPtZ-2];
  //    float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

  //    g_stat->SetMarkerStyle (markerStyle);
  //    g_stat->SetMarkerSize (markerSize);
  //    g_stat->SetLineWidth (3);
  //    g_stat->SetMarkerColor (finalColors[iPtZ-1]);
  //    g_stat->SetLineColor (finalColors[iPtZ-1]);

  //    ((TGAE*) g_stat->Clone ())->Draw ("P");

  //    markerStyle = FullToOpenMarker (markerStyle);
  //    markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

  //    g_stat->SetMarkerStyle (markerStyle);
  //    g_stat->SetMarkerSize (markerSize);
  //    g_stat->SetLineWidth (0);
  //    g_stat->SetMarkerColor (kBlack);

  //    ((TGAE*) g_stat->Clone ())->Draw ("P");

  //    SaferDelete (&g_stat);
  //  }

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextAlign (11);

  //  tl->SetTextSize (32);
  //  tl->DrawLatexNDC (0.31, 0.890, "#bf{#it{ATLAS}} Internal");
  //  tl->SetTextSize (30);
  //  tl->DrawLatexNDC (0.34, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.34, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

  //  tl->SetTextSize (27);
  //  tl->SetTextAlign (21);
  //  tl->DrawLatexNDC (0.418, 0.74, "15-30");
  //  tl->DrawLatexNDC (0.518, 0.74, "30-60");
  //  tl->DrawLatexNDC (0.610, 0.74, "60+");

  //  tl->SetTextAlign (11);
  //  tl->DrawLatexNDC (0.66, 0.74, "#it{p}_{T}^{#it{Z}} [GeV]");
  //  tl->DrawLatexNDC (0.66, 0.685+0.007, "ATLAS 0-10\%");
  //  tl->DrawLatexNDC (0.66, 0.630+0.014, "JEWEL");
  //  myMarkerAndBoxAndLineText (0.654, 0.685+0.007, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.560, 0.685+0.007, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.8, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.464, 0.685+0.007, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.654, 0.630+0.014, 1.4, 1001, finalModelFillColors[3], 0.40, finalModelFillColors[3], -1, 1, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.560, 0.630+0.014, 1.4, 1001, finalModelFillColors[2], 0.40, finalModelFillColors[2], -1, 1, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.464, 0.630+0.014, 1.4, 1001, finalModelFillColors[1], 0.40, finalModelFillColors[1], -1, 1, "", 0.045);

  //  l->SetLineStyle (1);
  //  l->SetLineWidth (4);
  //  l->SetLineColor (finalModelFillColors[3]);
  //  l->DrawLineNDC (0.654 - (0.8*0.045) + 0.02 - (0.04*1.4), 0.630+0.014+(0.25*0.045), 0.654 - (0.8*0.045) + 0.02, 0.630+0.014+(0.25*0.045));
  //  l->SetLineColor (finalModelFillColors[2]);
  //  l->DrawLineNDC (0.560 - (0.8*0.045) + 0.02 - (0.04*1.4), 0.630+0.014+(0.25*0.045), 0.560 - (0.8*0.045) + 0.02, 0.630+0.014+(0.25*0.045));
  //  l->SetLineColor (finalModelFillColors[1]);
  //  l->DrawLineNDC (0.464 - (0.8*0.045) + 0.02 - (0.04*1.4), 0.630+0.014+(0.25*0.045), 0.464 - (0.8*0.045) + 0.02, 0.630+0.014+(0.25*0.045));

  //  c9->SaveAs (Form ("%s/FinalPlots/iaa_ptch_jewelComp.pdf", plotPath.Data ()));
  //}




  //{
  //  TCanvas* c10 = new TCanvas ("c10", "", 800, 800);

  //  const double lMargin = 0.15;
  //  const double rMargin = 0.04;
  //  const double bMargin = 0.15;
  //  const double tMargin = 0.04;

  //  c10->SetLeftMargin (lMargin);
  //  c10->SetRightMargin (rMargin);
  //  c10->SetBottomMargin (bMargin);
  //  c10->SetTopMargin (tMargin);

  //  c10->SetLogx ();

  //  short iPtZ = 4;
  //  const short iCent = 3;

  //  const double xmin = pTchBins[iPtZ][0];
  //  const double xmax = pTchBins[iPtZ][nPtchBins[iPtZ]];
  //  const double ymin = 0.;
  //  const double ymax = 3.6;

  //  {
  //    TH1D* h = new TH1D ("", "", 1, xmin, xmax);

  //    TAxis* xax = h->GetXaxis ();
  //    TAxis* yax = h->GetYaxis ();

  //    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
  //    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
  //    xax->SetRangeUser (xmin, xmax);
  //    xax->SetLabelSize (0);

  //    yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
  //    yax->SetTitleOffset (1.0 * yax->GetTitleOffset ());
  //    yax->SetLabelSize (0);
  //    yax->SetRangeUser (ymin, ymax);

  //    h->SetLineWidth (0);

  //    h->DrawCopy ("");
  //    SaferDelete (&h);

  //    tl->SetTextFont (43);
  //    tl->SetTextSize (36);
  //    tl->SetTextAlign (21);

  //    const double yoff = ymin - 0.05 * (ymax-ymin) / (1.-tMargin-bMargin);
  //    tl->DrawLatex (1,  yoff, "1");
  //    tl->DrawLatex (2,  yoff, "2");
  //    tl->DrawLatex (3,  yoff, "3");
  //    tl->DrawLatex (4,  yoff, "4");
  //    tl->DrawLatex (5,  yoff, "5");
  //    tl->DrawLatex (6,  yoff, "6");
  //    tl->DrawLatex (7,  yoff, "7");
  //    tl->DrawLatex (10, yoff, "10");
  //    tl->DrawLatex (20, yoff, "20");
  //    tl->DrawLatex (30, yoff, "30");
  //    tl->DrawLatex (40, yoff, "40");
  //    tl->DrawLatex (60, yoff, "60");

  //    tl->SetTextAlign (32);

  //    const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-lMargin-rMargin));
  //    tl->DrawLatex (xoff, 0, "0");
  //    tl->DrawLatex (xoff, 0.5, "0.5");
  //    tl->DrawLatex (xoff, 1, "1");
  //    tl->DrawLatex (xoff, 1.5, "1.5");
  //    tl->DrawLatex (xoff, 2, "2");
  //    tl->DrawLatex (xoff, 2.5, "2.5");
  //    tl->DrawLatex (xoff, 3, "3");
  //    tl->DrawLatex (xoff, 3.5, "3.5");

  //    l->SetLineStyle (2);
  //    l->SetLineWidth (2);
  //    l->SetLineColor (kBlack);
  //    l->DrawLine (xmin, 1, xmax, 1);
  //  }

  //  for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
  //    gPad->SetLogx ();

  //    TGAE* g = (TGAE*) g_scetg_pth[iPtZ]->Clone ();
  //    ResetXErrors (g);
  //    deltaize (g, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, xmin, xmax);
  //    ResetXErrors (g);

  //    g->SetFillColorAlpha (finalModelFillColors[iPtZ-1], 0.4);
  //    ((TGAE*) g->Clone ())->Draw ("3");
  //    ResetTGAEErrors (g);
  //    ResetXErrors (g);
  //    g->SetLineColor (finalModelFillColors[iPtZ-1]);
  //    g->SetLineWidth (4);
  //    ((TGAE*) g->Clone ())->Draw ("L");
  //    SaferDelete (&g);
  //  }

  //  for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
  //    TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

  //    RecenterGraph (g_syst);
  //    ResetXErrors (g_syst);
  //    deltaize (g_syst, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, xmin, xmax);
  //    ResetXErrors (g_syst);
  //    SetConstantXErrors (g_syst, 0.04, true, pTchBins[4][0], pTchBins[4][nPtchBins[4]]);

  //    g_syst->SetMarkerSize (0);
  //    g_syst->SetLineWidth (1);
  //    g_syst->SetLineColor (finalColors[iPtZ-1]);
  //    g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

  //    g_syst->Draw ("5P");
  //  }

  //  for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
  //    TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

  //    RecenterGraph (g_stat);
  //    ResetXErrors (g_stat);
  //    deltaize (g_stat, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, xmin, xmax);
  //    ResetXErrors (g_stat);

  //    Style_t markerStyle = markerStyles[iPtZ-2];
  //    float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

  //    g_stat->SetMarkerStyle (markerStyle);
  //    g_stat->SetMarkerSize (markerSize);
  //    g_stat->SetLineWidth (3);
  //    g_stat->SetMarkerColor (finalColors[iPtZ-1]);
  //    g_stat->SetLineColor (finalColors[iPtZ-1]);

  //    ((TGAE*) g_stat->Clone ())->Draw ("P");

  //    markerStyle = FullToOpenMarker (markerStyle);
  //    markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

  //    g_stat->SetMarkerStyle (markerStyle);
  //    g_stat->SetMarkerSize (markerSize);
  //    g_stat->SetLineWidth (0);
  //    g_stat->SetMarkerColor (kBlack);

  //    ((TGAE*) g_stat->Clone ())->Draw ("P");

  //    SaferDelete (&g_stat);
  //  }

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextAlign (11);

  //  tl->SetTextSize (32);
  //  tl->DrawLatexNDC (0.31, 0.890, "#bf{#it{ATLAS}} Internal");
  //  tl->SetTextSize (30);
  //  tl->DrawLatexNDC (0.34, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.34, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

  //  tl->SetTextSize (27);
  //  tl->SetTextAlign (21);
  //  tl->DrawLatexNDC (0.388, 0.74, "15-30");
  //  tl->DrawLatexNDC (0.488, 0.74, "30-60");
  //  tl->DrawLatexNDC (0.580, 0.74, "60+");

  //  tl->SetTextAlign (11);
  //  tl->DrawLatexNDC (0.63, 0.74, "#it{p}_{T}^{#it{Z}} [GeV]");
  //  tl->DrawLatexNDC (0.63, 0.685+0.007, "ATLAS 0-10\%");
  //  tl->DrawLatexNDC (0.63, 0.630+0.014, "SCET_{G} (#it{g}^{ }=^{ }2.0#pm0.2)");
  //  myMarkerAndBoxAndLineText (0.624, 0.685+0.007, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.530, 0.685+0.007, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.8, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.434, 0.685+0.007, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.624, 0.630+0.014, 1.4, 1001, finalModelFillColors[3], 0.40, finalModelFillColors[3], -1, 1, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.530, 0.630+0.014, 1.4, 1001, finalModelFillColors[2], 0.40, finalModelFillColors[2], -1, 1, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.434, 0.630+0.014, 1.4, 1001, finalModelFillColors[1], 0.40, finalModelFillColors[1], -1, 1, "", 0.045);

  //  l->SetLineStyle (1);
  //  l->SetLineWidth (4);
  //  l->SetLineColor (finalModelFillColors[3]);
  //  l->DrawLineNDC (0.624 - (0.8*0.045) + 0.02 - (0.04*1.4), 0.630+0.014+(0.25*0.045), 0.624 - (0.8*0.045) + 0.02, 0.630+0.014+(0.25*0.045));
  //  l->SetLineColor (finalModelFillColors[2]);
  //  l->DrawLineNDC (0.530 - (0.8*0.045) + 0.02 - (0.04*1.4), 0.630+0.014+(0.25*0.045), 0.530 - (0.8*0.045) + 0.02, 0.630+0.014+(0.25*0.045));
  //  l->SetLineColor (finalModelFillColors[1]);
  //  l->DrawLineNDC (0.434 - (0.8*0.045) + 0.02 - (0.04*1.4), 0.630+0.014+(0.25*0.045), 0.434 - (0.8*0.045) + 0.02, 0.630+0.014+(0.25*0.045));

  //  c10->SaveAs (Form ("%s/FinalPlots/iaa_ptch_scetgComp.pdf", plotPath.Data ()));
  //}





  //{
  //  ////////////////////////////////////////////////////////////////////////////////////////////////
  //  // Plots mean track <pTch> vs. <pT^Z>
  //  ////////////////////////////////////////////////////////////////////////////////////////////////
  //  TCanvas* c11 = new TCanvas ("c_mean_pTch_comb", "", 800, 800);
  //  c11->cd ();

  //  TH1D* h = new TH1D ("htemp", "", 1, 5, 100);

  //  TAxis* xax = h->GetXaxis ();
  //  TAxis* yax = h->GetYaxis ();

  //  xax->SetTitle ("#LT#it{p}_{T}^{#it{Z}}#GT [GeV]");
  //  xax->SetLabelFont (43);
  //  xax->SetLabelSize (34);
  //  xax->SetRangeUser (5, 100);
  //  xax->SetMoreLogLabels ();

  //  yax->SetTitle ("#LT#it{p}_{T}^{ch}#GT [GeV]");
  //  yax->SetLabelFont (43);
  //  yax->SetLabelSize (34);
  //  yax->SetRangeUser (0, 10);

  //  h->SetLineWidth (0);
  //  h->SetMarkerSize (0);

  //  h->DrawCopy ("hist ][");
  //  SaferDelete (&h);

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextAlign (11);
  //  tl->SetTextSize (32);
  //  tl->DrawLatexNDC (0.22, 0.87, "#bf{#it{ATLAS}} Internal");
  //  myMarkerAndBoxAndLineText (0.31, 0.820, 2.0, 1001, finalFillColors[0], 0.30, finalColors[0], kOpenCircle,     1.5, "#it{pp}", 0.036);
  //  myMarkerAndBoxAndLineText (0.31, 0.765, 2.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.5, "30-80\%", 0.036);
  //  myMarkerAndBoxAndLineText (0.56, 0.820, 2.0, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.5, "10-30\%", 0.036);
  //  myMarkerAndBoxAndLineText (0.56, 0.765, 2.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "0-10\%", 0.036);

  //  tl->SetTextSize (28);
  //  tl->DrawLatexNDC (0.22, 0.26, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");
  //  tl->DrawLatexNDC (0.22, 0.215, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.22, 0.715, "#it{p}_{T}^{ch} > 2 GeV");

  //  for (short iCent = 0; iCent < numCentBins; iCent++) {
  //    TGAE* g = (TGAE*) g_trk_avg_pt_ptz_syst[iCent]->Clone ();

  //    g->SetMarkerSize (0);
  //    //g->SetLineWidth (0);
  //    g->SetLineWidth (1);
  //    g->SetLineColor (finalColors[iCent]);
  //    g->SetFillColorAlpha (finalFillColors[iCent], 0.3);
  //    ((TGAE*)g->Clone ())->Draw ("5P");
  //    g->Draw ("2P");
  //  } // end loop over iCent

  //  for (short iCent = 0; iCent < numCentBins; iCent++) {
  //    TGAE* g = (TGAE*) g_trk_avg_pt_ptz_stat[iCent]->Clone ();

  //    Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
  //    float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

  //    g->SetMarkerStyle (markerStyle);
  //    g->SetMarkerSize (markerSize);
  //    g->SetMarkerColor (finalColors[iCent]);
  //    g->SetLineColor (finalColors[iCent]);
  //    g->SetLineWidth (3);
  //    ((TGAE*) g->Clone ())->Draw ("P");

  //    if (iCent != 0) {
  //      markerStyle = FullToOpenMarker (markerStyle);
  //      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
  //        
  //      g->SetMarkerStyle (markerStyle);
  //      g->SetMarkerSize (markerSize);
  //      g->SetLineWidth (0);
  //      g->SetMarkerColor (kBlack);

  //      ((TGAE*) g->Clone ())->Draw ("P");
  //    }
  //    else {
  //      TLine* g_line = new TLine ();
  //      g_line->SetLineWidth (2);
  //      g_line->SetLineColor (kBlack);
  //      for (int i = 0; i < g->GetN (); i++) {
  //        double x, y, yhi, ylo;
  //        g->GetPoint (i, x, y);
  //        yhi = y + g->GetErrorYhigh (i);
  //        ylo = y - g->GetErrorYlow (i);
  //        g_line->DrawLine (x, ylo, x, yhi);
  //      }
  //      SaferDelete (&g_line);
  //    }

  //    SaferDelete (&g);
  //  } // end loop over iCent

  //  c11->SaveAs (Form ("%s/FinalPlots/mean_ptch.pdf", plotPath.Data ()));
  //}




  //{
  //  TCanvas* c12 = new TCanvas ("c12", "", 800, 800);

  //  const double lMargin = 0.15;
  //  const double rMargin = 0.04;
  //  const double bMargin = 0.15;
  //  const double tMargin = 0.04;

  //  c12->SetLeftMargin (lMargin);
  //  c12->SetRightMargin (rMargin);
  //  c12->SetBottomMargin (bMargin);
  //  c12->SetTopMargin (tMargin);

  //  c12->SetLogx ();

  //  short iPtZ = 4;
  //  const short iCent = 3;

  //  const double xmin = pTchBins[iPtZ][0];
  //  const double xmax = pTchBins[iPtZ][nPtchBins[iPtZ]];
  //  const double ymin = 0.;
  //  const double ymax = 3.6;

  //  {
  //    TH1D* h = new TH1D ("", "", 1, xmin, xmax);

  //    TAxis* xax = h->GetXaxis ();
  //    TAxis* yax = h->GetYaxis ();

  //    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
  //    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
  //    xax->SetRangeUser (xmin, xmax);
  //    xax->SetLabelSize (0);

  //    yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
  //    yax->SetTitleOffset (1.0 * yax->GetTitleOffset ());
  //    yax->SetLabelSize (0);
  //    yax->SetRangeUser (ymin, ymax);

  //    h->SetLineWidth (0);

  //    h->DrawCopy ("");
  //    SaferDelete (&h);

  //    tl->SetTextFont (43);
  //    tl->SetTextSize (36);
  //    tl->SetTextAlign (21);

  //    const double yoff = ymin - 0.05 * (ymax-ymin) / (1.-tMargin-bMargin);
  //    tl->DrawLatex (1,  yoff, "1");
  //    tl->DrawLatex (2,  yoff, "2");
  //    tl->DrawLatex (3,  yoff, "3");
  //    tl->DrawLatex (4,  yoff, "4");
  //    tl->DrawLatex (5,  yoff, "5");
  //    tl->DrawLatex (6,  yoff, "6");
  //    tl->DrawLatex (7,  yoff, "7");
  //    tl->DrawLatex (10, yoff, "10");
  //    tl->DrawLatex (20, yoff, "20");
  //    tl->DrawLatex (30, yoff, "30");
  //    tl->DrawLatex (40, yoff, "40");
  //    tl->DrawLatex (60, yoff, "60");

  //    tl->SetTextAlign (32);

  //    const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-lMargin-rMargin));
  //    tl->DrawLatex (xoff, 0, "0");
  //    tl->DrawLatex (xoff, 0.5, "0.5");
  //    tl->DrawLatex (xoff, 1, "1");
  //    tl->DrawLatex (xoff, 1.5, "1.5");
  //    tl->DrawLatex (xoff, 2, "2");
  //    tl->DrawLatex (xoff, 2.5, "2.5");
  //    tl->DrawLatex (xoff, 3, "3");
  //    tl->DrawLatex (xoff, 3.5, "3.5");

  //    l->SetLineStyle (2);
  //    l->SetLineWidth (2);
  //    l->SetLineColor (kBlack);
  //    l->DrawLine (xmin, 1, xmax, 1);
  //  }

  //  for (iPtZ = nPtZBins-1; iPtZ >= 3; iPtZ--) {
  //    gPad->SetLogx ();

  //    TGAE* g = (TGAE*) g_colbt_pth[iPtZ]->Clone ();
  //    ResetXErrors (g);
  //    deltaize (g, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, xmin, xmax);
  //    ResetXErrors (g);

  //    g->SetFillColorAlpha (finalModelFillColors[iPtZ-1], 0.4);
  //    ((TGAE*) g->Clone ())->Draw ("3");
  //    ResetTGAEErrors (g);
  //    ResetXErrors (g);
  //    g->SetLineColor (finalModelFillColors[iPtZ-1]);
  //    g->SetLineWidth (4);
  //    ((TGAE*) g->Clone ())->Draw ("L");
  //    SaferDelete (&g);
  //  }

  //  for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
  //    TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

  //    RecenterGraph (g_syst);
  //    ResetXErrors (g_syst);
  //    deltaize (g_syst, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, xmin, xmax);
  //    ResetXErrors (g_syst);
  //    SetConstantXErrors (g_syst, 0.04, true, pTchBins[4][0], pTchBins[4][nPtchBins[4]]);

  //    g_syst->SetMarkerSize (0);
  //    g_syst->SetLineWidth (1);
  //    g_syst->SetLineColor (finalColors[iPtZ-1]);
  //    g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

  //    g_syst->Draw ("5P");
  //  }

  //  for (iPtZ = nPtZBins-1; iPtZ >= 2; iPtZ--) {
  //    TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

  //    RecenterGraph (g_stat);
  //    ResetXErrors (g_stat);
  //    deltaize (g_stat, 0.06 * (iPtZ-2 - 0.5*(nPtZBins-3)), true, xmin, xmax);
  //    ResetXErrors (g_stat);

  //    Style_t markerStyle = markerStyles[iPtZ-2];
  //    float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

  //    g_stat->SetMarkerStyle (markerStyle);
  //    g_stat->SetMarkerSize (markerSize);
  //    g_stat->SetLineWidth (3);
  //    g_stat->SetMarkerColor (finalColors[iPtZ-1]);
  //    g_stat->SetLineColor (finalColors[iPtZ-1]);

  //    ((TGAE*) g_stat->Clone ())->Draw ("P");

  //    markerStyle = FullToOpenMarker (markerStyle);
  //    markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

  //    g_stat->SetMarkerStyle (markerStyle);
  //    g_stat->SetMarkerSize (markerSize);
  //    g_stat->SetLineWidth (0);
  //    g_stat->SetMarkerColor (kBlack);

  //    ((TGAE*) g_stat->Clone ())->Draw ("P");

  //    SaferDelete (&g_stat);
  //  }

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextAlign (11);

  //  tl->SetTextSize (32);
  //  tl->DrawLatexNDC (0.31, 0.890, "#bf{#it{ATLAS}} Internal");
  //  tl->SetTextSize (30);
  //  tl->DrawLatexNDC (0.34, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.34, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

  //  tl->SetTextSize (27);
  //  tl->SetTextAlign (21);
  //  //tl->DrawLatexNDC (0.388, 0.74, "15-30");
  //  tl->DrawLatexNDC (0.488, 0.74, "30-60");
  //  tl->DrawLatexNDC (0.580, 0.74, "60+");

  //  tl->SetTextAlign (11);
  //  tl->DrawLatexNDC (0.63, 0.74, "#it{p}_{T}^{#it{Z}} [GeV]");
  //  tl->DrawLatexNDC (0.63, 0.685+0.007, "ATLAS 0-10\%");
  //  tl->DrawLatexNDC (0.63, 0.630+0.014, "CoLBT-hydro");
  //  myMarkerAndBoxAndLineText (0.624, 0.685+0.007, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.530, 0.685+0.007, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.8, "", 0.045);
  //  //myMarkerAndBoxAndLineText (0.434, 0.685+0.007, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.624, 0.630+0.014, 1.4, 1001, finalModelFillColors[3], 0.40, finalModelFillColors[3], -1, 1, "", 0.045);
  //  myMarkerAndBoxAndLineText (0.530, 0.630+0.014, 1.4, 1001, finalModelFillColors[2], 0.40, finalModelFillColors[2], -1, 1, "", 0.045);
  //  //myMarkerAndBoxAndLineText (0.434, 0.630+0.014, 1.4, 1001, finalModelFillColors[1], 0.40, finalModelFillColors[1], -1, 1, "", 0.045);

  //  l->SetLineStyle (1);
  //  l->SetLineWidth (4);
  //  l->SetLineColor (finalModelFillColors[3]);
  //  l->DrawLineNDC (0.624 - (0.8*0.045) + 0.02 - (0.04*1.4), 0.630+0.014+(0.25*0.045), 0.624 - (0.8*0.045) + 0.02, 0.630+0.014+(0.25*0.045));
  //  l->SetLineColor (finalModelFillColors[2]);
  //  l->DrawLineNDC (0.530 - (0.8*0.045) + 0.02 - (0.04*1.4), 0.630+0.014+(0.25*0.045), 0.530 - (0.8*0.045) + 0.02, 0.630+0.014+(0.25*0.045));
  //  //l->SetLineColor (finalModelFillColors[1]);
  //  //l->DrawLineNDC (0.434 - (0.8*0.045) + 0.02 - (0.04*1.4), 0.630+0.014+(0.25*0.045), 0.434 - (0.8*0.045) + 0.02, 0.630+0.014+(0.25*0.045));

  //  c12->SaveAs (Form ("%s/FinalPlots/iaa_ptch_colbtComp.pdf", plotPath.Data ()));
  //}




  //{
  //  ////////////////////////////////////////////////////////////////////////////////////////////////
  //  // Plots mean number of charged particles opposite the Z
  //  ////////////////////////////////////////////////////////////////////////////////////////////////
  //  TCanvas* c13 = new TCanvas ("c_mean_ntrk_comb", "", 800, 800);
  //  c13->cd ();

  //  TH1D* h = new TH1D ("htemp", "", 1, 5, 100);

  //  TAxis* xax = h->GetXaxis ();
  //  TAxis* yax = h->GetYaxis ();

  //  xax->SetTitle ("#LT#it{p}_{T}^{#it{Z}}#GT [GeV]");
  //  xax->SetLabelFont (43);
  //  xax->SetLabelSize (34);
  //  xax->SetRangeUser (5, 100);
  //  xax->SetMoreLogLabels ();

  //  yax->SetTitle ("#LTN_{ch}#GT");
  //  yax->SetLabelFont (43);
  //  yax->SetLabelSize (34);
  //  yax->SetRangeUser (0, 8);

  //  h->SetLineWidth (0);
  //  h->SetMarkerSize (0);

  //  h->DrawCopy ("hist ][");
  //  SaferDelete (&h);

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextAlign (11);
  //  tl->SetTextSize (32);
  //  tl->DrawLatexNDC (0.22, 0.87, "#bf{#it{ATLAS}} Internal");
  //  myMarkerAndBoxAndLineText (0.31, 0.710, 2.0, 1001, finalFillColors[0], 0.30, finalColors[0], kOpenCircle,     1.5, "#it{pp}", 0.036);
  //  myMarkerAndBoxAndLineText (0.31, 0.655, 2.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.5, "30-80\%", 0.036);
  //  myMarkerAndBoxAndLineText (0.56, 0.710, 2.0, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.5, "10-30\%", 0.036);
  //  myMarkerAndBoxAndLineText (0.56, 0.655, 2.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "0-10\%", 0.036);

  //  tl->SetTextSize (28);
  //  tl->DrawLatexNDC (0.22, 0.815, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");
  //  tl->DrawLatexNDC (0.22, 0.765, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.22, 0.605, "#it{p}_{T}^{ch} > 2 GeV");

  //  for (short iCent = 0; iCent < numCentBins; iCent++) {
  //    TGAE* g = (TGAE*) g_avg_ntrk_ptz_syst[iCent]->Clone ();

  //    g->SetMarkerSize (0);
  //    //g->SetLineWidth (0);
  //    g->SetLineWidth (1);
  //    g->SetLineColor (finalColors[iCent]);
  //    g->SetFillColorAlpha (finalFillColors[iCent], 0.3);
  //    ((TGAE*)g->Clone ())->Draw ("5P");
  //    g->Draw ("2P");
  //  } // end loop over iCent

  //  for (short iCent = 0; iCent < numCentBins; iCent++) {
  //    TGAE* g = (TGAE*) g_avg_ntrk_ptz_stat[iCent]->Clone ();

  //    Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
  //    float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

  //    g->SetMarkerStyle (markerStyle);
  //    g->SetMarkerSize (markerSize);
  //    g->SetMarkerColor (finalColors[iCent]);
  //    g->SetLineColor (finalColors[iCent]);
  //    g->SetLineWidth (3);
  //    ((TGAE*) g->Clone ())->Draw ("P");

  //    if (iCent != 0) {
  //      markerStyle = FullToOpenMarker (markerStyle);
  //      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
  //        
  //      g->SetMarkerStyle (markerStyle);
  //      g->SetMarkerSize (markerSize);
  //      g->SetLineWidth (0);
  //      g->SetMarkerColor (kBlack);

  //      ((TGAE*) g->Clone ())->Draw ("P");
  //    }
  //    else {
  //      TLine* g_line = new TLine ();
  //      g_line->SetLineWidth (2);
  //      g_line->SetLineColor (kBlack);
  //      for (int i = 0; i < g->GetN (); i++) {
  //        double x, y, yhi, ylo;
  //        g->GetPoint (i, x, y);
  //        yhi = y + g->GetErrorYhigh (i);
  //        ylo = y - g->GetErrorYlow (i);
  //        g_line->DrawLine (x, ylo, x, yhi);
  //      }
  //      SaferDelete (&g_line);
  //    }

  //    SaferDelete (&g);
  //  } // end loop over iCent

  //  c13->SaveAs (Form ("%s/FinalPlots/mean_ntrk.pdf", plotPath.Data ()));
  //}





  //{
  //  ////////////////////////////////////////////////////////////////////////////////////////////////
  //  // Plots mean track <xhZ> vs. <pT^Z>
  //  ////////////////////////////////////////////////////////////////////////////////////////////////
  //  TCanvas* c14 = new TCanvas ("c_mean_xhZ_comb", "", 800, 800);
  //  c14->cd ();

  //  TH1D* h = new TH1D ("htemp", "", 1, 5, 100);

  //  TAxis* xax = h->GetXaxis ();
  //  TAxis* yax = h->GetYaxis ();

  //  xax->SetTitle ("#LT#it{p}_{T}^{#it{Z}}#GT [GeV]");
  //  xax->SetLabelFont (43);
  //  xax->SetLabelSize (34);
  //  xax->SetRangeUser (5, 100);
  //  xax->SetMoreLogLabels ();

  //  yax->SetTitle ("#LT#it{x}_{h#it{Z}}#GT");
  //  yax->SetLabelFont (43);
  //  yax->SetLabelSize (34);
  //  yax->SetRangeUser (0, 0.30);

  //  h->SetLineWidth (0);
  //  h->SetMarkerSize (0);

  //  h->DrawCopy ("hist ][");
  //  SaferDelete (&h);

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextAlign (11);
  //  tl->SetTextSize (32);
  //  tl->DrawLatexNDC (0.22, 0.87, "#bf{#it{ATLAS}} Internal");
  //  myMarkerAndBoxAndLineText (0.31, 0.820, 2.0, 1001, finalFillColors[0], 0.30, finalColors[0], kOpenCircle,     1.8, "#it{pp}", 0.036);
  //  myMarkerAndBoxAndLineText (0.31, 0.765, 2.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "30-80\%", 0.036);
  //  myMarkerAndBoxAndLineText (0.56, 0.820, 2.0, 1001, finalFillColors[2], 0.30, finalColors[2], markerStyles[1], 1.8, "10-30\%", 0.036);
  //  myMarkerAndBoxAndLineText (0.56, 0.765, 2.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "0-10\%", 0.036);

  //  tl->SetTextSize (28);
  //  tl->DrawLatexNDC (0.22, 0.26, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");
  //  tl->DrawLatexNDC (0.22, 0.215, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.22, 0.715, "#it{x}_{h#it{Z}} > 1/15");

  //  for (short iCent = 0; iCent < numCentBins; iCent++) {
  //    TGAE* g = (TGAE*) g_trk_avg_xhz_ptz_syst[iCent]->Clone ();

  //    g->SetMarkerSize (0);
  //    //g->SetLineWidth (0);
  //    g->SetLineWidth (1);
  //    g->SetLineColor (finalColors[iCent]);
  //    g->SetFillColorAlpha (finalFillColors[iCent], 0.3);
  //    ((TGAE*)g->Clone ())->Draw ("5P");
  //    g->Draw ("2P");
  //  } // end loop over iCent

  //  for (short iCent = 0; iCent < numCentBins; iCent++) {
  //    TGAE* g = (TGAE*) g_trk_avg_xhz_ptz_stat[iCent]->Clone ();

  //    Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iCent-1]);
  //    float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

  //    g->SetMarkerStyle (markerStyle);
  //    g->SetMarkerSize (markerSize);
  //    g->SetMarkerColor (finalColors[iCent]);
  //    g->SetLineColor (finalColors[iCent]);
  //    g->SetLineWidth (3);
  //    ((TGAE*) g->Clone ())->Draw ("P");

  //    if (iCent != 0) {
  //      markerStyle = FullToOpenMarker (markerStyle);
  //      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
  //        
  //      g->SetMarkerStyle (markerStyle);
  //      g->SetMarkerSize (markerSize);
  //      g->SetLineWidth (0);
  //      g->SetMarkerColor (kBlack);

  //      ((TGAE*) g->Clone ())->Draw ("P");
  //    }
  //    else {
  //      TLine* g_line = new TLine ();
  //      g_line->SetLineWidth (2);
  //      g_line->SetLineColor (kBlack);
  //      for (int i = 0; i < g->GetN (); i++) {
  //        double x, y, yhi, ylo;
  //        g->GetPoint (i, x, y);
  //        yhi = y + g->GetErrorYhigh (i);
  //        ylo = y - g->GetErrorYlow (i);
  //        g_line->DrawLine (x, ylo, x, yhi);
  //      }
  //      SaferDelete (&g_line);
  //    }

  //    SaferDelete (&g);
  //  } // end loop over iCent

  //  c14->SaveAs (Form ("%s/FinalPlots/mean_xhz.pdf", plotPath.Data ()));
  //}




  //{
  //  const char* canvasName = "c15";
  //  TCanvas* c15 = new TCanvas (canvasName, "", 800, 880);
  //  c15->Draw ();

  //  const double lMargin = 0.15;
  //  const double rMargin = 0.04;
  //  const double utMargin = 0.04;
  //  const double ubMargin = 0;//0.01;
  //  const double dtMargin = 0;//0.02;
  //  const double dbMargin = 0.35;

  //  TPad* uPad = new TPad (Form ("%s_uPad", canvasName), "", 0, 0.3, 1, 1);
  //  TPad* dPad = new TPad (Form ("%s_dPad", canvasName), "", 0, 0, 1, 0.3);

  //  uPad->SetLeftMargin (lMargin);
  //  uPad->SetRightMargin (rMargin);
  //  dPad->SetLeftMargin (lMargin);
  //  dPad->SetRightMargin (rMargin);
  //  uPad->SetTopMargin (utMargin);
  //  uPad->SetBottomMargin (ubMargin);
  //  dPad->SetTopMargin (dtMargin);
  //  dPad->SetBottomMargin (dbMargin);

  //  uPad->Draw ();
  //  dPad->Draw ();

  //  uPad->cd ();
  //  uPad->SetLogx ();
  //  uPad->SetLogy ();

  //  const double xmin = pTchBins[nPtZBins-1][0];
  //  const double xmax =  pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]];

  //  {
  //    TH1D* h = new TH1D ("", "", 1, xmin, xmax);

  //    TAxis* xax = h->GetXaxis ();
  //    TAxis* yax = h->GetYaxis ();

  //    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
  //    xax->SetTitleSize (0);
  //    xax->SetTickLength (0.03 * (1.-utMargin-ubMargin) / uPad->GetHNDC ());
  //    xax->SetTickLength (0.03);
  //    xax->SetLabelSize (0);

  //    xax->SetRangeUser (xmin, xmax);

  //    yax->SetTitle ("(1/N_{Z}) (d^{2}N_{ch} / d#it{p}_{T} d#Delta#phi) [GeV^{-1}]");
  //    yax->SetTitleFont (43);
  //    yax->SetTitleSize (36);
  //    yax->SetTickLength (0.02 * (1.-utMargin-ubMargin) / uPad->GetHNDC ());
  //    yax->SetLabelFont (43);
  //    yax->SetLabelSize (32);
  //    double ymin = 2e-3;
  //    double ymax = 8e3;
  //    yax->SetRangeUser (ymin, ymax);

  //    h->SetLineWidth (0);

  //    h->DrawCopy ("");
  //    SaferDelete (&h);
  //  }

  //  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
  //    TGAE* g = (TGAE*) g_pythia_finepth_ptz[iPtZ]->Clone ();
  //    SetMinErrors (g, 0.10, true);

  //    OffsetYAxis (g, pow (10, iPtZ-2), true);
  //    RecenterGraph (g);
  //    //ResetXErrors (g);

  //    //g->SetMarkerSize (0);
  //    //g->SetLineWidth (1);
  //    //g->SetLineColor (finalColors[iPtZ-1]);
  //    //g->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

  //    //((TGAE*) g->Clone ())->Draw ("3");
  //    //SaferDelete (&g);

  //    g->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);
  //    ((TGAE*) g->Clone ())->Draw ("3");
  //    ResetTGAEErrors (g);
  //    ResetXErrors (g);
  //    g->SetLineColor (finalFillColors[iPtZ-1]);
  //    g->SetLineWidth (2);
  //    ((TGAE*) g->Clone ())->Draw ("L");
  //    SaferDelete (&g);
  //  }

  //  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
  //    const int iCent = 0;
  //    TGAE* g_syst = (TGAE*) g_trk_pt_ptz_sub_syst[iPtZ][iCent]->Clone ();

  //    OffsetYAxis (g_syst, pow (10, iPtZ-2), true);
  //    RecenterGraph (g_syst);
  //    //ResetXErrors (g_syst);
  //    //deltaize (g_syst, 0.05*(0.5*(numCentBins-1)-iCent), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
  //    //ResetXErrors (g_syst);
  //    //SetConstantXErrors (g_syst, 0.05, true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);

  //    g_syst->SetMarkerSize (0);
  //    g_syst->SetLineWidth (1);
  //    g_syst->SetMarkerColor (finalColors[iPtZ-1]);
  //    g_syst->SetLineColor (finalColors[iPtZ-1]);
  //    g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

  //    ((TGAE*) g_syst->Clone ())->Draw ("5P");

  //    SaferDelete (&g_syst);

  //    TGAE* g_stat = make_graph (h_trk_pt_ptz_sub_stat[iPtZ][iCent]);

  //    Style_t markerStyle = FullToOpenMarker (markerStyles[iPtZ-2]);
  //    float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);

  //    OffsetYAxis (g_stat, pow (10, iPtZ-2), true);
  //    RecenterGraph (g_stat);
  //    ResetXErrors (g_stat);
  //    //deltaize (g_stat, 0.05*(0.5*(numCentBins-1)-iCent), true, pTchBins[nPtZBins-1][0], pTchBins[nPtZBins-1][nPtchBins[nPtZBins-1]]);
  //    //ResetXErrors (g_stat);

  //    g_stat->SetMarkerStyle (markerStyle);
  //    g_stat->SetMarkerSize (markerSize);
  //    g_stat->SetLineWidth (3);
  //    g_stat->SetMarkerColor (kBlack);
  //    g_stat->SetLineColor (finalColors[iPtZ-1]);

  //    ((TGAE*) g_stat->Clone ())->Draw ("P");

  //    TLine* g_stat_line = new TLine ();
  //    g_stat_line->SetLineWidth (2);
  //    g_stat_line->SetLineColor (finalColors[iPtZ-1]);
  //    for (int i = 0; i < g_stat->GetN (); i++) {
  //      double x, y, yhi, ylo;
  //      g_stat->GetPoint (i, x, y);
  //      yhi = y + g_stat->GetErrorYhigh (i);
  //      ylo = y - g_stat->GetErrorYlow (i);
  //      g_stat_line->DrawLine (x, ylo, x, yhi);
  //    }

  //    //markerStyle = kDot;
  //    
  //    //g_stat->SetMarkerStyle (markerStyle);
  //    //g_stat->SetMarkerSize (0);
  //    //g_stat->SetMarkerColor (kBlack);

  //    //((TGAE*) g_stat->Clone ())->Draw ("P");

  //    SaferDelete (&g_stat_line);

  //    SaferDelete (&g_stat);
  //  } // end loop over iPtZ

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextAlign (11);

  //  tl->SetTextSize (28);
  //  tl->DrawLatexNDC (0.20, 0.14, "#bf{#it{ATLAS}} Internal");
  //  tl->SetTextSize (28);
  //  tl->DrawLatexNDC (0.20, 0.09, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.20, 0.04, "Powheg + Pythia 8.186");

  //  tl->SetTextSize (22);
  //  tl->DrawLatexNDC (0.730, 0.885, "#it{p}_{T}^{#it{Z}} [GeV]");
  //  tl->DrawLatexNDC (0.730, 0.835, "60+ (#times 10^{2})");
  //  tl->DrawLatexNDC (0.730, 0.785, "30-60 (#times 10)");
  //  tl->DrawLatexNDC (0.730, 0.735, "15-30 (#times 1)");
  //  myMarkerAndBoxAndLineText (0.72, 0.840-0.002, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], -1, 1.6, "", 0.036);
  //  myMarkerAndBoxAndLineText (0.63, 0.840-0.002, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], FullToOpenMarker (markerStyles[2]), 2.2, "", 0.036);
  //  myMarkerAndBoxAndLineText (0.72, 0.790-0.002, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], -1, 1.6, "", 0.036);
  //  myMarkerAndBoxAndLineText (0.63, 0.790-0.002, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], FullToOpenMarker (markerStyles[1]), 1.6, "", 0.036);
  //  myMarkerAndBoxAndLineText (0.72, 0.740-0.002, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], -1, 1.6, "", 0.036);
  //  myMarkerAndBoxAndLineText (0.63, 0.740-0.002, 1.4, 1001, finalFillColors[1], 0.30, finalColors[1], FullToOpenMarker (markerStyles[0]), 1.6, "", 0.036);

  //  l->SetLineStyle (1);
  //  l->SetLineWidth (2);
  //  l->SetLineColor (finalModelFillColors[3]);
  //  l->DrawLineNDC (0.72 - (0.8*0.036) + 0.02 - (0.04*1.4), 0.840-0.002+(0.25*0.036), 0.72 - (0.8*0.036) + 0.02, 0.840-0.002+(0.25*0.036));
  //  l->SetLineColor (finalModelFillColors[2]);
  //  l->DrawLineNDC (0.72 - (0.8*0.036) + 0.02 - (0.04*1.4), 0.790-0.002+(0.25*0.036), 0.72 - (0.8*0.036) + 0.02, 0.790-0.002+(0.25*0.036));
  //  l->SetLineColor (finalModelFillColors[1]);
  //  l->DrawLineNDC (0.72 - (0.8*0.036) + 0.02 - (0.04*1.4), 0.740-0.002+(0.25*0.036), 0.72 - (0.8*0.036) + 0.02, 0.740-0.002+(0.25*0.036));

  //  tl->SetTextSize (22);
  //  tl->DrawLatexNDC (0.565, 0.885, "Data");
  //  tl->SetTextSize (22);
  //  tl->DrawLatexNDC (0.663, 0.885, "MC");


  //  dPad->cd ();
  //  dPad->SetLogx ();

  //  {
  //    TH1D* h = new TH1D ("", "", 1, xmin, xmax);

  //    TAxis* xax = h->GetXaxis ();
  //    TAxis* yax = h->GetYaxis ();

  //    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
  //    xax->SetTitleFont (43);
  //    xax->SetTitleSize (36);
  //    xax->SetTitleOffset (3.6);
  //    xax->SetTickLength (0.03 * (1.-dtMargin-dbMargin) / dPad->GetHNDC ());
  //    xax->SetLabelSize (0);
  //    xax->SetRangeUser (xmin, xmax);

  //    yax->SetTitle ("MC / data");
  //    yax->SetTitleFont (43);
  //    yax->SetTitleSize (36);
  //    yax->CenterTitle (true);
  //    yax->SetTickLength (0.02 * (1.-dtMargin-dbMargin) / dPad->GetHNDC ());
  //    yax->SetLabelFont (43);
  //    yax->SetLabelSize (32);
  //    double ymin = 0.88;
  //    double ymax = 1.12;
  //    yax->SetRangeUser (ymin, ymax);
  //    yax->SetNdivisions (504);

  //    h->SetLineWidth (0);

  //    h->DrawCopy ("");
  //    SaferDelete (&h);

  //    tl->SetTextFont (43);
  //    tl->SetTextSize (32);
  //    tl->SetTextAlign (21);
  //    const double yoff = ymin - (0.12 * (ymax - ymin) / (1.-dtMargin-dbMargin));
  //    tl->DrawLatex (1,  yoff, "1");
  //    tl->DrawLatex (2,  yoff, "2");
  //    tl->DrawLatex (3,  yoff, "3");
  //    tl->DrawLatex (4,  yoff, "4");
  //    tl->DrawLatex (5,  yoff, "5");
  //    tl->DrawLatex (6,  yoff, "6");
  //    tl->DrawLatex (7,  yoff, "7");
  //    tl->DrawLatex (10, yoff, "10");
  //    tl->DrawLatex (20, yoff, "20");
  //    tl->DrawLatex (30, yoff, "30");
  //    tl->DrawLatex (40, yoff, "40");
  //    tl->DrawLatex (60, yoff, "60");

  //    TLine* line = new TLine (xmin, 1, xmax, 1);
  //    line->SetLineColor (kBlack);
  //    line->SetLineStyle (2);
  //    line->Draw ("same");
  //  }

  //  for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
  //    const int iCent = 0;
  //    TGAE* g_syst = (TGAE*) g_trk_pt_ptz_sub_syst[iPtZ][iCent]->Clone ();
  //    TGAE* g_stat = make_graph (h_trk_pt_ptz_sub_stat[iPtZ][iCent]);
  //    TGAE* g_pyth = (TGAE*) g_pythia_pth_ptz[iPtZ]->Clone ();

  //    RecenterGraph (g_syst);
  //    RecenterGraph (g_stat);
  //    RecenterGraph (g_pyth);

  //    TGAE* g_ratio_stat = new TGAE ();
  //    TGAE* g_ratio_syst = new TGAE ();

  //    double x, y_data, y_pyth;
  //    double y_stat_hi, y_stat_lo, y_syst_hi, y_syst_lo, y_pyth_hi, y_pyth_lo, x_syst_hi, x_syst_lo;
  //    for (int iX = 0; iX < g_syst->GetN (); iX++) {
  //      g_syst->GetPoint (iX, x, y_data);
  //      y_syst_hi = g_syst->GetErrorYhigh (iX);
  //      y_syst_lo = g_syst->GetErrorYlow (iX);
  //      x_syst_hi = g_syst->GetErrorXhigh (iX);
  //      x_syst_lo = g_syst->GetErrorXlow (iX);
  //      g_stat->GetPoint (iX, x, y_data);
  //      y_stat_hi = g_stat->GetErrorYhigh (iX);
  //      y_stat_lo = g_stat->GetErrorYlow (iX);
  //      g_pyth->GetPoint (iX, x, y_pyth);
  //      y_pyth_hi = g_pyth->GetErrorYhigh (iX);
  //      y_pyth_lo = g_pyth->GetErrorYlow (iX);

  //      const double ratio = y_pyth / y_data;

  //      g_ratio_syst->SetPoint (g_ratio_syst->GetN (), x, ratio);
  //      g_ratio_syst->SetPointEXhigh (g_ratio_syst->GetN () - 1, x_syst_hi);
  //      g_ratio_syst->SetPointEXlow (g_ratio_syst->GetN () - 1, x_syst_lo);
  //      g_ratio_syst->SetPointEYhigh (g_ratio_syst->GetN () - 1, fabs (ratio) * fabs (y_syst_hi / y_data));
  //      g_ratio_syst->SetPointEYlow (g_ratio_syst->GetN () - 1, fabs (ratio) * fabs (y_syst_lo / y_data));

  //      g_ratio_stat->SetPoint (g_ratio_stat->GetN (), x, ratio);
  //      g_ratio_stat->SetPointEXhigh (g_ratio_stat->GetN () - 1, 0);
  //      g_ratio_stat->SetPointEXlow (g_ratio_stat->GetN () - 1, 0);
  //      g_ratio_stat->SetPointEYhigh (g_ratio_stat->GetN () - 1, fabs (ratio) * sqrt (pow (y_pyth_hi / y_pyth, 2) + pow (y_stat_hi / y_data, 2)));
  //      g_ratio_stat->SetPointEYlow (g_ratio_stat->GetN () - 1, fabs (ratio) * sqrt (pow (y_pyth_lo / y_pyth, 2) + pow (y_stat_lo / y_data, 2)));
  //    }

  //    g_ratio_syst->SetMarkerSize (0);
  //    g_ratio_syst->SetLineWidth (1);
  //    g_ratio_syst->SetMarkerColor (finalColors[iPtZ-1]);
  //    g_ratio_syst->SetLineColor (finalColors[iPtZ-1]);
  //    g_ratio_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

  //    ((TGAE*) g_ratio_syst->Clone ())->Draw ("5P");

  //    SaferDelete (&g_ratio_syst);

  //    Style_t markerStyle = FullToOpenMarker (markerStyles[iPtZ-2]);
  //    float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.2 : 1.6);

  //    g_ratio_stat->SetMarkerStyle (markerStyle);
  //    g_ratio_stat->SetMarkerSize (markerSize);
  //    g_ratio_stat->SetLineWidth (3);
  //    g_ratio_stat->SetMarkerColor (kBlack);
  //    g_ratio_stat->SetLineColor (finalColors[iPtZ-1]);

  //    ((TGAE*) g_ratio_stat->Clone ())->Draw ("P");

  //    TLine* g_ratio_stat_line = new TLine ();
  //    g_ratio_stat_line->SetLineWidth (2);
  //    g_ratio_stat_line->SetLineColor (finalColors[iPtZ-1]);
  //    for (int i = 0; i < g_stat->GetN (); i++) {
  //      double x, y, yhi, ylo;
  //      g_ratio_stat->GetPoint (i, x, y);
  //      yhi = y + g_ratio_stat->GetErrorYhigh (i);
  //      ylo = y - g_ratio_stat->GetErrorYlow (i);
  //      g_ratio_stat_line->DrawLine (x, ylo, x, yhi);
  //    }

  //    SaferDelete (&g_ratio_stat);
  //    SaferDelete (&g_ratio_stat_line);

  //    SaferDelete (&g_syst);
  //    SaferDelete (&g_stat);
  //    SaferDelete (&g_pyth);
  //  }

  //  c15->SaveAs (Form ("%s/FinalPlots/yield_allptz_pTch_pythiaComp_onePlot.pdf", plotPath.Data ()));
  //}




  {
    const char* canvasName = "c16";
    TCanvas* c16 = new TCanvas (canvasName, "", 2000, 765);

    const double llMargin = 0.12;
    const double lrMargin = 0.0;
    const double clMargin = 0.0;
    const double crMargin = 0.0;
    const double rlMargin = 0.0;
    const double rrMargin = 0.040;
    const double bMargin = 0.12;
    const double tMargin = 0.04;

    const double deltaL = (1. - llMargin - lrMargin);
    const double deltaC = (1. - clMargin - crMargin);
    const double deltaR = (1. - rlMargin - rrMargin);

    const double a = (double) (deltaR * deltaC / (deltaL*deltaR + deltaC*deltaR + deltaL*deltaC));
    const double b = (double) (deltaR * deltaL / (deltaL*deltaR + deltaC*deltaR + deltaL*deltaC));

    const double xPadLCMiddle = a;
    const double xPadCRMiddle = a+b;

    TPad* lPad = nullptr;
    TPad* cPad = nullptr;
    TPad* rPad = nullptr;

    lPad = new TPad (Form ("%s_lPad", canvasName), "", 0, 0, xPadLCMiddle, 1);
    cPad = new TPad (Form ("%s_cPad", canvasName), "", xPadLCMiddle, 0, xPadCRMiddle, 1);
    rPad = new TPad (Form ("%s_rPad", canvasName), "", xPadCRMiddle, 0, 1, 1);

    lPad->SetLeftMargin (llMargin);
    lPad->SetRightMargin (lrMargin);
    cPad->SetLeftMargin (clMargin);
    cPad->SetRightMargin (crMargin);
    rPad->SetLeftMargin (rlMargin);
    rPad->SetRightMargin (rrMargin);
    lPad->SetBottomMargin (bMargin);
    lPad->SetTopMargin (tMargin);
    cPad->SetBottomMargin (bMargin);
    cPad->SetTopMargin (tMargin);
    rPad->SetBottomMargin (bMargin);
    rPad->SetTopMargin (tMargin);

    TPad* pads[3] = {lPad, cPad, rPad};

    for (int i = 0; i < 3; i++) {
      pads[i]->Draw ();
      pads[i]->SetLogx ();
      pads[i]->SetLogy ();
    }

    const short iCent = numCentBins-1;

    for (short iPtZ = 2; iPtZ < 5; iPtZ++) {
      pads[iPtZ-2]->cd ();

      const double xmin = pTchBins[iPtZ][0];
      const double xmax = pTchBins[iPtZ][nPtchBins[iPtZ]];

      TH1D* h = new TH1D ("", "", 1, xmin, xmax);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetTitleOffset (0.8 * xax->GetTitleOffset ());
      xax->SetTitleFont (43);
      xax->SetTitleSize (36);
      xax->SetRangeUser (xmin, xmax);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
      yax->SetTitleOffset (0.8 * yax->GetTitleOffset ());
      yax->SetTitleFont (43);
      yax->SetTitleSize (36);
      const double ymin = 0.18;
      const double ymax = 7;
      yax->SetRangeUser (ymin, ymax);
      yax->SetLabelSize (0);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);
      const double yoff = ymin / exp (0.044 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      if (iPtZ == 2)
        tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (3,  yoff, "3");
      tl->DrawLatex (4,  yoff, "4");
      tl->DrawLatex (5,  yoff, "5");
      tl->DrawLatex (6,  yoff, "6");
      tl->DrawLatex (7,  yoff, "7");
      tl->DrawLatex (10, yoff, "10");
      if (iPtZ > 2) {
        tl->DrawLatex (20, yoff, "20");
      }
      if (iPtZ > 3) {
        tl->DrawLatex (30, yoff, "30");
        tl->DrawLatex (40, yoff, "40");
      }

      if (iPtZ == 2) {
        tl->SetTextAlign (32);

        const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-gPad->GetLeftMargin()-gPad->GetRightMargin()));
        tl->DrawLatex (xoff, 0.2, "0.2");
        tl->DrawLatex (xoff, 0.3, "0.3");
        tl->DrawLatex (xoff, 0.4, "0.4");
        tl->DrawLatex (xoff, 0.5, "0.5");
        tl->DrawLatex (xoff, 0.7, "0.7");
        tl->DrawLatex (xoff, 1, "1");
        tl->DrawLatex (xoff, 2, "2");
        tl->DrawLatex (xoff, 3, "3");
        tl->DrawLatex (xoff, 4, "4");
        tl->DrawLatex (xoff, 5, "5");
      }

      {
        TGAE* g = (TGAE*) g_hybrid_pth[iPtZ]->Clone ();
        RecenterGraph (g);
        SetMinErrors (g, minModelUnc, true);

        g->SetFillColorAlpha (hybridColor, hybridAlpha);
        ((TGAE*) g->Clone ())->Draw ("3");
        ResetTGAEErrors (g);
        ResetXErrors (g);
        g->SetLineColor (hybridColor);
        g->SetLineWidth (4);
        ((TGAE*) g->Clone ())->Draw ("L");
        SaferDelete (&g);
      }

      {
        TGAE* g = (TGAE*) g_scetg_pth[iPtZ]->Clone ();
  
        g->SetFillColorAlpha (scetgColor, scetgAlpha);
        ((TGAE*) g->Clone ())->Draw ("3");
        ResetTGAEErrors (g);
        ResetXErrors (g);
        g->SetLineColor (scetgColor);
        g->SetLineWidth (4);
        ((TGAE*) g->Clone ())->Draw ("L");
        SaferDelete (&g);
      }

      {
        TGAE* g = (TGAE*) g_jewel_pth[iPtZ]->Clone ();
        SetMinErrors (g, minModelUnc, true);
  
        g->SetFillColorAlpha (jewelColor, jewelAlpha);
        ((TGAE*) g->Clone ())->Draw ("3");
        ResetTGAEErrors (g);
        ResetXErrors (g);
        g->SetLineColor (jewelColor);
        g->SetLineWidth (4);
        ((TGAE*) g->Clone ())->Draw ("L");
        SaferDelete (&g);
      }

      if (iPtZ > 2) {
        TGAE* g = (TGAE*) g_colbt_pth[iPtZ]->Clone ();
        RecenterGraph (g);
        ResetXErrors (g);

        g->SetFillColorAlpha (colbtColor, colbtAlpha);
        ((TGAE*) g->Clone ())->Draw ("3");
        ResetTGAEErrors (g);
        ResetXErrors (g);
        g->SetLineColor (colbtColor);
        g->SetLineWidth (4);
        ((TGAE*) g->Clone ())->Draw ("L");
        SaferDelete (&g);
      }

      {
        TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

        RecenterGraph (g_syst);
        ResetXErrors (g_syst);
        SetConstantXErrors (g_syst, 0.060, true, xmin, xmax);

        g_syst->SetMarkerSize (0);
        g_syst->SetLineWidth (1);
        g_syst->SetLineColor (finalColors[iPtZ-1]);
        g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

        g_syst->Draw ("5P");
      }

      {
        TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

        RecenterGraph (g_stat);
        ResetXErrors (g_stat);
        //deltaize (g_stat, 0.95, true);
        //ResetXErrors (g_stat);

        Style_t markerStyle = markerStyles[iPtZ-2];
        float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        g_stat->SetLineWidth (3);
        g_stat->SetMarkerColor (finalColors[iPtZ-1]);
        g_stat->SetLineColor (finalColors[iPtZ-1]);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        markerStyle = FullToOpenMarker (markerStyle);
        markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
        
        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        g_stat->SetLineWidth (0);
        g_stat->SetMarkerColor (kBlack);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        SaferDelete (&g_stat);
      }

      if (iPtZ == 2)      MakeDataBox   (0.46, 0.785, finalFillColors[iPtZ-1], 0.30, markerStyles[iPtZ-2], (markerStyles[iPtZ-2] == kFullDiamond || markerStyles[iPtZ-2] == kOpenDiamond ? 2.4 : 1.8), 1.2, 1.2);
      else if (iPtZ == 3) MakeDataBox   (0.44, 0.785, finalFillColors[iPtZ-1], 0.30, markerStyles[iPtZ-2], (markerStyles[iPtZ-2] == kFullDiamond || markerStyles[iPtZ-2] == kOpenDiamond ? 2.4 : 1.8), 1.2*xPadLCMiddle/(xPadCRMiddle-xPadLCMiddle), 1.2);
      else if (iPtZ == 4) MakeDataBox   (0.44, 0.885, finalFillColors[iPtZ-1], 0.30, markerStyles[iPtZ-2], (markerStyles[iPtZ-2] == kFullDiamond || markerStyles[iPtZ-2] == kOpenDiamond ? 2.4 : 1.8), 1.2*xPadLCMiddle/(1.-xPadCRMiddle), 1.2);

      tl->SetTextColor (kBlack);
      tl->SetTextSize (32);
      tl->SetTextAlign (11);
      if (iPtZ == 2)      tl->DrawLatexNDC (0.45-0.007, 0.785-0.013, "Data, 15 < #it{p}_{T}^{#it{Z}} < 30 GeV");
      else if (iPtZ == 3) tl->DrawLatexNDC (0.43-0.007, 0.785-0.013, "Data, 30 < #it{p}_{T}^{#it{Z}} < 60 GeV");
      else if (iPtZ == 4) tl->DrawLatexNDC (0.43-0.007, 0.885-0.013, "Data, #it{p}_{T}^{#it{Z}} > 60 GeV");
    } // end loop over iPtZ

    lPad->cd ();
    tl->SetTextSize (36);
    tl->DrawLatexNDC (0.30, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (34);
    tl->DrawLatexNDC (0.34, 0.835, "0-10\% Pb+Pb #/#it{pp}");

    cPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.18, 0.890, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.18, 0.840, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    rPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.43-0.007, 0.828-0.013, "Hybrid Model");
    tl->DrawLatexNDC (0.43-0.007, 0.773-0.013, "CoLBT-hydro");
    tl->DrawLatexNDC (0.43-0.007, 0.718-0.013, "SCET_{G} (#it{g}^{ }=^{ }2.0#pm0.2)");
    tl->DrawLatexNDC (0.43-0.007, 0.663-0.013, "JEWEL");

    MakeTheoryBox (0.44, 0.828, hybridColor, hybridAlpha, 1.2*xPadLCMiddle/(1.-xPadCRMiddle), 1.2);
    MakeTheoryLine (0.44, 0.828, hybridColor, 1.2*xPadLCMiddle/(1.-xPadCRMiddle));
    MakeTheoryBox (0.44, 0.773, colbtColor, colbtAlpha, 1.2*xPadLCMiddle/(1.-xPadCRMiddle), 1.2);
    MakeTheoryLine (0.44, 0.773, colbtColor, 1.2*xPadLCMiddle/(1.-xPadCRMiddle));
    MakeTheoryBox (0.44, 0.718, scetgColor, scetgAlpha, 1.2*xPadLCMiddle/(1.-xPadCRMiddle), 1.2);
    MakeTheoryLine (0.44, 0.718, scetgColor, 1.2*xPadLCMiddle/(1.-xPadCRMiddle));
    MakeTheoryBox (0.44, 0.663, jewelColor, jewelAlpha, 1.2*xPadLCMiddle/(1.-xPadCRMiddle), 1.2);
    MakeTheoryLine (0.44, 0.663, jewelColor, 1.2*xPadLCMiddle/(1.-xPadCRMiddle));

    c16->SaveAs (Form ("%s/FinalPlots/iaa_pTch_theoryComp.pdf", plotPath.Data ()));
  }




  //{
  //  TCanvas* c16 = new TCanvas ("c16", "", 800, 800);
  //  const double lMargin = 0.12;
  //  const double rMargin = 0.04;
  //  const double bMargin = 0.15;
  //  const double tMargin = 0.04;

  //  c16->SetLeftMargin (lMargin);
  //  c16->SetRightMargin (rMargin);
  //  c16->SetBottomMargin (bMargin);
  //  c16->SetTopMargin (tMargin);

  //  c16->SetLogx ();
  //  c16->SetLogy ();

  //  const short iPtZ = 2;
  //  const short iCent = 3;

  //  const double xmin = pTchBins[iPtZ][0];
  //  const double xmax = pTchBins[iPtZ][nPtchBins[iPtZ]];

  //  {
  //    TH1D* h = new TH1D ("", "", nPtchBins[iPtZ], pTchBins[iPtZ]);

  //    TAxis* xax = h->GetXaxis ();
  //    TAxis* yax = h->GetYaxis ();

  //    xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
  //    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
  //    xax->SetRangeUser (xmin, xmax);
  //    xax->SetLabelSize (0);

  //    yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
  //    yax->SetTitleOffset (0.8 * yax->GetTitleOffset ());
  //    yax->SetLabelSize (0);
  //    const double ymin = 0.18;
  //    const double ymax = 7;
  //    yax->SetRangeUser (ymin, ymax);

  //    h->SetLineWidth (0);

  //    h->DrawCopy ("");
  //    SaferDelete (&h);

  //    tl->SetTextFont (43);
  //    tl->SetTextSize (36);
  //    tl->SetTextAlign (21);

  //    const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
  //    //const double yoff = ymin - 0.05 * (ymax-ymin) / (1.-tMargin-bMargin);
  //    tl->DrawLatex (1,  yoff, "1");
  //    tl->DrawLatex (2,  yoff, "2");
  //    tl->DrawLatex (3,  yoff, "3");
  //    tl->DrawLatex (4,  yoff, "4");
  //    tl->DrawLatex (5,  yoff, "5");
  //    tl->DrawLatex (6,  yoff, "6");
  //    tl->DrawLatex (7,  yoff, "7");
  //    tl->DrawLatex (10, yoff, "10");
  //    tl->DrawLatex (15, yoff, "15");

  //    tl->SetTextAlign (32);

  //    const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-lMargin-rMargin));
  //    tl->DrawLatex (xoff, 0.2, "0.2");
  //    tl->DrawLatex (xoff, 0.3, "0.3");
  //    tl->DrawLatex (xoff, 0.4, "0.4");
  //    tl->DrawLatex (xoff, 0.5, "0.5");
  //    tl->DrawLatex (xoff, 0.7, "0.7");
  //    tl->DrawLatex (xoff, 1, "1");
  //    tl->DrawLatex (xoff, 2, "2");
  //    tl->DrawLatex (xoff, 3, "3");
  //    tl->DrawLatex (xoff, 4, "4");
  //    tl->DrawLatex (xoff, 5, "5");

  //    l->SetLineStyle (2);
  //    l->SetLineWidth (2);
  //    l->SetLineColor (kBlack);
  //    l->DrawLine (xmin, 1, xmax, 1);
  //  }

  //  {
  //    TGAE* g = (TGAE*) g_hybrid_pth[iPtZ]->Clone ();
  //    RecenterGraph (g);
  //    SetMinErrors (g, minModelUnc, true);

  //    g->SetFillColorAlpha (hybridColor, hybridAlpha);
  //    ((TGAE*) g->Clone ())->Draw ("3");
  //    ResetTGAEErrors (g);
  //    ResetXErrors (g);
  //    g->SetLineColor (hybridColor);
  //    g->SetLineWidth (4);
  //    ((TGAE*) g->Clone ())->Draw ("L");
  //    SaferDelete (&g);
  //  }

  //  {
  //    TGAE* g = (TGAE*) g_scetg_pth[iPtZ]->Clone ();
  //
  //    g->SetFillColorAlpha (scetgColor, scetgAlpha);
  //    ((TGAE*) g->Clone ())->Draw ("3");
  //    ResetTGAEErrors (g);
  //    ResetXErrors (g);
  //    g->SetLineColor (scetgColor);
  //    g->SetLineWidth (4);
  //    ((TGAE*) g->Clone ())->Draw ("L");
  //    SaferDelete (&g);
  //  }

  //  {
  //    TGAE* g = (TGAE*) g_jewel_pth[iPtZ]->Clone ();
  //    SetMinErrors (g, minModelUnc, true);
  //
  //    g->SetFillColorAlpha (jewelColor, jewelAlpha);
  //    ((TGAE*) g->Clone ())->Draw ("3");
  //    ResetTGAEErrors (g);
  //    ResetXErrors (g);
  //    g->SetLineColor (jewelColor);
  //    g->SetLineWidth (4);
  //    ((TGAE*) g->Clone ())->Draw ("L");
  //    SaferDelete (&g);
  //  }

  //  {
  //    TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

  //    RecenterGraph (g_syst);
  //    ResetXErrors (g_syst);
  //    SetConstantXErrors (g_syst, 0.060, true, xmin, xmax);

  //    g_syst->SetMarkerSize (0);
  //    g_syst->SetLineWidth (1);
  //    g_syst->SetLineColor (finalColors[iPtZ-1]);
  //    g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

  //    g_syst->Draw ("5P");
  //  }

  //  {
  //    TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

  //    RecenterGraph (g_stat);
  //    ResetXErrors (g_stat);
  //    //deltaize (g_stat, 0.95, true);
  //    //ResetXErrors (g_stat);

  //    Style_t markerStyle = markerStyles[iPtZ-2];
  //    float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
  //    g_stat->SetMarkerStyle (markerStyle);
  //    g_stat->SetMarkerSize (markerSize);
  //    g_stat->SetLineWidth (3);
  //    g_stat->SetMarkerColor (finalColors[iPtZ-1]);
  //    g_stat->SetLineColor (finalColors[iPtZ-1]);

  //    ((TGAE*) g_stat->Clone ())->Draw ("P");

  //    markerStyle = FullToOpenMarker (markerStyle);
  //    markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
  //    
  //    g_stat->SetMarkerStyle (markerStyle);
  //    g_stat->SetMarkerSize (markerSize);
  //    g_stat->SetLineWidth (0);
  //    g_stat->SetMarkerColor (kBlack);

  //    ((TGAE*) g_stat->Clone ())->Draw ("P");

  //    SaferDelete (&g_stat);
  //  }
  //  //{
  //  //  TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat_2015proj[iPtZ][iCent]);

  //  //  RecenterGraph (g_stat);
  //  //  ResetXErrors (g_stat);
  //  //  deltaize (g_stat, 1.05, true);
  //  //  ResetXErrors (g_stat);

  //  //  g_stat->SetMarkerStyle (kFullCircle);
  //  //  g_stat->SetMarkerSize (2.3);
  //  //  g_stat->SetLineWidth (3);
  //  //  //g_stat->SetMarkerColor (kBlack);
  //  //  //g_stat->SetLineColor (kBlack);
  //  //  g_stat->SetMarkerColor (finalColors[2]);
  //  //  g_stat->SetLineColor (finalColors[2]);

  //  //  ((TGAE*) g_stat->Clone ())->Draw ("P");

  //  //  g_stat->SetMarkerStyle (kOpenCircle);
  //  //  g_stat->SetMarkerSize (2.3);
  //  //  g_stat->SetLineWidth (0);
  //  //  g_stat->SetMarkerColor (kBlack);

  //  //  ((TGAE*) g_stat->Clone ())->Draw ("P");

  //  //  SaferDelete (&g_stat);
  //  //}

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextAlign (11);

  //  tl->SetTextSize (32);
  //  tl->DrawLatexNDC (0.30, 0.890, "#bf{#it{ATLAS}} Internal");
  //  tl->SetTextSize (30);
  //  tl->DrawLatexNDC (0.33, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
  //  tl->DrawLatexNDC (0.33, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

  //  tl->SetTextSize (25);
  //  tl->DrawLatexNDC (0.43-0.009, 0.760-0.010, "ATLAS 0-10\%, 15 < #it{p}_{T}^{#it{Z}} < 30 GeV");
  //  tl->DrawLatexNDC (0.43-0.009, 0.715-0.010, "Hybrid Model");
  //  tl->DrawLatexNDC (0.43-0.009, 0.670-0.010, "SCET_{G} (#it{g}^{ }=^{ }2.0#pm0.2)");
  //  tl->DrawLatexNDC (0.78-0.009, 0.670-0.010, "JEWEL");

  //  MakeDataBox   (0.44, 0.760, finalFillColors[1], 0.30, markerStyles[0], 1.8);
  //  MakeTheoryBox (0.44, 0.715, hybridColor, hybridAlpha);
  //  MakeTheoryLine (0.44, 0.715, hybridColor);
  //  MakeTheoryBox (0.44, 0.670, scetgColor, scetgAlpha);
  //  MakeTheoryLine (0.44, 0.670, scetgColor);
  //  MakeTheoryBox (0.79, 0.670, jewelColor, jewelAlpha);
  //  MakeTheoryLine (0.79, 0.670, jewelColor);

  //  c16->SaveAs (Form ("%s/FinalPlots/iaa_pTch_theoryComp_iPtZ2.pdf", plotPath.Data ()));
  //}




  {
    TCanvas* c17 = new TCanvas ("c17", "", 800, 800);
    const double lMargin = 0.12;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c17->SetLeftMargin (lMargin);
    c17->SetRightMargin (rMargin);
    c17->SetBottomMargin (bMargin);
    c17->SetTopMargin (tMargin);

    c17->SetLogx ();
    c17->SetLogy ();

    const short iPtZ = 2;
    const short iCent = 3;

    const double xmin = pTchBins[iPtZ][0];
    const double xmax = pTchBins[iPtZ][nPtchBins[iPtZ]];

    {
      TH1D* h = new TH1D ("", "", nPtchBins[iPtZ], pTchBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (xmin, xmax);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
      yax->SetTitleOffset (0.8 * yax->GetTitleOffset ());
      yax->SetLabelSize (0);
      const double ymin = 0.18;
      const double ymax = 7;
      yax->SetRangeUser (ymin, ymax);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);

      const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      //const double yoff = ymin - 0.05 * (ymax-ymin) / (1.-tMargin-bMargin);
      tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (3,  yoff, "3");
      tl->DrawLatex (4,  yoff, "4");
      tl->DrawLatex (5,  yoff, "5");
      tl->DrawLatex (6,  yoff, "6");
      tl->DrawLatex (7,  yoff, "7");
      tl->DrawLatex (10, yoff, "10");
      tl->DrawLatex (15, yoff, "15");

      tl->SetTextAlign (32);

      const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-lMargin-rMargin));
      tl->DrawLatex (xoff, 0.2, "0.2");
      tl->DrawLatex (xoff, 0.3, "0.3");
      tl->DrawLatex (xoff, 0.4, "0.4");
      tl->DrawLatex (xoff, 0.5, "0.5");
      tl->DrawLatex (xoff, 0.7, "0.7");
      tl->DrawLatex (xoff, 1, "1");
      tl->DrawLatex (xoff, 2, "2");
      tl->DrawLatex (xoff, 3, "3");
      tl->DrawLatex (xoff, 4, "4");
      tl->DrawLatex (xoff, 5, "5");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);
      l->DrawLine (xmin, 1, xmax, 1);
    }

    {
      TGAE* g = (TGAE*) g_hybrid_pth[iPtZ]->Clone ();
      RecenterGraph (g);
      SetMinErrors (g, minModelUnc, true);

      g->SetFillColorAlpha (hybridColor, hybridAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (hybridColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_scetg_pth[iPtZ]->Clone ();
  
      g->SetFillColorAlpha (scetgColor, scetgAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (scetgColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_jewel_pth[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);
  
      g->SetFillColorAlpha (jewelColor, jewelAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (jewelColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

      RecenterGraph (g_syst);
      ResetXErrors (g_syst);
      SetConstantXErrors (g_syst, 0.060, true, xmin, xmax);

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      g_syst->SetLineColor (finalColors[iPtZ-1]);
      g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      g_syst->Draw ("5P");
    }

    {
      TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      //deltaize (g_stat, 0.95, true);
      //ResetXErrors (g_stat);

      Style_t markerStyle = markerStyles[iPtZ-2];
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (3);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      markerStyle = FullToOpenMarker (markerStyle);
      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (0);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    }
    //{
    //  TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat_2015proj[iPtZ][iCent]);

    //  RecenterGraph (g_stat);
    //  ResetXErrors (g_stat);
    //  deltaize (g_stat, 1.05, true);
    //  ResetXErrors (g_stat);

    //  g_stat->SetMarkerStyle (kFullCircle);
    //  g_stat->SetMarkerSize (2.3);
    //  g_stat->SetLineWidth (3);
    //  //g_stat->SetMarkerColor (kBlack);
    //  //g_stat->SetLineColor (kBlack);
    //  g_stat->SetMarkerColor (finalColors[2]);
    //  g_stat->SetLineColor (finalColors[2]);

    //  ((TGAE*) g_stat->Clone ())->Draw ("P");

    //  g_stat->SetMarkerStyle (kOpenCircle);
    //  g_stat->SetMarkerSize (2.3);
    //  g_stat->SetLineWidth (0);
    //  g_stat->SetMarkerColor (kBlack);

    //  ((TGAE*) g_stat->Clone ())->Draw ("P");

    //  SaferDelete (&g_stat);
    //}

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.33, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.33, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (25);
    tl->DrawLatexNDC (0.43-0.009, 0.760-0.010, "ATLAS 0-10\%, 15 < #it{p}_{T}^{#it{Z}} < 30 GeV");
    tl->DrawLatexNDC (0.43-0.009, 0.715-0.010, "Hybrid Model");
    tl->DrawLatexNDC (0.43-0.009, 0.670-0.010, "SCET_{G} (#it{g}^{ }=^{ }2.0#pm0.2)");
    tl->DrawLatexNDC (0.78-0.009, 0.670-0.010, "JEWEL");

    MakeDataBox   (0.44, 0.760, finalFillColors[1], 0.30, markerStyles[0], 1.8);
    MakeTheoryBox (0.44, 0.715, hybridColor, hybridAlpha);
    MakeTheoryLine (0.44, 0.715, hybridColor);
    MakeTheoryBox (0.44, 0.670, scetgColor, scetgAlpha);
    MakeTheoryLine (0.44, 0.670, scetgColor);
    MakeTheoryBox (0.79, 0.670, jewelColor, jewelAlpha);
    MakeTheoryLine (0.79, 0.670, jewelColor);

    c17->SaveAs (Form ("%s/FinalPlots/iaa_pTch_theoryComp_iPtZ2.pdf", plotPath.Data ()));
  }




  {
    TCanvas* c18 = new TCanvas ("c18", "", 800, 800);
    const double lMargin = 0.12;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c18->SetLeftMargin (lMargin);
    c18->SetRightMargin (rMargin);
    c18->SetBottomMargin (bMargin);
    c18->SetTopMargin (tMargin);

    c18->SetLogx ();
    c18->SetLogy ();

    const short iPtZ = 3;
    const short iCent = 3;

    const double xmin = pTchBins[iPtZ][0];
    const double xmax = pTchBins[iPtZ][nPtchBins[iPtZ]];

    {
      TH1D* h = new TH1D ("", "", nPtchBins[iPtZ], pTchBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (xmin, xmax);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
      yax->SetTitleOffset (0.8 * yax->GetTitleOffset ());
      yax->SetLabelSize (0);
      const double ymin = 0.18;
      const double ymax = 7;
      yax->SetRangeUser (ymin, ymax);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);

      const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      //const double yoff = ymin - 0.05 * (ymax-ymin) / (1.-tMargin-bMargin);
      tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (3,  yoff, "3");
      tl->DrawLatex (4,  yoff, "4");
      tl->DrawLatex (5,  yoff, "5");
      tl->DrawLatex (6,  yoff, "6");
      tl->DrawLatex (7,  yoff, "7");
      tl->DrawLatex (10, yoff, "10");
      tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");

      tl->SetTextAlign (32);

      const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-lMargin-rMargin));
      tl->DrawLatex (xoff, 0.2, "0.2");
      tl->DrawLatex (xoff, 0.3, "0.3");
      tl->DrawLatex (xoff, 0.4, "0.4");
      tl->DrawLatex (xoff, 0.5, "0.5");
      tl->DrawLatex (xoff, 0.7, "0.7");
      tl->DrawLatex (xoff, 1, "1");
      tl->DrawLatex (xoff, 2, "2");
      tl->DrawLatex (xoff, 3, "3");
      tl->DrawLatex (xoff, 4, "4");
      tl->DrawLatex (xoff, 5, "5");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);
      l->DrawLine (xmin, 1, xmax, 1);
    }

    {
      TGAE* g = (TGAE*) g_hybrid_pth[iPtZ]->Clone ();
      RecenterGraph (g);
      SetMinErrors (g, minModelUnc, true);

      g->SetFillColorAlpha (hybridColor, hybridAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (hybridColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_scetg_pth[iPtZ]->Clone ();
  
      g->SetFillColorAlpha (scetgColor, scetgAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (scetgColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_jewel_pth[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);
  
      g->SetFillColorAlpha (jewelColor, jewelAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (jewelColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_colbt_pth[iPtZ]->Clone ();
      RecenterGraph (g);
      ResetXErrors (g);

      g->SetFillColorAlpha (colbtColor, colbtAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (colbtColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

      RecenterGraph (g_syst);
      ResetXErrors (g_syst);
      SetConstantXErrors (g_syst, 0.060, true, xmin, xmax);

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      g_syst->SetLineColor (finalColors[iPtZ-1]);
      g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      g_syst->Draw ("5P");
    }

    {
      TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      //deltaize (g_stat, 0.95, true);
      //ResetXErrors (g_stat);

      Style_t markerStyle = markerStyles[iPtZ-2];
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (3);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      markerStyle = FullToOpenMarker (markerStyle);
      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (0);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    }
    //{
    //  TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat_2015proj[iPtZ][iCent]);

    //  RecenterGraph (g_stat);
    //  ResetXErrors (g_stat);
    //  deltaize (g_stat, 1.05, true);
    //  ResetXErrors (g_stat);

    //  g_stat->SetMarkerStyle (kFullCircle);
    //  g_stat->SetMarkerSize (2.3);
    //  g_stat->SetLineWidth (3);
    //  //g_stat->SetMarkerColor (kBlack);
    //  //g_stat->SetLineColor (kBlack);
    //  g_stat->SetMarkerColor (finalColors[2]);
    //  g_stat->SetLineColor (finalColors[2]);

    //  ((TGAE*) g_stat->Clone ())->Draw ("P");

    //  g_stat->SetMarkerStyle (kOpenCircle);
    //  g_stat->SetMarkerSize (2.3);
    //  g_stat->SetLineWidth (0);
    //  g_stat->SetMarkerColor (kBlack);

    //  ((TGAE*) g_stat->Clone ())->Draw ("P");

    //  SaferDelete (&g_stat);
    //}

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.33, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.33, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (25);
    tl->DrawLatexNDC (0.43-0.009, 0.760-0.010, "ATLAS 0-10\%, 30 < #it{p}_{T}^{#it{Z}} < 60 GeV");
    tl->DrawLatexNDC (0.43-0.009, 0.715-0.010, "Hybrid Model");
    tl->DrawLatexNDC (0.78-0.009, 0.715-0.010, "CoLBT-hydro");
    tl->DrawLatexNDC (0.43-0.009, 0.670-0.010, "SCET_{G} (#it{g}^{ }=^{ }2.0#pm0.2)");
    tl->DrawLatexNDC (0.78-0.009, 0.670-0.010, "JEWEL");

    MakeDataBox   (0.44, 0.760, finalFillColors[2], 0.30, markerStyles[1], 1.8);
    MakeTheoryBox (0.44, 0.715, hybridColor, hybridAlpha);
    MakeTheoryLine (0.44, 0.715, hybridColor);
    MakeTheoryBox (0.79, 0.715, colbtColor, colbtAlpha);
    MakeTheoryLine (0.79, 0.715, colbtColor);
    MakeTheoryBox (0.44, 0.670, scetgColor, scetgAlpha);
    MakeTheoryLine (0.44, 0.670, scetgColor);
    MakeTheoryBox (0.79, 0.670, jewelColor, jewelAlpha);
    MakeTheoryLine (0.79, 0.670, jewelColor);

    c18->SaveAs (Form ("%s/FinalPlots/iaa_pTch_theoryComp_iPtZ3.pdf", plotPath.Data ()));
  }




  {
    TCanvas* c19 = new TCanvas ("c19", "", 800, 800);
    const double lMargin = 0.12;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c19->SetLeftMargin (lMargin);
    c19->SetRightMargin (rMargin);
    c19->SetBottomMargin (bMargin);
    c19->SetTopMargin (tMargin);

    c19->SetLogx ();
    c19->SetLogy ();

    const short iPtZ = 4;
    const short iCent = 3;

    const double xmin = pTchBins[iPtZ][0];
    const double xmax = pTchBins[iPtZ][nPtchBins[iPtZ]];

    {
      TH1D* h = new TH1D ("", "", 1, xmin, xmax);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (xmin, xmax);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
      yax->SetTitleOffset (0.8 * yax->GetTitleOffset ());
      yax->SetLabelSize (0);
      const double ymin = 0.18;
      const double ymax = 7;
      yax->SetRangeUser (ymin, ymax);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);

      const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (3,  yoff, "3");
      tl->DrawLatex (4,  yoff, "4");
      tl->DrawLatex (5,  yoff, "5");
      tl->DrawLatex (6,  yoff, "6");
      tl->DrawLatex (7,  yoff, "7");
      tl->DrawLatex (10, yoff, "10");
      tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");

      tl->SetTextAlign (32);

      const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-lMargin-rMargin));
      tl->DrawLatex (xoff, 0.2, "0.2");
      tl->DrawLatex (xoff, 0.3, "0.3");
      tl->DrawLatex (xoff, 0.4, "0.4");
      tl->DrawLatex (xoff, 0.5, "0.5");
      tl->DrawLatex (xoff, 0.7, "0.7");
      tl->DrawLatex (xoff, 1, "1");
      tl->DrawLatex (xoff, 2, "2");
      tl->DrawLatex (xoff, 3, "3");
      tl->DrawLatex (xoff, 4, "4");
      tl->DrawLatex (xoff, 5, "5");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);
      l->DrawLine (xmin, 1, xmax, 1);
    }

    {
      TGAE* g = (TGAE*) g_hybrid_pth[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);

      g->SetFillColorAlpha (hybridColor, hybridAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (hybridColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_scetg_pth[iPtZ]->Clone ();
  
      g->SetFillColorAlpha (scetgColor, scetgAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (scetgColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_jewel_pth[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);
  
      g->SetFillColorAlpha (jewelColor, jewelAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (jewelColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_colbt_pth[iPtZ]->Clone ();
      RecenterGraph (g);
      ResetXErrors (g);

      g->SetFillColorAlpha (colbtColor, colbtAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (colbtColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();

      RecenterGraph (g_syst);
      ResetXErrors (g_syst);
      SetConstantXErrors (g_syst, 0.060, true, xmin, xmax);

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      g_syst->SetLineColor (finalColors[iPtZ-1]);
      g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      g_syst->Draw ("5P");
    }

    {
      TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      //deltaize (g_stat, 0.95, true);
      //ResetXErrors (g_stat);

      Style_t markerStyle = markerStyles[iPtZ-2];
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (3);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      markerStyle = FullToOpenMarker (markerStyle);
      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (0);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    }
    //{
    //  TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat_2015proj[iPtZ][iCent]);

    //  RecenterGraph (g_stat);
    //  ResetXErrors (g_stat);
    //  deltaize (g_stat, 1.05, true);
    //  ResetXErrors (g_stat);

    //  g_stat->SetMarkerStyle (kFullCircle);
    //  g_stat->SetMarkerSize (2.3);
    //  g_stat->SetLineWidth (3);
    //  //g_stat->SetMarkerColor (kBlack);
    //  //g_stat->SetLineColor (kBlack);
    //  g_stat->SetMarkerColor (finalColors[2]);
    //  g_stat->SetLineColor (finalColors[2]);

    //  ((TGAE*) g_stat->Clone ())->Draw ("P");

    //  g_stat->SetMarkerStyle (kOpenCircle);
    //  g_stat->SetMarkerSize (2.3);
    //  g_stat->SetLineWidth (0);
    //  g_stat->SetMarkerColor (kBlack);

    //  ((TGAE*) g_stat->Clone ())->Draw ("P");

    //  SaferDelete (&g_stat);
    //}

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.33, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.33, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (25);
    tl->DrawLatexNDC (0.43-0.007, 0.760-0.010, "ATLAS 0-10\%, #it{p}_{T}^{#it{Z}} > 60 GeV");
    tl->DrawLatexNDC (0.43-0.007, 0.715-0.010, "Hybrid Model");
    tl->DrawLatexNDC (0.78-0.007, 0.715-0.010, "CoLBT-hydro");
    tl->DrawLatexNDC (0.43-0.007, 0.670-0.010, "SCET_{G} (#it{g}^{ }=^{ }2.0#pm0.2)");
    tl->DrawLatexNDC (0.78-0.007, 0.670-0.010, "JEWEL");

    MakeDataBox   (0.44, 0.760, finalFillColors[3], 0.30, markerStyles[2], 2.4);
    MakeTheoryBox (0.44, 0.715, hybridColor, hybridAlpha);
    MakeTheoryLine (0.44, 0.715, hybridColor);
    MakeTheoryBox (0.79, 0.715, colbtColor, colbtAlpha);
    MakeTheoryLine (0.79, 0.715, colbtColor);
    MakeTheoryBox (0.44, 0.670, scetgColor, scetgAlpha);
    MakeTheoryLine (0.44, 0.670, scetgColor);
    MakeTheoryBox (0.79, 0.670, jewelColor, jewelAlpha);
    MakeTheoryLine (0.79, 0.670, jewelColor);

    c19->SaveAs (Form ("%s/FinalPlots/iaa_pTch_theoryComp_iPtZ4.pdf", plotPath.Data ()));
  }




  {
    TCanvas* c20 = new TCanvas ("c20", "", 800, 800);
    const double lMargin = 0.12;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c20->SetLeftMargin (lMargin);
    c20->SetRightMargin (rMargin);
    c20->SetBottomMargin (bMargin);
    c20->SetTopMargin (tMargin);

    c20->SetLogx ();
    c20->SetLogy ();

    const short iPtZ = 2;
    const short iCent = 3;

    const double xmin = (iPtZ == 2 ? 1./8. : xhZBins[iPtZ][0]);
    const double xmax = xhZBins[iPtZ][nXhZBins[iPtZ]];

    {
      TH1D* h = new TH1D ("", "", 1, xmin, xmax);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{x}_{h#it{Z}} = #it{p}_{T}^{ch} #/#it{p}_{T}^{#it{Z}}");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (xmin, xmax);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{x}_{h#it{Z}})");
      yax->SetTitleOffset (0.8 * yax->GetTitleOffset ());
      yax->SetLabelSize (0);
      const double ymin = 0.05;
      const double ymax = 7;
      yax->SetRangeUser (ymin, ymax);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);

      const double yoff = ymin / exp (0.054 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      tl->DrawLatex (2e-1,  yoff, "2#times10^{-1}");
      tl->DrawLatex (5e-1,  yoff, "5#times10^{-1}");
      tl->DrawLatex (1,     yoff, "1");

      tl->SetTextAlign (32);

      const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-lMargin-rMargin));
      tl->DrawLatex (xoff, 0.1, "0.1");
      tl->DrawLatex (xoff, 0.2, "0.2");
      tl->DrawLatex (xoff, 0.3, "0.3");
      tl->DrawLatex (xoff, 0.5, "0.5");
      tl->DrawLatex (xoff, 1, "1");
      tl->DrawLatex (xoff, 2, "2");
      tl->DrawLatex (xoff, 3, "3");
      tl->DrawLatex (xoff, 5, "5");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);
      l->DrawLine (xmin, 1, xmax, 1);
    }

    {
      TGAE* g = (TGAE*) g_hybrid_xhz[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);

      g->SetFillColorAlpha (hybridColor, hybridAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (hybridColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_scetg_xhz[iPtZ]->Clone ();
  
      g->SetFillColorAlpha (scetgColor, scetgAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (scetgColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_jewel_xhz[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);
  
      g->SetFillColorAlpha (jewelColor, jewelAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (jewelColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_iaa_syst[iPtZ][iCent]->Clone ();

      RecenterGraph (g_syst);
      ResetXErrors (g_syst);
      SetConstantXErrors (g_syst, 0.060, true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      //g_syst->SetLineColor (kBlack);
      //g_syst->SetFillColorAlpha (kGray, 0.3);
      g_syst->SetLineColor (finalColors[iPtZ-1]);
      g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      g_syst->Draw ("5P");
    }

    {
      TGAE* g_stat = make_graph (h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      //deltaize (g_stat, 0.95, true);
      //ResetXErrors (g_stat);

      Style_t markerStyle = markerStyles[iPtZ-2];
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (3);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      markerStyle = FullToOpenMarker (markerStyle);
      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (0);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.33, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.33, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (25);
    tl->DrawLatexNDC (0.32-0.007, 0.294-0.010, "ATLAS 0-10\%, 15 < #it{p}_{T}^{#it{Z}} < 30 GeV");
    tl->DrawLatexNDC (0.32-0.007, 0.247-0.010, "Hybrid Model");
    tl->DrawLatexNDC (0.32-0.007, 0.200-0.010, "SCET_{G} (#it{g}^{ }=^{ }2.0#pm0.2)");
    tl->DrawLatexNDC (0.69-0.007, 0.200-0.010, "JEWEL");

    MakeDataBox   (0.33, 0.294, finalFillColors[1], 0.30, markerStyles[0], 1.8);
    MakeTheoryBox (0.33, 0.247, hybridColor, hybridAlpha);
    MakeTheoryLine (0.33, 0.247, hybridColor);
    MakeTheoryBox (0.33, 0.200, scetgColor, scetgAlpha);
    MakeTheoryLine (0.33, 0.200, scetgColor);
    MakeTheoryBox (0.70, 0.200, jewelColor, jewelAlpha);
    MakeTheoryLine (0.70, 0.200, jewelColor);

    c20->SaveAs (Form ("%s/FinalPlots/iaa_xhZ_theoryComp_iPtZ2.pdf", plotPath.Data ()));
  }




  {
    TCanvas* c21 = new TCanvas ("c21", "", 800, 800);
    const double lMargin = 0.12;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c21->SetLeftMargin (lMargin);
    c21->SetRightMargin (rMargin);
    c21->SetBottomMargin (bMargin);
    c21->SetTopMargin (tMargin);

    c21->SetLogx ();
    c21->SetLogy ();

    const short iPtZ = 3;
    const short iCent = 3;

    const double xmin = xhZBins[iPtZ][0];
    const double xmax = xhZBins[iPtZ][nXhZBins[iPtZ]];

    {
      TH1D* h = new TH1D ("", "", 1, xmin, xmax);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{x}_{h#it{Z}} = #it{p}_{T}^{ch} #/#it{p}_{T}^{#it{Z}}");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (xmin, xmax);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{x}_{h#it{Z}})");
      yax->SetTitleOffset (0.8 * yax->GetTitleOffset ());
      yax->SetLabelSize (0);
      const double ymin = 0.05;
      const double ymax = 7;
      yax->SetRangeUser (ymin, ymax);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);

      const double yoff = ymin / exp (0.054 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
      tl->DrawLatex (1e-1,  yoff, "10^{-1}");
      tl->DrawLatex (2e-1,  yoff, "2#times10^{-1}");
      tl->DrawLatex (5e-1,  yoff, "5#times10^{-1}");
      tl->DrawLatex (1,     yoff, "1");

      tl->SetTextAlign (32);

      const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-lMargin-rMargin));
      tl->DrawLatex (xoff, 0.1, "0.1");
      tl->DrawLatex (xoff, 0.2, "0.2");
      tl->DrawLatex (xoff, 0.3, "0.3");
      tl->DrawLatex (xoff, 0.5, "0.5");
      tl->DrawLatex (xoff, 1, "1");
      tl->DrawLatex (xoff, 2, "2");
      tl->DrawLatex (xoff, 3, "3");
      tl->DrawLatex (xoff, 5, "5");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);
      l->DrawLine (xmin, 1, xmax, 1);
    }

    {
      TGAE* g = (TGAE*) g_hybrid_xhz[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);

      g->SetFillColorAlpha (hybridColor, hybridAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (hybridColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_scetg_xhz[iPtZ]->Clone ();
  
      g->SetFillColorAlpha (scetgColor, scetgAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (scetgColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_jewel_xhz[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);
  
      g->SetFillColorAlpha (jewelColor, jewelAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (jewelColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_colbt_xhz[iPtZ]->Clone ();
      RecenterGraph (g);
      ResetXErrors (g);

      g->SetFillColorAlpha (colbtColor, colbtAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (colbtColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_iaa_syst[iPtZ][iCent]->Clone ();

      RecenterGraph (g_syst);
      ResetXErrors (g_syst);
      SetConstantXErrors (g_syst, 0.060, true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      //g_syst->SetLineColor (kBlack);
      //g_syst->SetFillColorAlpha (kGray, 0.3);
      g_syst->SetLineColor (finalColors[iPtZ-1]);
      g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      g_syst->Draw ("5P");
    }

    {
      TGAE* g_stat = make_graph (h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      //deltaize (g_stat, 0.95, true);
      //ResetXErrors (g_stat);

      Style_t markerStyle = markerStyles[iPtZ-2];
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (3);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      markerStyle = FullToOpenMarker (markerStyle);
      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (0);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.33, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.33, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (25);
    tl->DrawLatexNDC (0.32-0.007, 0.294-0.010, "ATLAS 0-10\%, 30 < #it{p}_{T}^{#it{Z}} < 60 GeV");
    tl->DrawLatexNDC (0.32-0.007, 0.247-0.010, "Hybrid Model");
    tl->DrawLatexNDC (0.69-0.007, 0.247-0.010, "CoLBT-hydro");
    tl->DrawLatexNDC (0.32-0.007, 0.200-0.010, "SCET_{G} (#it{g}^{ }=^{ }2.0#pm0.2)");
    tl->DrawLatexNDC (0.69-0.007, 0.200-0.010, "JEWEL");

    MakeDataBox   (0.33, 0.294, finalFillColors[2], 0.30, markerStyles[1], 1.8);
    MakeTheoryBox (0.33, 0.247, hybridColor, hybridAlpha);
    MakeTheoryLine (0.33, 0.247, hybridColor);
    MakeTheoryBox (0.70, 0.247, colbtColor, colbtAlpha);
    MakeTheoryLine (0.70, 0.247, colbtColor);
    MakeTheoryBox (0.33, 0.200, scetgColor, scetgAlpha);
    MakeTheoryLine (0.33, 0.200, scetgColor);
    MakeTheoryBox (0.70, 0.200, jewelColor, jewelAlpha);
    MakeTheoryLine (0.70, 0.200, jewelColor);

    c21->SaveAs (Form ("%s/FinalPlots/iaa_xhZ_theoryComp_iPtZ3.pdf", plotPath.Data ()));
  }




  {
    TCanvas* c22 = new TCanvas ("c22", "", 800, 800);
    const double lMargin = 0.12;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c22->SetLeftMargin (lMargin);
    c22->SetRightMargin (rMargin);
    c22->SetBottomMargin (bMargin);
    c22->SetTopMargin (tMargin);

    c22->SetLogx ();
    c22->SetLogy ();

    const short iPtZ = 4;
    const short iCent = 3;

    const double xmin = xhZBins[iPtZ][0];
    const double xmax = xhZBins[iPtZ][nXhZBins[iPtZ]];

    {
      TH1D* h = new TH1D ("", "", 1, xmin, xmax);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{x}_{h#it{Z}} = #it{p}_{T}^{ch} #/#it{p}_{T}^{#it{Z}}");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (xmin, xmax);
      xax->SetLabelSize (0);

      yax->SetTitle ("#it{I}_{AA} (#it{x}_{h#it{Z}})");
      yax->SetTitleOffset (0.8 * yax->GetTitleOffset ());
      yax->SetLabelSize (0);
      const double ymin = 0.05;
      const double ymax = 7;
      yax->SetRangeUser (ymin, ymax);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);

      const double yoff = ymin / exp (0.054 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
      tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
      tl->DrawLatex (1e-1,  yoff, "10^{-1}");
      tl->DrawLatex (2e-1,  yoff, "2#times10^{-1}");
      tl->DrawLatex (5e-1,  yoff, "5#times10^{-1}");
      tl->DrawLatex (1,     yoff, "1");

      tl->SetTextAlign (32);

      const double xoff = xmin / exp (0.01 * (log (xmax) - log (xmin)) / (1.-lMargin-rMargin));
      tl->DrawLatex (xoff, 0.1, "0.1");
      tl->DrawLatex (xoff, 0.2, "0.2");
      tl->DrawLatex (xoff, 0.3, "0.3");
      tl->DrawLatex (xoff, 0.5, "0.5");
      tl->DrawLatex (xoff, 1, "1");
      tl->DrawLatex (xoff, 2, "2");
      tl->DrawLatex (xoff, 3, "3");
      tl->DrawLatex (xoff, 5, "5");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);
      l->DrawLine (xmin, 1, xmax, 1);
    }

    {
      TGAE* g = (TGAE*) g_hybrid_xhz[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);

      g->SetFillColorAlpha (hybridColor, hybridAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (hybridColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_scetg_xhz[iPtZ]->Clone ();
  
      g->SetFillColorAlpha (scetgColor, scetgAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (scetgColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_jewel_xhz[iPtZ]->Clone ();
      SetMinErrors (g, minModelUnc, true);
  
      g->SetFillColorAlpha (jewelColor, jewelAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (jewelColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g = (TGAE*) g_colbt_xhz[iPtZ]->Clone ();
      RecenterGraph (g);
      ResetXErrors (g);

      g->SetFillColorAlpha (colbtColor, colbtAlpha);
      ((TGAE*) g->Clone ())->Draw ("3");
      ResetTGAEErrors (g);
      ResetXErrors (g);
      g->SetLineColor (colbtColor);
      g->SetLineWidth (4);
      ((TGAE*) g->Clone ())->Draw ("L");
      SaferDelete (&g);
    }

    {
      TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_iaa_syst[iPtZ][iCent]->Clone ();

      RecenterGraph (g_syst);
      ResetXErrors (g_syst);
      SetConstantXErrors (g_syst, 0.060, true, xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);

      g_syst->SetMarkerSize (0);
      g_syst->SetLineWidth (1);
      g_syst->SetLineColor (finalColors[iPtZ-1]);
      g_syst->SetFillColorAlpha (finalFillColors[iPtZ-1], 0.3);

      g_syst->Draw ("5P");
    }

    {
      TGAE* g_stat = make_graph (h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]);

      RecenterGraph (g_stat);
      ResetXErrors (g_stat);
      //deltaize (g_stat, 0.95, true);
      //ResetXErrors (g_stat);

      Style_t markerStyle = markerStyles[iPtZ-2];
      float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (3);
      g_stat->SetMarkerColor (finalColors[iPtZ-1]);
      g_stat->SetLineColor (finalColors[iPtZ-1]);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      markerStyle = FullToOpenMarker (markerStyle);
      markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
      
      g_stat->SetMarkerStyle (markerStyle);
      g_stat->SetMarkerSize (markerSize);
      g_stat->SetLineWidth (0);
      g_stat->SetMarkerColor (kBlack);

      ((TGAE*) g_stat->Clone ())->Draw ("P");

      SaferDelete (&g_stat);
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.33, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.33, 0.800, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (25);
    tl->DrawLatexNDC (0.32-0.007, 0.294-0.010, "ATLAS 0-10\%, #it{p}_{T}^{#it{Z}} > 60 GeV");
    tl->DrawLatexNDC (0.32-0.007, 0.247-0.010, "Hybrid Model");
    tl->DrawLatexNDC (0.69-0.007, 0.247-0.010, "CoLBT-hydro");
    tl->DrawLatexNDC (0.32-0.007, 0.200-0.010, "SCET_{G} (#it{g}^{ }=^{ }2.0#pm0.2)");
    tl->DrawLatexNDC (0.69-0.007, 0.200-0.010, "JEWEL");

    MakeDataBox   (0.33, 0.294, finalFillColors[3], 0.30, markerStyles[2], 2.4);
    MakeTheoryBox (0.33, 0.247, hybridColor, hybridAlpha);
    MakeTheoryLine (0.33, 0.247, hybridColor);
    MakeTheoryBox (0.70, 0.247, colbtColor, colbtAlpha);
    MakeTheoryLine (0.70, 0.247, colbtColor);
    MakeTheoryBox (0.33, 0.200, scetgColor, scetgAlpha);
    MakeTheoryLine (0.33, 0.200, scetgColor);
    MakeTheoryBox (0.70, 0.200, jewelColor, jewelAlpha);
    MakeTheoryLine (0.70, 0.200, jewelColor);

    c22->SaveAs (Form ("%s/FinalPlots/iaa_xhZ_theoryComp_iPtZ4.pdf", plotPath.Data ()));
  }



  {
    TCanvas* c23 = new TCanvas ("c23", "", 800, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;
  
    c23->SetLeftMargin (lMargin);
    c23->SetRightMargin (rMargin);
    c23->SetBottomMargin (bMargin);
    c23->SetTopMargin (tMargin);
  
    c23->SetLogx ();
    c23->SetLogy ();

    short iPtZ = nPtZBins-1;
 
    const double xmin = pTchBins[iPtZ][0];
    const double xmax = pTchBins[iPtZ][nPtchBins[iPtZ]];

    {
      TH1D* h = new TH1D ("", "", 1, xmin, xmax);
  
      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();
  
      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetRangeUser (xmin, xmax);
      xax->SetLabelSize (0);
  
      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      yax->SetMoreLogLabels ();
      const double ymin = 0.18;
      const double ymax = 7;
      yax->SetRangeUser (ymin, ymax);
  
      h->SetLineWidth (0);
  
      h->DrawCopy ("");
      SaferDelete (&h);
  
      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);
  
      const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
      tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (3,  yoff, "3");
      tl->DrawLatex (4,  yoff, "4");
      tl->DrawLatex (5,  yoff, "5");
      tl->DrawLatex (6,  yoff, "6");
      tl->DrawLatex (7,  yoff, "7");
      tl->DrawLatex (10, yoff, "10");
      tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");
  
      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);
      l->DrawLine (xmin, 1, xmax, 1);
    }
  
    TGAE* g = nullptr;
  
    g = (TGAE*) tg_CMS_Zh_IAA_0_30_syst->Clone ();
    g->SetMarkerSize (0);
    g->SetMarkerColor (cmsColor);
    g->SetLineColor (kGray+2);
    g->SetLineWidth (1);
    g->SetFillColorAlpha (kGray+2, 0.3);
    ((TGE*) g->Clone ())->Draw ("5P");
  
    SaferDelete (&g);
  
    g = (TGAE*) tg_CMS_Zh_IAA_0_30_stat->Clone ();
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (3);
    ((TGE*) g->Clone ())->Draw ("P");

    SaferDelete (&g);

  
    TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][3]->Clone ();
    RecenterGraph (g_syst);
    ResetXErrors (g_syst);
    SetConstantXErrors (g_syst, 0.060, true, xmin, xmax);
  
    g_syst->SetMarkerSize (0);
    g_syst->SetMarkerColor (finalColors[3]);
    g_syst->SetLineColor (finalColors[3]);
    g_syst->SetLineWidth (1);
    g_syst->SetFillColorAlpha (finalFillColors[3], 0.3);
  
    ((TGAE*) g_syst->Clone ())->Draw ("5P");

    SaferDelete (&g_syst);
  
    TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][3]);
  
    RecenterGraph (g_stat);
    ResetXErrors (g_stat);
  
    g_stat->SetMarkerStyle (kFullDiamond);
    g_stat->SetMarkerSize (2.4);
    g_stat->SetLineWidth (3);
    g_stat->SetMarkerColor (finalColors[3]);
    g_stat->SetLineColor (finalColors[3]);
  
    ((TGAE*) g_stat->Clone ())->Draw ("P");
  
    g_stat->SetMarkerStyle (kOpenDiamond);
    g_stat->SetMarkerSize (2.4);
    g_stat->SetLineWidth (0);
    g_stat->SetMarkerColor (kBlack);
  
    ((TGAE*) g_stat->Clone ())->Draw ("P");
  
    SaferDelete (&g_stat);


    g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][2]->Clone ();
    RecenterGraph (g_syst);
    ResetXErrors (g_syst);
    SetConstantXErrors (g_syst, 0.060, true, xmin, xmax);
  
    g_syst->SetMarkerSize (0);
    g_syst->SetMarkerColor (finalColors[2]);
    g_syst->SetLineColor (finalColors[2]);
    g_syst->SetLineWidth (1);
    g_syst->SetFillColorAlpha (finalFillColors[2], 0.3);
  
    ((TGAE*) g_syst->Clone ())->Draw ("5P");

    SaferDelete (&g_syst);
  
    g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][2]);
  
    RecenterGraph (g_stat);
    ResetXErrors (g_stat);
  
    g_stat->SetMarkerStyle (kFullSquare);
    g_stat->SetMarkerSize (1.8);
    g_stat->SetLineWidth (3);
    g_stat->SetMarkerColor (finalColors[2]);
    g_stat->SetLineColor (finalColors[2]);
  
    ((TGAE*) g_stat->Clone ())->Draw ("P");
  
    g_stat->SetMarkerStyle (kOpenSquare);
    g_stat->SetMarkerSize (1.8);
    g_stat->SetLineWidth (0);
    g_stat->SetMarkerColor (kBlack);
  
    ((TGAE*) g_stat->Clone ())->Draw ("P");
  
    SaferDelete (&g_stat);
  
    myMarkerAndBoxAndLineText (0.33, 0.90-0.012, 1.4, 1001, finalFillColors[3], 0.30, finalColors[3], kFullDiamond,  2.4, "ATLAS #it{p}_{T}^{#it{Z}} > 60 GeV, 0-10%", 0.032);
    myMarkerAndBoxAndLineText (0.33, 0.85-0.012, 1.4, 1001, finalFillColors[2], 0.30, finalColors[2], kFullSquare,  1.8, "ATLAS #it{p}_{T}^{#it{Z}} > 60 GeV, 10-30%", 0.032);
    myMarkerAndBoxAndLineText (0.33, 0.80-0.012, 1.4, 1001, kGray+2, 0.30, kBlack,  kFullCircle, 1.8, "CMS #it{p}_{T}^{#it{Z}} > 30 GeV, 0-30% (#it{preliminary})", 0.032);
  
    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);
  
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.24, 0.310, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.26, 0.260, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.26, 0.210, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    //TPad* p_hamburglar = new TPad ("p_hamburglar", "", 0.65, 0.50, 0.90, 0.75);
    //p_hamburglar->Draw ();
    //p_hamburglar->cd ();

    //TImage* image = TImage::Open ("../DataThieving/CMS_Zh_IAA_HP2020/Hamburglar.jpg");
    //image->Draw ();


    c23->SaveAs (Form ("%s/FinalPlots/iaa_pTch_cmsComp.pdf", plotPath.Data ()));
  }




  {
    const char* canvasName = "c24";
    TCanvas* c24 = new TCanvas (canvasName, "", 1600, 1200);
    const double llMargin = 0.17;
    const double lrMargin = 0.032;
    const double clMargin = 0.032;
    const double crMargin = 0.032;
    const double rlMargin = 0.032;
    const double rrMargin = 0.040;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    const double deltaL = (1. - llMargin - lrMargin);
    const double deltaC = (1. - clMargin - crMargin);
    const double deltaR = (1. - rlMargin - rrMargin);

    const double a = (double) (deltaR * deltaC / (deltaL*deltaR + deltaC*deltaR + deltaL*deltaC));
    const double b = (double) (deltaR * deltaL / (deltaL*deltaR + deltaC*deltaR + deltaL*deltaC));

    const double xPadLCMiddle = a;
    const double xPadCRMiddle = a+b;

    const double yPadMiddle = 0.50;

    TPad* luPad = nullptr;
    TPad* cuPad = nullptr;
    TPad* ruPad = nullptr;
    TPad* ldPad = nullptr;
    TPad* cdPad = nullptr;
    TPad* rdPad = nullptr;

    luPad = new TPad (Form ("%s_luPad", canvasName), "", 0, yPadMiddle, xPadLCMiddle, 1);
    cuPad = new TPad (Form ("%s_cuPad", canvasName), "", xPadLCMiddle, yPadMiddle, xPadCRMiddle, 1);
    ruPad = new TPad (Form ("%s_ruPad", canvasName), "", xPadCRMiddle, yPadMiddle, 1, 1);
    ldPad = new TPad (Form ("%s_ldPad", canvasName), "", 0, 0, xPadLCMiddle, yPadMiddle);
    cdPad = new TPad (Form ("%s_cdPad", canvasName), "", xPadLCMiddle, 0, xPadCRMiddle, yPadMiddle);
    rdPad = new TPad (Form ("%s_rdPad", canvasName), "", xPadCRMiddle, 0, 1, yPadMiddle);

    luPad->SetLeftMargin (llMargin);
    luPad->SetRightMargin (lrMargin);
    cuPad->SetLeftMargin (clMargin);
    cuPad->SetRightMargin (crMargin);
    ruPad->SetLeftMargin (rlMargin);
    ruPad->SetRightMargin (rrMargin);
    luPad->SetBottomMargin (bMargin);
    luPad->SetTopMargin (tMargin);
    cuPad->SetBottomMargin (bMargin);
    cuPad->SetTopMargin (tMargin);
    ruPad->SetBottomMargin (bMargin);
    ruPad->SetTopMargin (tMargin);
    ldPad->SetLeftMargin (llMargin);
    ldPad->SetRightMargin (lrMargin);
    cdPad->SetLeftMargin (clMargin);
    cdPad->SetRightMargin (crMargin);
    rdPad->SetLeftMargin (rlMargin);
    rdPad->SetRightMargin (rrMargin);
    ldPad->SetBottomMargin (bMargin);
    ldPad->SetTopMargin (tMargin);
    cdPad->SetBottomMargin (bMargin);
    cdPad->SetTopMargin (tMargin);
    rdPad->SetBottomMargin (bMargin);
    rdPad->SetTopMargin (tMargin);

    TPad* pads[6] = {luPad, cuPad, ruPad, ldPad, cdPad, rdPad};
    for (int i = 0; i < 6; i++)
      pads[i]->Draw ();

    for (short iPtZ : {3, 4}) {
      const double xmin = pTchBins[iPtZ][0];
      const double xmax = pTchBins[iPtZ][nPtchBins[iPtZ]];

      for (short iCent : {1, 2, 3}) {
  
        pads[iCent-1 + 3*(iPtZ-3)]->cd ();
        {
          gPad->SetLogx ();
          gPad->SetLogy ();

          TH1D* h = new TH1D ("", "", 1, xmin, xmax);

          TAxis* xax = h->GetXaxis ();
          TAxis* yax = h->GetYaxis ();

          xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
          xax->SetRangeUser (xmin, xmax);
          xax->SetTitleFont (43);
          xax->SetTitleSize (30);
          xax->SetTitleOffset (2.5);
          xax->SetLabelSize (0);

          yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
          yax->SetMoreLogLabels ();
          const double ymin = 0.18;
          const double ymax = 7;
          yax->SetRangeUser (ymin, ymax);
          yax->SetTitleFont (43);
          yax->SetTitleSize (30);
          yax->SetTitleOffset (2.6);
          yax->SetLabelFont (43);
          if (iCent == 1) yax->SetLabelSize (28);
          else yax->SetLabelSize (0);

          h->SetLineWidth (0);

          h->DrawCopy ("");
          SaferDelete (&h);

          tl->SetTextFont (43);
          tl->SetTextSize (28);
          tl->SetTextAlign (21);

          const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
          tl->DrawLatex (1,  yoff, "1");
          tl->DrawLatex (2,  yoff, "2");
          tl->DrawLatex (3,  yoff, "3");
          tl->DrawLatex (4,  yoff, "4");
          tl->DrawLatex (5,  yoff, "5");
          tl->DrawLatex (6,  yoff, "6");
          tl->DrawLatex (7,  yoff, "7");
          tl->DrawLatex (10, yoff, "10");
          tl->DrawLatex (20, yoff, "20");
          tl->DrawLatex (30, yoff, "30");
          tl->DrawLatex (40, yoff, "40");
          tl->DrawLatex (60, yoff, "60");

          l->SetLineStyle (2);
          l->SetLineWidth (2);
          l->SetLineColor (kBlack);
          l->DrawLine (xmin, 1, xmax, 1);
        }


        // plot final result
        TGAE* g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst[iPtZ][iCent]->Clone ();
        RecenterGraph (g_syst);
        ResetXErrors (g_syst);
        deltaize (g_syst, 0.06, true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
        ResetXErrors (g_syst);
        SetConstantXErrors (g_syst, 0.060, true, pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
  
        g_syst->SetMarkerSize (0);
        g_syst->SetMarkerColor (finalColors[3]);
        g_syst->SetLineColor (finalColors[3]);
        g_syst->SetLineWidth (1);
        g_syst->SetFillColorAlpha (finalFillColors[3], 0.3);
  
        ((TGAE*) g_syst->Clone ())->Draw ("5P");

        SaferDelete (&g_syst);
  
        TGAE* g_stat = make_graph (h_trk_pt_ptz_iaa_stat[iPtZ][iCent]);
  
        RecenterGraph (g_stat);
        ResetXErrors (g_stat);
        deltaize (g_stat, 0.06, true, xmin, xmax);
        ResetXErrors (g_stat);
  
        g_stat->SetMarkerStyle (markerStyles[2]);
        g_stat->SetMarkerSize (2.4);
        g_stat->SetLineWidth (3);
        g_stat->SetMarkerColor (finalColors[3]);
        g_stat->SetLineColor (finalColors[3]);
  
        ((TGAE*) g_stat->Clone ())->Draw ("P");
  
        g_stat->SetMarkerStyle (FullToOpenMarker (markerStyles[2]));
        g_stat->SetMarkerSize (2.4);
        g_stat->SetLineWidth (0);
        g_stat->SetMarkerColor (kBlack);
  
        ((TGAE*) g_stat->Clone ())->Draw ("P");


        // plot CONF note result
        g_syst = (TGAE*) g_trk_pt_ptz_iaa_syst_conf[iPtZ][iCent]->Clone ();
        RecenterGraph (g_syst);
        ResetXErrors (g_syst);
        deltaize (g_syst, -0.06, true, xmin, xmax);
        ResetXErrors (g_syst);
        SetConstantXErrors (g_syst, 0.060, true, xmin, xmax);
  
        g_syst->SetMarkerSize (0);
        g_syst->SetMarkerColor (finalColors[1]);
        g_syst->SetLineColor (finalColors[1]);
        g_syst->SetLineWidth (1);
        g_syst->SetFillColorAlpha (finalFillColors[1], 0.3);
  
        ((TGAE*) g_syst->Clone ())->Draw ("5P");

        SaferDelete (&g_syst);
  
        g_stat = make_graph (h_trk_pt_ptz_iaa_stat_conf[iPtZ][iCent]);
  
        RecenterGraph (g_stat);
        ResetXErrors (g_stat);
        deltaize (g_stat, -0.06, true, xmin, xmax);
        ResetXErrors (g_stat);
  
        g_stat->SetMarkerStyle (markerStyles[0]);
        g_stat->SetMarkerSize (1.8);
        g_stat->SetLineWidth (3);
        g_stat->SetMarkerColor (finalColors[1]);
        g_stat->SetLineColor (finalColors[1]);
  
        ((TGAE*) g_stat->Clone ())->Draw ("P");
  
        g_stat->SetMarkerStyle (FullToOpenMarker (markerStyles[0]));
        g_stat->SetMarkerSize (1.8);
        g_stat->SetLineWidth (0);
        g_stat->SetMarkerColor (kBlack);
  
        ((TGAE*) g_stat->Clone ())->Draw ("P");
  
        SaferDelete (&g_stat);

      } // end loop over iCent
    } // end loop over iPtZ

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    luPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.21, 0.31, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.25, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.22, 0.19, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    ldPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.21, 0.31, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.25, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.22, 0.19, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (30);
    luPad->cd ();
    tl->DrawLatexNDC (0.50, 0.86, "30 < #it{p}_{T}^{#it{Z}} < 60 GeV");
    tl->DrawLatexNDC (0.50, 0.80, "30-80\%");
    cuPad->cd ();
    tl->DrawLatexNDC (0.50, 0.86, "30 < #it{p}_{T}^{#it{Z}} < 60 GeV");
    tl->DrawLatexNDC (0.50, 0.80, "10-30\%");
    ruPad->cd ();
    tl->DrawLatexNDC (0.50, 0.86, "30 < #it{p}_{T}^{#it{Z}} < 60 GeV");
    tl->DrawLatexNDC (0.50, 0.80, "0-10\%");

    ldPad->cd ();
    tl->DrawLatexNDC (0.50, 0.86, "#it{p}_{T}^{#it{Z}} > 60 GeV");
    tl->DrawLatexNDC (0.50, 0.80, "30-80\%");
    cdPad->cd ();
    tl->DrawLatexNDC (0.50, 0.86, "#it{p}_{T}^{#it{Z}} > 60 GeV");
    tl->DrawLatexNDC (0.50, 0.80, "10-30\%");
    rdPad->cd ();
    tl->DrawLatexNDC (0.50, 0.86, "#it{p}_{T}^{#it{Z}} > 60 GeV");
    tl->DrawLatexNDC (0.50, 0.80, "0-10\%");

    cuPad->cd ();
    myMarkerAndBoxAndLineText (0.22, 0.220, 3.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "Preliminary", 0.016 / (gPad->GetWNDC ()));
    ruPad->cd ();
    myMarkerAndBoxAndLineText (0.22, 0.220, 3.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "Final", 0.016 / (gPad->GetWNDC ()));

    cdPad->cd ();
    myMarkerAndBoxAndLineText (0.22, 0.220, 3.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "Preliminary", 0.016 / (gPad->GetWNDC ()));
    rdPad->cd ();
    myMarkerAndBoxAndLineText (0.22, 0.220, 3.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "Final", 0.016 / (gPad->GetWNDC ()));

    c24->SaveAs (Form ("%s/FinalPlots/iaa_pTch_CONFComp.pdf", plotPath.Data ()));
  }




  {
    const char* canvasName = "c25";
    TCanvas* c25 = new TCanvas (canvasName, "", 1600, 1200);
    const double llMargin = 0.17;
    const double lrMargin = 0.032;
    const double clMargin = 0.032;
    const double crMargin = 0.032;
    const double rlMargin = 0.032;
    const double rrMargin = 0.040;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    const double deltaL = (1. - llMargin - lrMargin);
    const double deltaC = (1. - clMargin - crMargin);
    const double deltaR = (1. - rlMargin - rrMargin);

    const double a = (double) (deltaR * deltaC / (deltaL*deltaR + deltaC*deltaR + deltaL*deltaC));
    const double b = (double) (deltaR * deltaL / (deltaL*deltaR + deltaC*deltaR + deltaL*deltaC));

    const double xPadLCMiddle = a;
    const double xPadCRMiddle = a+b;

    const double yPadMiddle = 0.5;

    TPad* luPad = nullptr;
    TPad* cuPad = nullptr;
    TPad* ruPad = nullptr;
    TPad* ldPad = nullptr;
    TPad* cdPad = nullptr;
    TPad* rdPad = nullptr;

    luPad = new TPad (Form ("%s_luPad", canvasName), "", 0, yPadMiddle, xPadLCMiddle, 1);
    cuPad = new TPad (Form ("%s_cuPad", canvasName), "", xPadLCMiddle, yPadMiddle, xPadCRMiddle, 1);
    ruPad = new TPad (Form ("%s_ruPad", canvasName), "", xPadCRMiddle, yPadMiddle, 1, 1);
    ldPad = new TPad (Form ("%s_ldPad", canvasName), "", 0, 0, xPadLCMiddle, yPadMiddle);
    cdPad = new TPad (Form ("%s_cdPad", canvasName), "", xPadLCMiddle, 0, xPadCRMiddle, yPadMiddle);
    rdPad = new TPad (Form ("%s_rdPad", canvasName), "", xPadCRMiddle, 0, 1, yPadMiddle);

    luPad->SetLeftMargin (llMargin);
    luPad->SetRightMargin (lrMargin);
    cuPad->SetLeftMargin (clMargin);
    cuPad->SetRightMargin (crMargin);
    ruPad->SetLeftMargin (rlMargin);
    ruPad->SetRightMargin (rrMargin);
    luPad->SetBottomMargin (bMargin);
    luPad->SetTopMargin (tMargin);
    cuPad->SetBottomMargin (bMargin);
    cuPad->SetTopMargin (tMargin);
    ruPad->SetBottomMargin (bMargin);
    ruPad->SetTopMargin (tMargin);
    ldPad->SetLeftMargin (llMargin);
    ldPad->SetRightMargin (lrMargin);
    cdPad->SetLeftMargin (clMargin);
    cdPad->SetRightMargin (crMargin);
    rdPad->SetLeftMargin (rlMargin);
    rdPad->SetRightMargin (rrMargin);
    ldPad->SetBottomMargin (bMargin);
    ldPad->SetTopMargin (tMargin);
    cdPad->SetBottomMargin (bMargin);
    cdPad->SetTopMargin (tMargin);
    rdPad->SetBottomMargin (bMargin);
    rdPad->SetTopMargin (tMargin);

    TPad* pads[6] = {luPad, cuPad, ruPad, ldPad, cdPad, rdPad};
    for (int i = 0; i < 6; i++)
      pads[i]->Draw ();

    for (short iPtZ : {3, 4}) {
      const double xmin = xhZBins[iPtZ][0];
      const double xmax = xhZBins[iPtZ][nXhZBins[iPtZ]];

      for (short iCent : {1, 2, 3}) {
  
        pads[iCent-1 + 3*(iPtZ-3)]->cd ();
        {
          gPad->SetLogx ();
          gPad->SetLogy ();

          TH1D* h = new TH1D ("", "", 1, xmin, xmax);

          TAxis* xax = h->GetXaxis ();
          TAxis* yax = h->GetYaxis ();

          xax->SetTitle ("#it{x}_{h#it{Z}}");
          xax->SetRangeUser (xmin, xmax);
          xax->SetTitleFont (43);
          xax->SetTitleSize (30);
          xax->SetTitleOffset (2.5);
          xax->SetLabelSize (0);

          yax->SetTitle ("#it{I}_{AA} (#it{x}_{h#it{Z}})");
          yax->SetMoreLogLabels ();
          const double ymin = 0.05;
          const double ymax = 7;
          yax->SetRangeUser (ymin, ymax);
          yax->SetTitleFont (43);
          yax->SetTitleSize (30);
          yax->SetTitleOffset (2.6);
          yax->SetLabelFont (43);
          if (iCent == 1) yax->SetLabelSize (28);
          else yax->SetLabelSize (0);

          h->SetLineWidth (0);

          h->DrawCopy ("");
          SaferDelete (&h);

          tl->SetTextFont (43);
          tl->SetTextSize (28);
          tl->SetTextAlign (21);

          const double yoff = ymin / exp (0.05 * (log (ymax) - log (ymin)) / (1.-tMargin-bMargin));
          if (iPtZ > 2) {
            if (iPtZ > 3) tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
            else tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
          }
          tl->DrawLatex (1e-1,  yoff, "10^{-1}");
          tl->DrawLatex (2e-1,  yoff, "2#times10^{-1}");
          tl->DrawLatex (5e-1,  yoff, "5#times10^{-1}");
          tl->DrawLatex (1,     yoff, "1");

          l->SetLineStyle (2);
          l->SetLineWidth (2);
          l->SetLineColor (kBlack);
          l->DrawLine (xmin, 1, xmax, 1);
        }


        // plot final result
        TGAE* g_syst = (TGAE*) g_trk_xhz_ptz_iaa_syst[iPtZ][iCent]->Clone ();
        RecenterGraph (g_syst);
        ResetXErrors (g_syst);
        deltaize (g_syst, 0.06, true, xmin, xmax);
        ResetXErrors (g_syst);
        SetConstantXErrors (g_syst, 0.060, true, xmin, xmax);
  
        g_syst->SetMarkerSize (0);
        g_syst->SetMarkerColor (finalColors[3]);
        g_syst->SetLineColor (finalColors[3]);
        g_syst->SetLineWidth (1);
        g_syst->SetFillColorAlpha (finalFillColors[3], 0.3);
  
        ((TGAE*) g_syst->Clone ())->Draw ("5P");

        SaferDelete (&g_syst);
  
        TGAE* g_stat = make_graph (h_trk_xhz_ptz_iaa_stat[iPtZ][iCent]);
  
        RecenterGraph (g_stat);
        ResetXErrors (g_stat);
        deltaize (g_stat, 0.06, true, xmin, xmax);
        ResetXErrors (g_stat);
  
        g_stat->SetMarkerStyle (markerStyles[2]);
        g_stat->SetMarkerSize (2.4);
        g_stat->SetLineWidth (3);
        g_stat->SetMarkerColor (finalColors[3]);
        g_stat->SetLineColor (finalColors[3]);
  
        ((TGAE*) g_stat->Clone ())->Draw ("P");
  
        g_stat->SetMarkerStyle (FullToOpenMarker (markerStyles[2]));
        g_stat->SetMarkerSize (2.4);
        g_stat->SetLineWidth (0);
        g_stat->SetMarkerColor (kBlack);
  
        ((TGAE*) g_stat->Clone ())->Draw ("P");


        // plot CONF note result
        g_syst = (TGAE*) g_trk_xhz_ptz_iaa_syst_conf[iPtZ][iCent]->Clone ();
        RecenterGraph (g_syst);
        ResetXErrors (g_syst);
        deltaize (g_syst, -0.06, true, xmin, xmax);
        ResetXErrors (g_syst);
        SetConstantXErrors (g_syst, 0.060, true, xmin, xmax);
  
        g_syst->SetMarkerSize (0);
        g_syst->SetMarkerColor (finalColors[1]);
        g_syst->SetLineColor (finalColors[1]);
        g_syst->SetLineWidth (1);
        g_syst->SetFillColorAlpha (finalFillColors[1], 0.3);
  
        ((TGAE*) g_syst->Clone ())->Draw ("5P");

        SaferDelete (&g_syst);
  
        g_stat = make_graph (h_trk_xhz_ptz_iaa_stat_conf[iPtZ][iCent]);
  
        RecenterGraph (g_stat);
        ResetXErrors (g_stat);
        deltaize (g_stat, -0.06, true, xmin, xmax);
        ResetXErrors (g_stat);
  
        g_stat->SetMarkerStyle (markerStyles[0]);
        g_stat->SetMarkerSize (1.8);
        g_stat->SetLineWidth (3);
        g_stat->SetMarkerColor (finalColors[1]);
        g_stat->SetLineColor (finalColors[1]);
  
        ((TGAE*) g_stat->Clone ())->Draw ("P");
  
        g_stat->SetMarkerStyle (FullToOpenMarker (markerStyles[0]));
        g_stat->SetMarkerSize (1.8);
        g_stat->SetLineWidth (0);
        g_stat->SetMarkerColor (kBlack);
  
        ((TGAE*) g_stat->Clone ())->Draw ("P");
  
        SaferDelete (&g_stat);

      } // end loop over iCent
    } // end loop over iPtZ

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    luPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.21, 0.31, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.25, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.22, 0.19, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    ldPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.21, 0.31, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.25, "#it{pp}, #sqrt{s} = 5.02 TeV, 260 pb^{-1}");
    tl->DrawLatexNDC (0.22, 0.19, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4-1.7 nb^{-1}");

    tl->SetTextSize (30);
    luPad->cd ();
    tl->DrawLatexNDC (0.50, 0.86, "30 < #it{p}_{T}^{#it{Z}} < 60 GeV");
    tl->DrawLatexNDC (0.50, 0.80, "30-80\%");
    cuPad->cd ();
    tl->DrawLatexNDC (0.50, 0.86, "30 < #it{p}_{T}^{#it{Z}} < 60 GeV");
    tl->DrawLatexNDC (0.50, 0.80, "10-30\%");
    ruPad->cd ();
    tl->DrawLatexNDC (0.50, 0.86, "30 < #it{p}_{T}^{#it{Z}} < 60 GeV");
    tl->DrawLatexNDC (0.50, 0.80, "0-10\%");

    ldPad->cd ();
    tl->DrawLatexNDC (0.50, 0.86, "#it{p}_{T}^{#it{Z}} > 60 GeV");
    tl->DrawLatexNDC (0.50, 0.80, "30-80\%");
    cdPad->cd ();
    tl->DrawLatexNDC (0.50, 0.86, "#it{p}_{T}^{#it{Z}} > 60 GeV");
    tl->DrawLatexNDC (0.50, 0.80, "10-30\%");
    rdPad->cd ();
    tl->DrawLatexNDC (0.50, 0.86, "#it{p}_{T}^{#it{Z}} > 60 GeV");
    tl->DrawLatexNDC (0.50, 0.80, "0-10\%");

    cuPad->cd ();
    myMarkerAndBoxAndLineText (0.22, 0.220, 3.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "Preliminary", 0.016 / (gPad->GetWNDC ()));
    ruPad->cd ();
    myMarkerAndBoxAndLineText (0.22, 0.220, 3.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "Final", 0.016 / (gPad->GetWNDC ()));

    cdPad->cd ();
    myMarkerAndBoxAndLineText (0.22, 0.220, 3.0, 1001, finalFillColors[1], 0.30, finalColors[1], markerStyles[0], 1.8, "Preliminary", 0.016 / (gPad->GetWNDC ()));
    rdPad->cd ();
    myMarkerAndBoxAndLineText (0.22, 0.220, 3.0, 1001, finalFillColors[3], 0.30, finalColors[3], markerStyles[2], 2.4, "Final", 0.016 / (gPad->GetWNDC ()));

    c25->SaveAs (Form ("%s/FinalPlots/iaa_xhZ_CONFComp.pdf", plotPath.Data ()));
  }




  {
    const char* canvasName = "c26";
    TCanvas* c26 = new TCanvas (canvasName, "", 1400, 1400);
    c26->cd ();

    const double lMargin = 0.15*4./7.;
    const double rMargin = 0.04*4./7.;
    const double bMargin = 0.22;
    const double tMargin = 0.13;

    const double deltaD = (1. - bMargin);
    const double deltaC = (1.);
    const double deltaU = (1. - tMargin);

    const double a = (double) (deltaU * deltaC / (deltaD*deltaU + deltaC*deltaU + deltaD*deltaC));
    const double b = (double) (deltaU * deltaD / (deltaD*deltaU + deltaC*deltaU + deltaD*deltaC));

    const double yPadDCMiddle = a;
    const double yPadCUMiddle = a+b;

    const double padX = 0.5*(1.+lMargin-rMargin);

    TPad* ulPad = new TPad (Form ("%s_ulPad", canvasName), "", 0, yPadCUMiddle, padX, 1);
    TPad* urPad = new TPad (Form ("%s_urPad", canvasName), "", padX, yPadCUMiddle, 1, 1);
    TPad* clPad = new TPad (Form ("%s_clPad", canvasName), "", 0, yPadDCMiddle, padX, yPadCUMiddle);
    TPad* crPad = new TPad (Form ("%s_crPad", canvasName), "", padX, yPadDCMiddle, 1, yPadCUMiddle);
    TPad* dlPad = new TPad (Form ("%s_dlPad", canvasName), "", 0, 0, padX, yPadDCMiddle);
    TPad* drPad = new TPad (Form ("%s_drPad", canvasName), "", padX, 0, 1, yPadDCMiddle);

    ulPad->SetLeftMargin (lMargin/padX);
    ulPad->SetRightMargin (0);
    urPad->SetLeftMargin (0);
    urPad->SetRightMargin (rMargin/(1.-padX));
    ulPad->SetBottomMargin (0);
    ulPad->SetTopMargin (tMargin);
    urPad->SetBottomMargin (0);
    urPad->SetTopMargin (tMargin);
    clPad->SetLeftMargin (lMargin/padX);
    clPad->SetRightMargin (0);
    crPad->SetLeftMargin (0);
    crPad->SetRightMargin (rMargin/(1.-padX));
    clPad->SetBottomMargin (0);
    clPad->SetTopMargin (0);
    crPad->SetBottomMargin (0);
    crPad->SetTopMargin (0);
    dlPad->SetLeftMargin (lMargin/padX);
    dlPad->SetRightMargin (0);
    drPad->SetLeftMargin (0);
    drPad->SetRightMargin (rMargin/(1.-padX));
    dlPad->SetBottomMargin (bMargin);
    dlPad->SetTopMargin (0);
    drPad->SetBottomMargin (bMargin);
    drPad->SetTopMargin (0);

    ulPad->Draw ();
    urPad->Draw ();
    clPad->Draw ();
    crPad->Draw ();
    dlPad->Draw ();
    drPad->Draw ();

    std::vector <TPad*> pads = {ulPad, urPad, clPad, crPad, dlPad, drPad};
    std::vector <TH1D***> stats = {h_trk_dphi_ptz_lt4_sub_stat, h_trk_dphi_ptz_gt4_sub_stat};
    std::vector <TGAE***> systs = {g_trk_dphi_ptz_lt4_sub_syst, g_trk_dphi_ptz_gt4_sub_syst};

    const double xmin = 0;
    const double xmax = pi;

    for (short iPad = 0; iPad < 6; iPad++) {
      const short iPtZ = (iPad < 2 ? 4 : (iPad < 4 ? 3 : 2));

      pads[iPad]->cd ();

      TH1D* h = new TH1D ("", "", 1, xmin, xmax);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetTitle ("#Delta#phi_{h#it{Z}}");
      xax->SetTitleFont (43);
      xax->SetTitleSize (40);
      xax->SetTitleOffset (2.0 * xax->GetTitleOffset ());
      xax->SetRangeUser (xmin, xmax);
      xax->SetLabelSize (0);

      if (pads[iPad] == ulPad) {
        yax->SetTitle ("(1/N_{Z}) (dN_{ch} / d#Delta#phi)");
        yax->SetTitleFont (43);
        yax->SetTitleSize (40);
        yax->SetTitleOffset (1.8 * yax->GetTitleOffset ());
      }
      else 
        yax->SetTitleSize (0);
      yax->SetLabelSize (0);
      const double ymin = (iPtZ == 4 ? -1.4 : (iPtZ == 3 ? -0.8 : -0.25));
      const double ymax = (iPtZ == 4 ? 11 : (iPtZ == 3 ? 3.8 : 1.35));
      yax->SetRangeUser (ymin, ymax);
      yax->SetNdivisions (805);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);

      const double yoff = ymin - 0.05 * (ymax-ymin) / (1.-tMargin-bMargin);
      if (pads[iPad] == ulPad || pads[iPad] == dlPad)
        tl->DrawLatex (0,  yoff, "0");
      tl->DrawLatex (0.5,  yoff, "0.5");
      tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (1.5,  yoff, "1.5");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (2.5,  yoff, "2.5");
      tl->DrawLatex (3,  yoff, "3");

      tl->SetTextAlign (32);
      const double xmin = 0;
      const double xmax = pi;
      const double xoff = xmin - 0.01 * (xmax - xmin) / (1.-lMargin);
      if (iPtZ == 2) {
        tl->DrawLatex (xoff, -0.2, "-0.2");
        tl->DrawLatex (xoff, 0, "0");
        tl->DrawLatex (xoff, 0.2, "0.2");
        tl->DrawLatex (xoff, 0.4, "0.4");
        tl->DrawLatex (xoff, 0.6, "0.6");
        tl->DrawLatex (xoff, 0.8, "0.8");
        tl->DrawLatex (xoff, 1.0, "1");
        tl->DrawLatex (xoff, 1.2, "1.2");
      }
      else if (iPtZ == 3) {
        tl->DrawLatex (xoff, 0, "0");
        tl->DrawLatex (xoff, 1, "1");
        tl->DrawLatex (xoff, 2, "2");
        tl->DrawLatex (xoff, 3, "3");
        tl->DrawLatex (xoff, 4, "4");
      }
      else {
        tl->DrawLatex (xoff, 0, "0");
        tl->DrawLatex (xoff, 2, "2");
        tl->DrawLatex (xoff, 4, "4");
        tl->DrawLatex (xoff, 6, "6");
        tl->DrawLatex (xoff, 8, "8");
        tl->DrawLatex (xoff, 10, "10");
      }

      l->SetLineStyle (7);
      l->SetLineWidth (1);
      l->SetLineColor (kBlack);
      l->DrawLine (3.*pi/4., ymin, 3.*pi/4., ymax);
      l->SetLineStyle (2);
      l->DrawLine (xmin, 0, xmax, 0);

      TBox* shadedBox = new TBox (3.*pi/4., ymin, pi, ymax);
      shadedBox->SetFillColorAlpha (kGray, 0.3);
      shadedBox->Draw ();
    } // end loop over iPad
    

    for (short iPad = 0; iPad < 6; iPad++) {
      const short iPtZ = (iPad < 2 ? 4 : (iPad < 4 ? 3 : 2));

      pads[iPad]->cd ();

      for (short iCent : {0, numCentBins}) {
        TGAE* g_syst = (TGAE*) systs[iPad%2][iPtZ][iCent]->Clone ();

        g_syst->SetMarkerSize (0);
        g_syst->SetLineWidth (1);
        g_syst->SetMarkerColor (iCent == 0 ? kBlack : finalColors[iPtZ-1]);
        g_syst->SetLineColor (iCent == 0 ? kBlack : finalColors[iPtZ-1]);
        g_syst->SetFillColorAlpha (iCent == 0 ? kBlack : finalFillColors[iPtZ-1], 0.3);

        ((TGAE*) g_syst->Clone ())->Draw ("5P");

        SaferDelete (&g_syst);
      } // end loop over iCent

      for (short iCent : {0, numCentBins}) {
        TGAE* g_stat = make_graph (stats[iPad%2][iPtZ][iCent]);

        Style_t markerStyle = (iCent == 0 ? kOpenCircle : markerStyles[iPtZ-2]);
        float markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);

        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        g_stat->SetLineWidth (3);
        g_stat->SetMarkerColor (iCent == 0 ? kBlack : finalColors[iPtZ-1]);
        g_stat->SetLineColor (iCent == 0 ? kBlack : finalColors[iPtZ-1]);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        markerStyle = FullToOpenMarker (markerStyle);
        markerSize = (markerStyle == kFullDiamond || markerStyle == kOpenDiamond ? 2.4 : 1.8);
        
        g_stat->SetMarkerStyle (markerStyle);
        g_stat->SetMarkerSize (markerSize);
        if (iCent > 0) g_stat->SetLineWidth (0);
        g_stat->SetMarkerColor (kBlack);

        ((TGAE*) g_stat->Clone ())->Draw ("P");

        SaferDelete (&g_stat);
      } // end loop over iCent
    } // end loop over iPad

    tl->SetTextColor (kBlack);

    tl->SetTextAlign (22);
    tl->SetTextSize (34);
    ulPad->cd ();
    tl->DrawLatexNDC (0.5*(1.+lMargin/padX), 0.938, "2 < #it{p}_{T}^{ch} < 4 GeV");
    urPad->cd ();
    tl->DrawLatexNDC (0.5*(1.-rMargin/(1.-padX)), 0.938, "#it{p}_{T}^{ch} > 4 GeV");

    tl->SetTextAlign (11);

    ulPad->cd ();
    tl->SetTextSize (40);
    tl->DrawLatexNDC (0.23, 0.750, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.23, 0.660, "#it{p}_{T}^{#it{Z}} > 60 GeV");
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.32, (0.580-0.020), "#it{pp}");
    tl->DrawLatexNDC (0.32, (0.500-0.020), "0-30% Pb+Pb");
    MakeDataBox   (0.33, (0.580), finalFillColors[0], 0.30, kOpenCircle, 1.8, 1., 0.6/(1.-yPadCUMiddle));
    MakeDataBox   (0.33, (0.500), finalFillColors[3], 0.30, markerStyles[2], 2.4, 1., 0.6/(1.-yPadCUMiddle));

    clPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.23, (0.660+0.09)/(1.-tMargin), "30 < #it{p}_{T}^{#it{Z}} < 60 GeV");
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.32, (0.580-0.020+0.09)/(1.-tMargin), "#it{pp}");
    tl->DrawLatexNDC (0.32, (0.500-0.020+0.09)/(1.-tMargin), "0-30% Pb+Pb");
    MakeDataBox   (0.33, (0.580+0.09)/(1.-tMargin), finalFillColors[0], 0.30, kOpenCircle, 1.8, 1., 0.6/(yPadCUMiddle-yPadDCMiddle));
    MakeDataBox   (0.33, (0.500+0.09)/(1.-tMargin), finalFillColors[2], 0.30, markerStyles[1], 1.8, 1., 0.6/(yPadCUMiddle-yPadDCMiddle));

    dlPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.23, (0.660+0.09)/(1.-tMargin) * (1.-bMargin) + bMargin, "15 < #it{p}_{T}^{#it{Z}} < 30 GeV");
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.32, (0.580-0.020+0.09)/(1.-tMargin) * (1.-bMargin) + bMargin, "#it{pp}");
    tl->DrawLatexNDC (0.32, (0.500-0.020+0.09)/(1.-tMargin) * (1.-bMargin) + bMargin, "0-30% Pb+Pb");
    MakeDataBox   (0.33, (0.580+0.09)/(1.-tMargin) * (1.-bMargin) + bMargin, finalFillColors[0], 0.30, kOpenCircle, 1.8, 1., 0.6/yPadDCMiddle);
    MakeDataBox   (0.33, (0.500+0.09)/(1.-tMargin) * (1.-bMargin) + bMargin, finalFillColors[1], 0.30, markerStyles[0], 1.8, 1., 0.6/yPadDCMiddle);

    urPad->cd ();
    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.07, 0.760, "#it{pp}, #sqrt{s} = 5.02 TeV");
    tl->DrawLatexNDC (0.07, 0.660, "260 pb^{-1}");

    crPad->cd ();
    tl->DrawLatexNDC (0.07, 0.760/(1.-tMargin), "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV");
    tl->DrawLatexNDC (0.07, 0.660/(1.-tMargin), "1.4-1.7 nb^{-1}");

    c26->SaveAs (Form ("%s/FinalPlots/yield_dphi_allptz.pdf", plotPath.Data ()));

  }




  {

    ofstream table1File;
    table1File.open (Form ("%s/table1.tex", tablesPath.Data ()));
  
    double x, y, x_err, y_stat, y_syst;

    table1File << "\\begin{table}[!ht]" << endl;
    table1File << "\\begin{center}" << endl;
    table1File << "\\hspace*{-0.75cm}" << endl;
    table1File << "\\renewcommand{\\arraystretch}{1.2}" << endl;
    table1File << "\\begin{tabular}{|l|l|l|l|l|}" << endl;

    table1File << "\\hline" << endl;
    table1File << "\\multicolumn{1}{|c|}{\\multirow{2}{*}{\\ptch [\\GeV]}} & \\multicolumn{4}{c|}{$(1/N_\\mathrm{Z})(d^2N_\\mathrm{ch}/d\\ptch d\\Delta\\phi)$ $\\pm$ (Stat. Unc.) $\\pm$ (Syst. Unc.) [\\GeV$^{-1}$]} \\\\ \\cline{2-5}" << endl;
    table1File << " & \\multicolumn{1}{c|}{\\small \\pp} & \\multicolumn{1}{c|}{\\small 30--80\\% \\PbPb} & \\multicolumn{1}{c|}{\\small 10--30\\% \\PbPb} & \\multicolumn{1}{c|}{\\small 0--10\\% \\PbPb} \\\\ \\hline \\hline" << endl;
      
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {

      const char* ptzStr = (iPtZ == nPtZBins-1 ? Form ("{$\\iptz > \\SI{%g}{\\GeV}$}", zPtBins[iPtZ]) : Form ("{%g < \\iptz < $\\SI{%g}{\\GeV}$}", zPtBins[iPtZ], zPtBins[iPtZ+1]));

      table1File << "\\multicolumn{5}{|c|}{\\small " << ptzStr << "} \\\\ \\hline" << endl;

      for (int i = 0; i < nPtchBins[iPtZ]; i++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          TGAE* g_syst = g_trk_pt_ptz_sub_syst[iPtZ][iCent];
          TH1D* h_stat = h_trk_pt_ptz_sub_stat[iPtZ][iCent];
          x_err = g_syst->GetErrorX (i);
          y_syst = g_syst->GetErrorY (i);
          x = h_stat->GetBinCenter (i+1);
          y = h_stat->GetBinContent (i+1);
          y_stat = h_stat->GetBinError (i+1);

          string s_y = to_string (y);
          string s_y_stat = to_string (y_stat);
          string s_y_syst = to_string (y_syst);

          FormatMeasurement (s_y, s_y_stat, s_y_syst, 2);

          if (iCent == 0) table1File << "\\footnotesize {$" << x-x_err << " - " << x+x_err << "$} & ";
          table1File << "\\scriptsize {$" << s_y << " \\pm " << s_y_stat << " \\pm " << s_y_syst << "$} ";
          if (iCent < numCentBins-1) table1File << " & ";
        } // end loop over iCent
        table1File << " \\\\" << endl;
      }
      table1File << "\\hline";
      if (iPtZ < nPtZBins-1) table1File << " \\hline" << endl;
    } // end loop over iPtZ

    table1File << endl << "\\end{tabular}" << endl;
    table1File << "\\caption{Summary of UE subtracted per-\\Zboson charged particle yields $(1/N_Z)(d^2N_\\mathrm{ch}/d\\ptch d\\Delta\\phi)$, with total statistical and systematic uncertainties.}" << endl;
    table1File << "\\label{tab1}" << endl;
    table1File << "\\end{center}" << endl;
    table1File << "\\end{table}" << endl;

    table1File.close ();
  }




  {

    ofstream table2File;
    table2File.open (Form ("%s/table2.tex", tablesPath.Data ()));
  
    double x, y, x_err, y_stat, y_syst;

    table2File << "\\begin{table}[!ht]" << endl;
    table2File << "\\begin{center}" << endl;
    table2File << "\\renewcommand{\\arraystretch}{1.2}" << endl;
    table2File << "\\begin{tabular}{|l|l|l|l|l|}" << endl;

    table2File << "\\hline" << endl;
    table2File << "\\multicolumn{1}{|c|}{\\multirow{2}{*}{\\xhz}} & \\multicolumn{4}{c|}{$(1/N_\\mathrm{Z})(d^2N_\\mathrm{ch}/d\\xhz d\\Delta\\phi)$ $\\pm$ (Stat. Unc.) $\\pm$ (Syst. Unc.)} \\\\ \\cline{2-5}" << endl;
    table2File << " & \\multicolumn{1}{c|}{\\small \\pp} & \\multicolumn{1}{c|}{\\small 30--80\\% \\PbPb} & \\multicolumn{1}{c|}{\\small 10--30\\% \\PbPb} & \\multicolumn{1}{c|}{\\small 0--10\\% \\PbPb} \\\\ \\hline \\hline" << endl;
      
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {

      const char* ptzStr = (iPtZ == nPtZBins-1 ? Form ("{$\\iptz > \\SI{%g}{\\GeV}$}", zPtBins[iPtZ]) : Form ("{%g < \\iptz < $\\SI{%g}{\\GeV}$}", zPtBins[iPtZ], zPtBins[iPtZ+1]));

      table2File << "\\multicolumn{5}{|c|}{\\small " << ptzStr << "} \\\\ \\hline" << endl;

      for (int i = 0; i < nXhZBins[iPtZ] - (iPtZ == 2 ? 1 : 0); i++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          TGAE* g_syst = g_trk_xhz_ptz_sub_syst[iPtZ][iCent];
          TH1D* h_stat = h_trk_xhz_ptz_sub_stat[iPtZ][iCent];

          x_err = g_syst->GetErrorX (i);
          y_syst = g_syst->GetErrorY (i);
          x = h_stat->GetBinCenter (i+1);
          y = h_stat->GetBinContent (i+1);
          y_stat = h_stat->GetBinError (i+1);

          string s_y = to_string (y);
          string s_y_stat = to_string (y_stat);
          string s_y_syst = to_string (y_syst);

          FormatMeasurement (s_y, s_y_stat, s_y_syst, 2);

          int xlo = round (1./(x-x_err));
          int xhi = round (1./(x+x_err));

          if (iCent == 0) {
            table2File << "\\footnotesize {$1/" << xlo;
            if (xhi > 1) table2File << " - 1/";
            else table2File << " - ";
            table2File << xhi << "$} & ";
          }

          table2File << "\\scriptsize {$" << s_y << " \\pm " << s_y_stat << " \\pm " << s_y_syst << "$} ";
          if (iCent < numCentBins-1) table2File << " & ";
        } // end loop over iCent
        table2File << " \\\\" << endl;
      }
      table2File << "\\hline";
      if (iPtZ < nPtZBins-1) table2File << " \\hline" << endl;
    } // end loop over iPtZ

    table2File << endl << "\\end{tabular}" << endl;
    table2File << "\\caption{Summary of UE subtracted per-\\Zboson charged particle yields $(1/N_Z)(d^2N_\\mathrm{ch}/d\\xhz d\\Delta\\phi)$, with total statistical and systematic uncertainties.}" << endl;
    table2File << "\\label{tab2}" << endl;
    table2File << "\\end{center}" << endl;
    table2File << "\\end{table}" << endl;

    table2File.close ();
  }




  {
    ofstream table3File;
    table3File.open (Form ("%s/table3.tex", tablesPath.Data ()));

    double x, y, x_err, y_stat, y_syst;

    table3File << "\\begin{table}[!ht]" << endl;
    table3File << "\\begin{center}" << endl;
    table3File << "\\renewcommand{\\arraystretch}{1.2}" << endl;
    table3File << "\\begin{tabular}{|l|l|l|l|}" << endl;

    table3File << "\\hline" << endl;
    table3File << "\\multicolumn{1}{|c|}{\\multirow{2}{*}{\\ptch [\\GeV]}} & \\multicolumn{3}{c|}{\\IAA $\\pm$ (Stat. Unc.) $\\pm$ (Syst. Unc.)} \\\\ \\cline{2-4}" << endl;
    table3File << " & \\multicolumn{1}{c|}{\\small 30--80\\% \\PbPb \\,/ \\pp} & \\multicolumn{1}{c|}{\\small 10--30\\% \\PbPb \\,/ \\pp} & \\multicolumn{1}{c|}{\\small 0--10\\% \\PbPb \\,/ \\pp} \\\\ \\hline \\hline" << endl;
      
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {

      const char* ptzStr = (iPtZ == nPtZBins-1 ? Form ("{$\\iptz > \\SI{%g}{\\GeV}$}", zPtBins[iPtZ]) : Form ("{%g < \\iptz < $\\SI{%g}{\\GeV}$}", zPtBins[iPtZ], zPtBins[iPtZ+1]));

      table3File << "\\multicolumn{4}{|c|}{\\small " << ptzStr << "} \\\\ \\hline" << endl;

      for (int i = 0; i < nPtchBins[iPtZ]; i++) {
        for (short iCent = 1; iCent < numCentBins; iCent++) {
          TGAE* g_syst = g_trk_pt_ptz_iaa_syst[iPtZ][iCent];
          TH1D* h_stat = h_trk_pt_ptz_iaa_stat[iPtZ][iCent];
          x_err = g_syst->GetErrorX (i);
          y_syst = g_syst->GetErrorY (i);
          x = h_stat->GetBinCenter (i+1);
          y = h_stat->GetBinContent (i+1);
          y_stat = h_stat->GetBinError (i+1);

          string s_y = to_string (y);
          string s_y_stat = to_string (y_stat);
          string s_y_syst = to_string (y_syst);

          FormatMeasurement (s_y, s_y_stat, s_y_syst, 2);

          if (iCent == 1) table3File << "\\footnotesize {$" << x-x_err << " - " << x+x_err << "$} & ";
          table3File << "\\footnotesize {$" << s_y << " \\pm " << s_y_stat << " \\pm " << s_y_syst << "$} ";
          if (iCent < numCentBins-1) table3File << " & ";
        } // end loop over iCent
        table3File << " \\\\" << endl;
      }
      table3File << "\\hline";
      if (iPtZ < nPtZBins-1) table3File << " \\hline" << endl;
    } // end loop over iPtZ

    table3File << endl << "\\end{tabular}" << endl;
    table3File << "\\caption{Summary of \\PbPb to \\pp ratios of UE subtracted per-\\Zboson charged particle yields $\\IAA (\\ptch)$, with total statistical and systematic uncertainties.}" << endl;
    table3File << "\\label{tab3}" << endl;
    table3File << "\\end{center}" << endl;
    table3File << "\\end{table}" << endl;

    table3File.close ();
  }




  {
    ofstream table4File;
    table4File.open (Form ("%s/table4.tex", tablesPath.Data ()));

    double x, y, x_err, y_stat, y_syst;

    table4File << "\\begin{table}[!ht]" << endl;
    table4File << "\\begin{center}" << endl;
    table4File << "\\renewcommand{\\arraystretch}{1.2}" << endl;
    table4File << "\\begin{tabular}{|l|l|l|l|}" << endl;

    table4File << "\\hline" << endl;
    table4File << "\\multicolumn{1}{|c|}{\\multirow{2}{*}{\\xhz}} & \\multicolumn{3}{c|}{\\IAA $\\pm$ (Stat. Unc.) $\\pm$ (Syst. Unc.)} \\\\ \\cline{2-4}" << endl;
    table4File << " & \\multicolumn{1}{c|}{\\small 30--80\\% \\PbPb \\,/ \\pp} & \\multicolumn{1}{c|}{\\small 10--30\\% \\PbPb \\,/ \\pp} & \\multicolumn{1}{c|}{\\small 0--10\\% \\PbPb \\,/ \\pp} \\\\ \\hline \\hline" << endl;
      
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {

      const char* ptzStr = (iPtZ == nPtZBins-1 ? Form ("{$\\iptz > \\SI{%g}{\\GeV}$}", zPtBins[iPtZ]) : Form ("{%g < \\iptz < $\\SI{%g}{\\GeV}$}", zPtBins[iPtZ], zPtBins[iPtZ+1]));

      table4File << "\\multicolumn{4}{|c|}{\\small " << ptzStr << "} \\\\ \\hline" << endl;

      for (int i = 0; i < nXhZBins[iPtZ] - (iPtZ == 2 ? 1 : 0); i++) {
        for (short iCent = 1; iCent < numCentBins; iCent++) {
          TGAE* g_syst = g_trk_xhz_ptz_iaa_syst[iPtZ][iCent];
          TH1D* h_stat = h_trk_xhz_ptz_iaa_stat[iPtZ][iCent];

          x_err = g_syst->GetErrorX (i);
          y_syst = g_syst->GetErrorY (i);
          x = h_stat->GetBinCenter (i+1);
          y = h_stat->GetBinContent (i+1);
          y_stat = h_stat->GetBinError (i+1);

          string s_y = to_string (y);
          string s_y_stat = to_string (y_stat);
          string s_y_syst = to_string (y_syst);

          FormatMeasurement (s_y, s_y_stat, s_y_syst, 2);

          int xlo = round (1./(x-x_err));
          int xhi = round (1./(x+x_err));

          if (iCent == 1) {
            table4File << "\\footnotesize {$1/" << xlo;
            if (xhi > 1) table4File << " - 1/";
            else table4File << " - ";
            table4File << xhi << "$} & ";
          }
          table4File << "\\footnotesize {$" << s_y << " \\pm " << s_y_stat << " \\pm " << s_y_syst << "$} ";
          if (iCent < numCentBins-1) table4File << " & ";
        } // end loop over iCent
        table4File << " \\\\" << endl;
      }
      table4File << "\\hline";
      if (iPtZ < nPtZBins-1) table4File << " \\hline" << endl;
    } // end loop over iPtZ

    table4File << endl << "\\end{tabular}" << endl;
    table4File << "\\caption{Summary of \\PbPb to \\pp ratios of UE subtracted per-\\Zboson charged particle yields $\\IAA (\\xhz)$, with total statistical and systematic uncertainties.}" << endl;
    table4File << "\\label{tab4}" << endl;
    table4File << "\\end{center}" << endl;
    table4File << "\\end{table}" << endl;

    table4File.close ();
  }




  {
    ofstream table5File;
    table5File.open (Form ("%s/table5.tex", tablesPath.Data ()));

    double x, ntrk, x_err, ntrk_stat, ntrk_syst, ptch, ptch_stat, ptch_syst, xhz, xhz_stat, xhz_syst;

    table5File << "\\begin{table}[!ht]" << endl;
    table5File << "\\begin{center}" << endl;
    table5File << "\\renewcommand{\\arraystretch}{1.2}" << endl;
    table5File << "\\begin{tabular}{|l|l|l|l|}" << endl;
    table5File << "\\hline" << endl;      
    table5File << "\\multirow{2}{*}{$\\left<\\iptz\\right>\\pm$ (Stat. Unc.) [\\GeV]} & \\multicolumn{3}{c|}{$\\left<\\text{Value}\\right>$ $\\pm$ (Stat. Unc.) $\\pm$ (Syst. Unc.)}" << " \\\\ \\cline{2-4}" << endl;
    table5File << "& \\multicolumn{1}{c|}{$N_\\mathrm{ch}$} & \\multicolumn{1}{c|}{\\ptch [\\GeV]} & \\multicolumn{1}{c|}{\\xhz} \\\\" << endl;
    for (short iCent = 0; iCent < numCentBins; iCent++) {

      const char* centStr = (iCent == 0 ? "\\pp" : Form ("%i--%i\\%% \\PbPb", (int)centCuts[iCent], (int)centCuts[iCent-1]));
      table5File << "\\hline \\hline \\multicolumn{4}{|c|}{" << centStr << "} \\\\ \\hline" << endl;
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {

        TGAE* g_ntrk_stat = g_avg_ntrk_ptz_stat[iCent];
        TGAE* g_ntrk_syst = g_avg_ntrk_ptz_syst[iCent];
        TGAE* g_ptch_stat = g_trk_avg_pt_ptz_stat[iCent];
        TGAE* g_ptch_syst = g_trk_avg_pt_ptz_syst[iCent];
        TGAE* g_xhz_stat = g_trk_avg_xhz_ptz_stat[iCent];
        TGAE* g_xhz_syst = g_trk_avg_xhz_ptz_syst[iCent];

        g_ntrk_stat->GetPoint (iPtZ-2, x, ntrk);
        x_err = g_ntrk_stat->GetErrorX (iPtZ-2);
        ntrk_stat = g_ntrk_stat->GetErrorY (iPtZ-2);
        ntrk_syst = g_ntrk_syst->GetErrorY (iPtZ-2);
        g_ptch_stat->GetPoint (iPtZ-2, x, ptch);
        ptch_stat = g_ptch_stat->GetErrorY (iPtZ-2);
        ptch_syst = g_ptch_syst->GetErrorY (iPtZ-2);
        g_xhz_stat->GetPoint (iPtZ-2, x, xhz);
        xhz_stat = g_xhz_stat->GetErrorY (iPtZ-2);
        xhz_syst = g_xhz_syst->GetErrorY (iPtZ-2);

        string s_x = to_string (x);
        string s_x_err = to_string (x_err);
        string dummy = to_string (x_err); FormatMeasurement (s_x, s_x_err, dummy, 2);

        string s_ntrk = to_string (ntrk);
        string s_ntrk_stat = to_string (ntrk_stat);
        string s_ntrk_syst = to_string (ntrk_syst);
        string s_ptch = to_string (ptch);
        string s_ptch_stat = to_string (ptch_stat);
        string s_ptch_syst = to_string (ptch_syst);
        string s_xhz = to_string (xhz);
        string s_xhz_stat = to_string (xhz_stat);
        string s_xhz_syst = to_string (xhz_syst);

        FormatMeasurement (s_ntrk, s_ntrk_stat, s_ntrk_syst, 2);
        FormatMeasurement (s_ptch, s_ptch_stat, s_ptch_syst, 2);
        FormatMeasurement (s_xhz, s_xhz_stat, s_xhz_syst, 2);

        table5File << "{\\footnotesize $" << s_x << " \\pm " << s_x_err << "$} & \\footnotesize {$" << s_ntrk << " \\pm " << s_ntrk_stat << " \\pm " << s_ntrk_syst << "$} ";
        table5File << "& \\footnotesize {$" << s_ptch << " \\pm " << s_ptch_stat << " \\pm " << s_ptch_syst << "$} ";
        table5File << "& \\footnotesize {$" << s_xhz << " \\pm " << s_xhz_stat << " \\pm " << s_xhz_syst << "$} ";
        table5File << " \\\\" << endl;
      } // end loop over iPtZ
    } // end loop over iCent

    table5File << "\\hline" << endl;
    table5File << "\\end{tabular}" << endl;
    table5File << "\\caption{Summary of mean number of charged particles $N_\\mathrm{ch}$, \\ptch, and \\xhz, with total statistical and systematic uncertainties, for each collision system and centrality.}" << endl;
    table5File << "\\label{tab5}" << endl;
    table5File << "\\end{center}" << endl;
    table5File << "\\end{table}" << endl;

    table5File.close ();
  }




  //{

  //  ofstream table6File;
  //  table6File.open (Form ("%s/table6.tex", tablesPath.Data ()));
  //
  //  double x, y, x_err, y_stat, y_syst;

  //  table6File << "\\hline" << endl;
  //  table6File << "\\multicolumn{1}{|c|}{\\multirow{2}{*}{$\\Delta\\phi_\\mathrm{h#it{Z}}$}} & \\multicolumn{2}{c|}{$(1/N_\\mathrm{Z})(dN_\\mathrm{ch}/dd\\Delta\\phi_\\mathrm{h#it{Z}})$ $\\pm$ (Stat. Unc.) $\\pm$ (Syst. Unc.) [\\GeV$^{-1}$]} \\\\ \\cline{2-5}" << endl;
  //  table6File << " & \\multicolumn{1}{c|}{\\small \\pp} & \\multicolumn{1}{c|}{\\small 0-30\\% \\PbPb} \\\\ \\hline \\hline" << endl;
  //    
  //  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {

  //    const char* ptzStr = (iPtZ == nPtZBins-1 ? Form ("{$\\iptz > \\SI{%g}{\\GeV}$}", zPtBins[iPtZ]) : Form ("{%g < \\iptz < $\\SI{%g}{\\GeV}$}", zPtBins[iPtZ], zPtBins[iPtZ+1]));

  //    table6File << "\\multicolumn{3}{|c|}{\\small " << ptzStr << "} \\\\ \\hline" << endl;

  //    for (int i = 0; i < g_trk_dphi_ptz_gt4_sub_syst[iPtZ][0]->GetN (); i++) {
  //      for (short iCent : {0, numCentBins}) {
  //        TGAE* g_syst = g_trk_dphi_ptz_gt4_sub_syst[iPtZ][iCent];
  //        TH1D* h_stat = h_trk_dphi_ptz_gt4_sub_stat[iPtZ][iCent];
  //        x_err = g_syst->GetErrorX (i);
  //        y_syst = g_syst->GetErrorY (i);
  //        x = h_stat->GetBinCenter (i+1);
  //        y = h_stat->GetBinContent (i+1);
  //        y_stat = h_stat->GetBinError (i+1);

  //        string s_y = to_string (y);
  //        string s_y_stat = to_string (y_stat);
  //        string s_y_syst = to_string (y_syst);

  //        double x_lo = round (20.*(x-x_err)/pi)/20.;
  //        double x_hi = round (20.*(x+x_err)/pi)/20.;

  //        string s_x_lo = (x_lo == 0 ? "0" : to_string (x_lo));
  //        string s_x_hi = (x_hi == 1. ? "" : to_string (x_hi));
  //        while (s_x_lo.back () == '0' && s_x_lo.length () > 1)
  //          s_x_lo.pop_back ();
  //        while (s_x_hi.back () == '0' && s_x_hi.length () > 1)
  //          s_x_hi.pop_back ();
  //        s_x_hi = s_x_hi + "\\pi";

  //        FormatMeasurement (s_y, s_y_stat, s_y_syst, 2);

  //        if (iCent == 0) table6File << "\\footnotesize {$" << s_x_lo << " - " << s_x_hi << "$} & ";
  //        table6File << "\\scriptsize {$" << s_y << " \\pm " << s_y_stat << " \\pm " << s_y_syst << "$} ";
  //        if (iCent < numCentBins-1) table6File << " & ";
  //      } // end loop over iCent
  //      table6File << " \\\\" << endl;
  //    }
  //    table6File << "\\hline";
  //    if (iPtZ < nPtZBins-1) table6File << " \\hline" << endl;
  //  } // end loop over iPtZ
  //  table6File.close ();
  //}

} // end of macro

#endif

