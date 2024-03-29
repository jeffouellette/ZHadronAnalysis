#ifndef __PlotTrackMomentumResolution_C__
#define __PlotTrackMomentumResolution_C__

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF1.h>
#include <TLatex.h>
#include <TLine.h>
#include <TCanvas.h>

#include <Utilities.h>

#include "Params.h"

typedef TGraphAsymmErrors TGAE;

const int numFinerEtachBins = 40;
const double* finerEtachBins = linspace (-2.5, 2.5, numFinerEtachBins);

const double etachBins[6] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};
const int numEtachBins = sizeof (etachBins) / sizeof (etachBins[0]) - 1;

const double pTchBinsPP[33] = {0.80, 0.84, 0.88, 0.92, 0.96, 1, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 60, 100};
const int numPtchBinsPP = sizeof (pTchBinsPP) / sizeof (pTchBinsPP[0]) - 1;
const double pTchBinsPbPb[31] = {0.80, 0.84, 0.88, 0.92, 0.96, 1, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 7, 8, 10, 15, 20, 30, 60, 100};
const int numPtchBinsPbPb = sizeof (pTchBinsPbPb) / sizeof (pTchBinsPbPb[0]) - 1;

const bool isPbPb = false;

TH2D** h2_avg_tms = nullptr;
TH2D** h2_avg_tmr = nullptr;
TH1D*** h_avg_tms = nullptr;
TH1D*** h_avg_tmr = nullptr;


void PlotTrackMomentumResolution () {

  TFile* inFile = new TFile (Form ("%s/TrackingMomentumResolution/Nominal/trackingMomentumResolutionFactors.root", rootPath.Data ()), "read");

  h2_avg_tms = new TH2D*[numCentBins];
  h2_avg_tmr = new TH2D*[numCentBins];
  h_avg_tms = new TH1D**[numCentBins];
  h_avg_tmr = new TH1D**[numCentBins];

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    h2_avg_tms[iCent] = (TH2D*) inFile->Get (Form ("h2_avg_tms_iCent%i", iCent));
    h2_avg_tmr[iCent] = (TH2D*) inFile->Get (Form ("h2_avg_tmr_iCent%i", iCent));

    h_avg_tms[iCent] = new TH1D*[numEtachBins+1];
    h_avg_tmr[iCent] = new TH1D*[numEtachBins+1];
    for (int iEta = 0; iEta <= numEtachBins; iEta++) {
      h_avg_tms[iCent][iEta] = (TH1D*) inFile->Get (Form ("h_avg_tms_iCent%i_iEta%i", iCent, iEta));
      h_avg_tmr[iCent][iEta] = (TH1D*) inFile->Get (Form ("h_avg_tmr_iCent%i_iEta%i", iCent, iEta));
    } // end loop over iEta
  } // end loop over iCent



  TLatex* tl = new TLatex ();
  TLine* l = new TLine ();


  for (short iCent : {0, 3}) {
    TCanvas* c = new TCanvas (Form ("c_tms_iCent%i", iCent), "", 800, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetBottomMargin (bMargin);
    c->SetTopMargin (tMargin);

    {
      gPad->SetLogx ();

      TH1D* htemp = new TH1D ("htemp", "", 1, 0.8, 100);
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelSize (0);

      yax->SetTitle ("TMS [%]");
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 95;
      const double ymax = 105;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
      SaferDelete (&htemp);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      
      const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
      //tl->DrawLatex (0.8,  yoff, "0.8");
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
      //tl->DrawLatex (80, yoff, "80");
      tl->DrawLatex (100, yoff, "100");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      //l->SetLineColor (kPink-8);
      l->DrawLine ((isPbPb ? pTchBinsPP : pTchBinsPbPb)[0], 100, (isPbPb ? pTchBinsPP : pTchBinsPbPb)[(isPbPb ? numPtchBinsPP : numPtchBinsPbPb)], 100);
    }
     

    for (int iEta = 0; iEta < numEtachBins; iEta++) {
      TGAE* g = make_graph (h_avg_tms[iCent][iEta]);

      g->SetMarkerColor (colors[iEta+1]);
      g->SetLineColor (colors[iEta+1]);
      g->SetMarkerStyle (kOpenCircle);
     
      g->Draw ("P"); 
    }
    {
      TGAE* g = make_graph (h_avg_tms[iCent][numEtachBins]);

      g->SetMarkerColor (colors[0]);
      g->SetLineColor (colors[0]);
      g->SetMarkerStyle (kFullCircle);
     
      g->Draw ("P"); 
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    if (iCent == 0) {
      tl->SetTextSize (30);
      tl->DrawLatexNDC (0.26, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV");
      tl->DrawLatexNDC (0.26, 0.800, "Powheg + Pythia #it{Z} #rightarrow #it{ll}");
    }
    else if (iCent == 3) {
      tl->SetTextSize (30);
      tl->DrawLatexNDC (0.26, 0.845, "Pb+Pb, #sqrt{s} = 5.02 TeV, 0-10%");
      tl->DrawLatexNDC (0.26, 0.800, "Powheg + Pythia #it{Z} #rightarrow #it{ll} + Hijing overlay");
    }
    for (int iEta = 0; iEta < numEtachBins; iEta++)
      myMarkerTextNoLine (0.3, 0.39-0.04*iEta, colors[iEta+1], kOpenCircle, Form ("%g < |#it{#eta}| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.04);
    myMarkerTextNoLine (0.3, 0.44, colors[0], kFullCircle, "All #it{#eta}", 1.2, 0.04);

    c->SaveAs (Form ("../Plots/TrackMomentumResolutionStudy/TrackMomentumScale_iCent%i.pdf", iCent));
  } // end loop over iCent



  for (short iCent : {0, 3}) {
    TCanvas* c = new TCanvas (Form ("c_tmr_iCent%i", iCent), "", 800, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetBottomMargin (bMargin);
    c->SetTopMargin (tMargin);

    {
      gPad->SetLogx ();

      TH1D* htemp = new TH1D ("htemp", "", 1, 0.8, 100);
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelSize (0);

      yax->SetTitle ("TMR [%]");
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 0;
      const double ymax = 20;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
      SaferDelete (&htemp);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      
      const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
      //tl->DrawLatex (0.8,  yoff, "0.8");
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
      //tl->DrawLatex (80, yoff, "80");
      tl->DrawLatex (100, yoff, "100");
    }
     

    for (int iEta = 0; iEta < numEtachBins; iEta++) {
      TGAE* g = make_graph (h_avg_tmr[iCent][iEta]);

      g->SetMarkerColor (colors[iEta+1]);
      g->SetLineColor (colors[iEta+1]);
      g->SetMarkerStyle (kOpenCircle);
     
      g->Draw ("P"); 
    }
    {
      TGAE* g = make_graph (h_avg_tmr[iCent][numEtachBins]);

      g->SetMarkerColor (colors[0]);
      g->SetLineColor (colors[0]);
      g->SetMarkerStyle (kFullCircle);
     
      g->Draw ("P"); 
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    if (iCent == 0) {
      tl->SetTextSize (30);
      tl->DrawLatexNDC (0.26, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV");
      tl->DrawLatexNDC (0.26, 0.800, "Powheg + Pythia #it{Z} #rightarrow #it{ll}");
    }
    else if (iCent == 3) {
      tl->SetTextSize (30);
      tl->DrawLatexNDC (0.26, 0.845, "Pb+Pb, #sqrt{s} = 5.02 TeV, 0-10%");
      tl->DrawLatexNDC (0.26, 0.800, "Powheg + Pythia #it{Z} #rightarrow #it{ll} + Hijing overlay");
    }

    for (int iEta = 0; iEta < numEtachBins; iEta++)
      myMarkerTextNoLine (0.5, 0.69-0.04*iEta, colors[iEta+1], kOpenCircle, Form ("%g < |#it{#eta}| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.04);
    myMarkerTextNoLine (0.5, 0.74, colors[0], kFullCircle, "All #it{#eta}", 1.2, 0.04);

    c->SaveAs (Form ("../Plots/TrackMomentumResolutionStudy/TrackMomentumResolution_iCent%i.pdf", iCent));
  }


}

#endif
