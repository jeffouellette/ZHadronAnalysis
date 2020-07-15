#ifndef __LowPtVariation_h__
#define __LowPtVariation_h__

#include "Params.h"
#include "PhysicsAnalysis.h"
#include "Systematic.h"
#include "ReweightingVariation.h"

#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;

class LowPtVariation : ReweightingVariation {

  public:

  LowPtVariation (const char* _name = "lowPtVariation") : ReweightingVariation (_name) { }

  virtual void GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* bkg) override;
  virtual void ApplyRelativeVariations (PhysicsAnalysis* nominal, const bool upVar = true) override;
};


void LowPtVariation :: GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* bkg) {

  cout << "Calculating low pT variation weights from subtracted yields." << endl;

  TGAE** g_pTch = Get1DArray <TGAE*> (numCentBins);
  //TGAE** g_xhZ = Get1DArray <TGAE*> (numCentBins);

  float max_pp_sys = 0;

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    g_pTch[iCent] = new TGAE ();
    g_pTch[iCent]->SetName (Form ("g_pTch_iCent%i", iCent));
    //g_xhZ[iCent] = new TGAE ();
    //g_xhZ[iCent]->SetName (Form ("g_xhZ_iCent%i", iCent));
  } // end loop over iCent

  TH1D* h_ee = nullptr, *h_mumu = nullptr, *h_comb = nullptr, *h_ee_bkg = nullptr, *h_mumu_bkg = nullptr;
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    const double wee = nominal->h_z_counts[0][2][iCent]->GetBinContent (2) / nominal->h_z_counts[2][2][iCent]->GetBinContent (2);
    const double wmm = nominal->h_z_counts[1][2][iCent]->GetBinContent (2) / nominal->h_z_counts[2][2][iCent]->GetBinContent (2);    

    {
      TGAE* g = g_pTch[iCent];

      h_ee = nominal->h_trk_pt_ptz_sub[0][2][iCent];
      h_mumu = nominal->h_trk_pt_ptz_sub[1][2][iCent];
      h_comb = nominal->h_trk_pt_ptz_sub[2][2][iCent];
      h_ee_bkg = bkg->h_trk_pt_ptz[0][2][iCent];
      h_mumu_bkg = bkg->h_trk_pt_ptz[1][2][iCent];
      assert (h_ee->GetNbinsX () == h_mumu->GetNbinsX () && h_ee->GetNbinsX () == h_ee_bkg->GetNbinsX () && h_ee->GetNbinsX () == h_mumu_bkg->GetNbinsX ());
      int nBinsX = h_ee->GetNbinsX ();

      for (int ix = 1; ix <= nBinsX; ix++) {
        if (h_ee->GetBinCenter (ix) < pTchBins[2][0] || pTchBins[2][nPtchBins[2]] < h_ee->GetBinCenter (ix))
          continue;

        const double yee = h_ee->GetBinContent (ix);
        const double yee_err = sqrt (pow (h_ee->GetBinError (ix), 2) + pow (h_ee_bkg->GetBinError (ix), 2));
        const double ymm = h_mumu->GetBinContent (ix);
        const double ymm_err = sqrt (pow (h_mumu->GetBinError (ix), 2) + pow (h_mumu_bkg->GetBinError (ix), 2));
        const double ycomb = h_comb->GetBinContent (ix);

        assert (yee != 0 && ymm != 0 && ycomb != 0);

        if (iCent == 0) {
          const double fee = (yee-ycomb) / ycomb;
          const double fee_err = sqrt (pow (wmm*ymm*yee_err / pow (ycomb, 2), 2) + pow (-wmm*yee*ymm_err / pow (ycomb, 2), 2));
          const double fmm = (ymm-ycomb) / ycomb;
          const double fmm_err = sqrt (pow (wee*yee*ymm_err / pow (ycomb, 2), 2) + pow (-wee*ymm*yee_err / pow (ycomb, 2), 2));
          if (iCent == 0)
            max_pp_sys = fmax (max_pp_sys, fmax (fabs (fee), fabs(fmm)));
        }
        const double f = (yee-ymm) / ycomb;
        const double f_err = sqrt (pow (ymm*yee_err/(ycomb*ycomb), 2) + pow (yee*ymm_err/(ycomb*ycomb), 2));

        g->SetPoint (g->GetN (), h_ee->GetBinCenter (ix), f);
        g->SetPointEYhigh (g->GetN () - 1, f_err);
        g->SetPointEYlow (g->GetN () - 1, f_err);
        g->SetPointEXhigh (g->GetN () - 1, 0.5 * h_ee->GetBinWidth (ix));
        g->SetPointEXlow (g->GetN () - 1, 0.5 * h_ee->GetBinWidth (ix));
      } // end loop over ix
    }

    //{
    //  TGAE* g = g_xhZ[iCent];

    //  h_ee = nominal->h_trk_xhz_ptz_sub[0][2][iCent];
    //  h_mumu = nominal->h_trk_xhz_ptz_sub[1][2][iCent];
    //  h_comb = nominal->h_trk_xhz_ptz_sub[2][2][iCent];
    //  h_ee_bkg = bkg->h_trk_xhz_ptz[0][2][iCent];
    //  h_mumu_bkg = bkg->h_trk_xhz_ptz[1][2][iCent];
    //  assert (h_ee->GetNbinsX () == h_mumu->GetNbinsX () && h_ee->GetNbinsX () == h_ee_bkg->GetNbinsX () && h_ee->GetNbinsX () == h_mumu_bkg->GetNbinsX ());
    //  int nBinsX = h_ee->GetNbinsX ();

    //  for (int ix = 1; ix <= nBinsX; ix++) {
    //    if (h_ee->GetBinCenter (ix) < xhZBins[2][0] || xhZBins[2][nXhZBins[2]] < h_ee->GetBinCenter (ix))
    //      continue;

    //    const double yee = h_ee->GetBinContent (ix);
    //    const double yee_err = sqrt (pow (h_ee->GetBinError (ix), 2) + pow (h_ee_bkg->GetBinError (ix), 2));
    //    const double ymm = h_mumu->GetBinContent (ix);
    //    const double ymm_err = sqrt (pow (h_mumu->GetBinError (ix), 2) + pow (h_mumu_bkg->GetBinError (ix), 2));
    //    const double ycomb = h_comb->GetBinContent (ix);

    //    assert (yee != 0 && ymm != 0 && ycomb != 0);

    //    const double f = (yee-ymm) / ycomb;
    //    const double f_err = sqrt (pow (ymm*yee_err/(ycomb*ycomb), 2) + pow (yee*ymm_err/(ycomb*ycomb), 2));

    //    g->SetPoint (g->GetN (), h_ee->GetBinCenter (ix), f);
    //    g->SetPointEYhigh (g->GetN () - 1, f_err);
    //    g->SetPointEYlow (g->GetN () - 1, f_err);
    //    g->SetPointEXhigh (g->GetN () - 1, 0.5 * h_ee->GetBinWidth (ix));
    //    g->SetPointEXlow (g->GetN () - 1, 0.5 * h_ee->GetBinWidth (ix));
    //  } // end loop over ix
    //}

  } // end loop over iCent

  TF1** f_ratio_pTch = Get1DArray <TF1*> (numCentBins);
  //TF1** f_ratio_xhZ = Get1DArray <TF1*> (numCentBins);

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    if (iCent == 0) {
      f_ratio_pTch[iCent] = new TF1 (Form ("f_ratio_pTch_iCent%i", iCent), "[0]+[1]*log(x)+[2]*log(x)*log(x)", pTchBins[2][0], pTchBins[2][nPtchBins[2]]);
      f_ratio_pTch[iCent]->SetParameter (0, 0);
      f_ratio_pTch[iCent]->SetParameter (1, 0);
      f_ratio_pTch[iCent]->SetParameter (2, 0);
      //f_ratio_xhZ[iCent] = new TF1 (Form ("f_ratio_xhZ_%s_iCent%i", iCent), "[0]+[1]*log(x)+[2]*log(x)*log(x)", xhZBins[2][0], xhZBins[2][nXhZBins[2]]);
      //f_ratio_xhZ[iCent]->SetParameter (0, 0);
      //f_ratio_xhZ[iCent]->SetParameter (1, 0);
      //f_ratio_xhZ[iCent]->SetParameter (2, 0);
    }
    else {
      f_ratio_pTch[iCent] = new TF1 (Form ("f_ratio_pTch_iCent%i", iCent), "[0]", pTchBins[2][0], pTchBins[2][nPtchBins[2]]);
      f_ratio_pTch[iCent]->SetParameter (0, 0);
      //f_ratio_xhZ[iCent] = new TF1 (Form ("f_ratio_xhZ_iCent%i", iCent), "[0]+[1]*log(x)", xhZBins[2][0], xhZBins[2][nXhZBins[2]]);
      //f_ratio_xhZ[iCent]->SetParameter (0, 0);
      //f_ratio_xhZ[iCent]->SetParameter (1, 0);
    }
    g_pTch[iCent]->Fit (f_ratio_pTch[iCent], "RN0Q");
    //g_xhZ[iCent]->Fit (f_ratio_xhZ[iCent], "RN0Q");
  } // end loop over iCent

  TH1D* h = nullptr;

  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    {
      const short iCent = 0;
      h = (TH1D*) nominal->h_trk_pt_ptz_sub[2][iPtZ][iCent]->Clone (Form ("h_relVarPt_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()));
      relVarPt[2][iPtZ][iCent] = h;
      for (short iX = 1; iX <= h->GetNbinsX (); iX++) {
        const double rUnc = (iPtZ != 2 ? 0 : max_pp_sys);
        h->SetBinContent (iX, rUnc * h->GetBinContent (iX));
        h->SetBinError (iX, 0);
      }

      h = (TH1D*) nominal->h_trk_xhz_ptz_sub[2][iPtZ][iCent]->Clone (Form ("h_relVarX_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()));
      relVarX[2][iPtZ][iCent] = h;
      for (short iX = 1; iX <= h->GetNbinsX (); iX++) {
        const double rUnc = (iPtZ != 2 ? 0 : max_pp_sys);
        //const double rUnc = 0;//(iPtZ != 2 ? 0 : 2.*fabs (f_ratio_pTch[iCent]->Eval (h->GetBinCenter (iX))) / sqrt (12.));
        h->SetBinContent (iX, rUnc * h->GetBinContent (iX));
        h->SetBinError (iX, 0);
      }
    } // end iCent = 0 scope
    for (short iCent = 1; iCent < numCentBins; iCent++) {
      h = (TH1D*) nominal->h_trk_pt_ptz_sub[2][iPtZ][iCent]->Clone (Form ("h_relVarPt_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()));
      relVarPt[2][iPtZ][iCent] = h;
      for (short iX = 1; iX <= h->GetNbinsX (); iX++) {
        const double rUnc = (iPtZ != 2 ? 0 : fabs (f_ratio_pTch[iCent]->Eval (h->GetBinCenter (iX))) / sqrt (12.));
        h->SetBinContent (iX, rUnc * h->GetBinContent (iX));
        h->SetBinError (iX, 0);
      }

      h = (TH1D*) nominal->h_trk_xhz_ptz_sub[2][iPtZ][iCent]->Clone (Form ("h_relVarX_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()));
      relVarX[2][iPtZ][iCent] = h;
      for (short iX = 1; iX <= h->GetNbinsX (); iX++) {
        // systematic uncertainty vs. x_hZ deemed not needed.
        const double rUnc = 0;//(iPtZ != 2 ? 0 : fabs (f_ratio_pTch[iCent]->Eval (h->GetBinCenter (iX))) / sqrt (12.));
        h->SetBinContent (iX, rUnc * h->GetBinContent (iX));
        h->SetBinError (iX, 0);
      }
    } // end loop over iCent
  } // end loop over iPtZ

  //Delete2DArray (f_ratio_pTch, 2, numCentBins);
  //Delete2DArray (f_ratio_xhZ, 2, numCentBins);
  //Delete2DArray (g_pTch, 2, numCentBins);
  //Delete2DArray (g_xhZ, 2, numCentBins);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void LowPtVariation :: ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar) {
  TH1D* rVar = nullptr, *h = nullptr;
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        rVar = relVarPt[2][iPtZ][iCent];
        h = a->h_trk_pt_ptz[iSpc][iPtZ][iCent];
        if (h) {
          assert (h->GetNbinsX () == rVar->GetNbinsX ());
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            h->SetBinContent (iX, h->GetBinContent (iX) + (upVar ? 1 : -1) * rVar->GetBinContent (iX));
          } // end loop over iX
        }
        h = a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent];
        if (h) {
          assert (h->GetNbinsX () == rVar->GetNbinsX ());
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            h->SetBinContent (iX, h->GetBinContent (iX) + (upVar ? 1 : -1) * rVar->GetBinContent (iX));
          } // end loop over iX
        }
        h = a->h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent];
        if (h) {
          TH1D* hpp = a->h_trk_pt_ptz_sub[iSpc][iPtZ][0];
          TH1D* varpp = relVarPt[2][iPtZ][0];
          TH1D* hPbPb = a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent];
          TH1D* varPbPb = relVarPt[2][iPtZ][iCent];

          assert (h->GetNbinsX () == rVar->GetNbinsX ());
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            const double unc = fabs (h->GetBinContent (iX)) * sqrt (pow (varpp->GetBinContent (iX) / hpp->GetBinContent (iX), 2) + pow (varPbPb->GetBinContent (iX) / hPbPb->GetBinContent (iX), 2));
            h->SetBinContent (iX, h->GetBinContent (iX) + (upVar ? 1 : -1) * unc);
          } // end loop over iX
        }

        rVar = relVarX[2][iPtZ][iCent];
        h = a->h_trk_xhz_ptz[iSpc][iPtZ][iCent];
        if (h) {
          assert (h->GetNbinsX () == rVar->GetNbinsX ());
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            h->SetBinContent (iX, h->GetBinContent (iX) + (upVar ? 1 : -1) * rVar->GetBinContent (iX));
          } // end loop over iX
        }
        h = a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent];
        if (h) {
          assert (h->GetNbinsX () == rVar->GetNbinsX ());
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            h->SetBinContent (iX, h->GetBinContent (iX) + (upVar ? 1 : -1) * rVar->GetBinContent (iX));
          } // end loop over iX
        }
        h = a->h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent];
        if (h) {
          TH1D* hpp = a->h_trk_xhz_ptz_sub[iSpc][iPtZ][0];
          TH1D* varpp = relVarX[2][iPtZ][0];
          TH1D* hPbPb = a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent];
          TH1D* varPbPb = relVarX[2][iPtZ][iCent];

          assert (h->GetNbinsX () == rVar->GetNbinsX ());
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            const double unc = fabs (h->GetBinContent (iX)) * sqrt (pow (varpp->GetBinContent (iX) / hpp->GetBinContent (iX), 2) + pow (varPbPb->GetBinContent (iX) / hPbPb->GetBinContent (iX), 2));
            h->SetBinContent (iX, h->GetBinContent (iX) + (upVar ? 1 : -1) * unc);
          } // end loop over iX
        }
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc

  return;
}


#endif
