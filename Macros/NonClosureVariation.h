#ifndef __NonClosureVariation_h__
#define __NonClosureVariation_h__

#include "Params.h"
#include "PhysicsAnalysis.h"

#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;

class NonClosureVariation {

  protected:
  string name;

  public:
  TH1D**** relVarPt = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D**** relVarX  = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D**** relVardPhi = Get3DArray <TH1D*> (nPtZBins, 2, numCentBins+1); // iPtZ, iPtch
  TH1D**** absVardPhi = Get3DArray <TH1D*> (nPtZBins, 2, numCentBins+1); // iPtZ, iPtch

  NonClosureVariation (const char* _name = "nonClosureVariation") {
    name = _name;
  }

  ~NonClosureVariation () {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          SaferDelete (&relVarPt[iSpc][iPtZ][iCent]);
          SaferDelete (&relVarX[iSpc][iPtZ][iCent]);
        } // end loop over iCent
      } // end loop over iPtZ
    } // end loop over iSpc
    Delete3DArray (&relVarPt, 3, nPtZBins, numCentBins);
    Delete3DArray (&relVarX, 3, nPtZBins, numCentBins);

    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
      for (short iCent : {0, numCentBins}) {
        SaferDelete (&relVardPhi[iPtZ][0][iCent]);
        SaferDelete (&relVardPhi[iPtZ][1][iCent]);
        SaferDelete (&absVardPhi[iPtZ][0][iCent]);
        SaferDelete (&absVardPhi[iPtZ][1][iCent]);
      } // end loop over iCent
    } // end loop over iPtZ

    Delete3DArray (&relVardPhi, nPtZBins, 2, numCentBins+1);
    Delete3DArray (&absVardPhi, nPtZBins, 2, numCentBins+1);
  }

  virtual void GetRelativeVariations (PhysicsAnalysis* nominal);
  virtual void ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar = true);
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Gets the relative variations from the MC closure study as a function of S/B ratio.
////////////////////////////////////////////////////////////////////////////////////////////////
void NonClosureVariation :: GetRelativeVariations (PhysicsAnalysis* nominal) {

  nominal->CalculateSigToBkg ();
 
  TDirectory* _gDirectory = gDirectory;

  TFile* f_pTch = new TFile (Form ("%s/ClosureSystematic/closureSystematic_pTch.root", rootPath.Data ()), "read");
  TF1* f_closure_sigToBkg_pTch = (TF1*) f_pTch->Get ("f_closure_sigToBkg_pTch");
  TFile* f_xhZ = new TFile (Form ("%s/ClosureSystematic/closureSystematic_xhZ.root", rootPath.Data ()), "read");
  TF1* f_closure_sigToBkg_xhZ = (TF1*) f_xhZ->Get ("f_closure_sigToBkg_xhZ");

  _gDirectory->cd ();

  TH1D* h = nullptr, *s2b = nullptr;

  cout << "Calculating non-closure variations on yields." << endl;

  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        h = (TH1D*) nominal->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_nonClosureVarPt_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        relVarPt[iSpc][iPtZ][iCent] = h;
        s2b = nominal->h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent];

        for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
          assert (h->GetBinCenter (iX) == s2b->GetBinCenter (iX));
          float rerr = 1;
          if (iCent == 0) rerr = 1.02; // constant 2% error in pp
          else if (!isnan (s2b->GetBinContent (iX)) && s2b->GetBinContent (iX) > 0)
            rerr = f_closure_sigToBkg_pTch->Eval (s2b->GetBinContent (iX));
          h->SetBinContent (iX, 1.+0.5*fabs(rerr-1.));
          h->SetBinError (iX, 0);
        }
        //cout << relVarPt[iSpc][iPtZ][iCent]->GetNbinsX () << endl;

        h = (TH1D*) nominal->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_nonClosureVarX_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        relVarX[iSpc][iPtZ][iCent] = h;
        s2b = nominal->h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent];

        for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
          assert (h->GetBinCenter (iX) == s2b->GetBinCenter (iX));
          float rerr = 1;
          if (iCent == 0) rerr = 1.02; // constant 2% error in pp
          else if (!isnan (s2b->GetBinContent (iX)) && s2b->GetBinContent (iX) > 0)
            rerr = f_closure_sigToBkg_xhZ->Eval (s2b->GetBinContent (iX));
          h->SetBinContent (iX, 1.+0.5*fabs(rerr-1.));
          h->SetBinError (iX, 0);
        }
        //cout << relVarX[iSpc][iPtZ][iCent]->GetNbinsX () << endl;
      } // end loop over iCent
    } // end loop over iSpc

    {
      const short iSpc = 2;
      for (short iCent : {0, numCentBins}) {

        TH1D* h_s2b = nominal->h_trk_dphi_sig_to_bkg[iSpc][iPtZ][maxNPtchBins][iCent];
        TH1D* h_sub = nominal->h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins][iCent];
        TH1D* h_bkg = nominal->bkg->h_trk_dphi[iSpc][iPtZ][maxNPtchBins][iCent];
        relVardPhi[iPtZ][0][iCent] = (TH1D*) h_sub->Clone (Form ("h_nonClosureRelVardPhi_comb_gt4_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()));
        absVardPhi[iPtZ][0][iCent] = (TH1D*) h_sub->Clone (Form ("h_nonClosureAbsVardPhi_comb_gt4_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()));

        float avg_rel_err = 1;
        float avg_abs_err = 0;
        float avg_s2b = 0;
        float avg_sub = 0;
        float avg_bkg = 0;
        short num_points = 0;
        for (int iX = 1; iX <= h_sub->GetNbinsX (); iX++) {
          if (h_sub->GetBinCenter (iX) > 3.*pi/4.)
            continue; // only average < 3pi/4 points
          avg_sub += h_sub->GetBinContent (iX);
          avg_bkg += h_bkg->GetBinContent (iX);
          num_points++;
        }
        assert (num_points > 0 && avg_bkg != 0);
        avg_s2b = avg_sub / avg_bkg;
        assert (avg_s2b > 0);

        avg_sub = avg_sub/num_points;
        avg_bkg = avg_bkg/num_points;

        if (iCent == 0) avg_rel_err = 1.02; // constant 2% error in pp
        else avg_rel_err = 1.+0.5*fabs (f_closure_sigToBkg_pTch->Eval (avg_s2b)-1.);
        avg_abs_err = (avg_rel_err-1.) * fabs (avg_sub);

        for (int iX = 1; iX <= relVardPhi[iPtZ][0][iCent]->GetNbinsX (); iX++) {
          if (h_sub->GetBinCenter (iX) < 3.*pi/4.) {
            relVardPhi[iPtZ][0][iCent]->SetBinContent (iX, avg_rel_err);
            relVardPhi[iPtZ][0][iCent]->SetBinError (iX, 0);
            absVardPhi[iPtZ][0][iCent]->SetBinContent (iX, avg_abs_err);
            absVardPhi[iPtZ][0][iCent]->SetBinError (iX, 0);
          }
          else {
            const float rel_err = (iCent == 0 ? 1.02 : 1.+0.5*fabs (f_closure_sigToBkg_pTch->Eval (h_s2b->GetBinContent (iX))-1.));
            relVardPhi[iPtZ][0][iCent]->SetBinContent (iX, rel_err);
            relVardPhi[iPtZ][0][iCent]->SetBinError (iX, 0);
            absVardPhi[iPtZ][0][iCent]->SetBinContent (iX, (rel_err-1.) * h_sub->GetBinContent (iX));
            absVardPhi[iPtZ][0][iCent]->SetBinError (iX, 0);
          }
        }


        h_s2b = nominal->h_trk_dphi_sig_to_bkg[iSpc][iPtZ][maxNPtchBins+1][iCent];
        h_sub = nominal->h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins+1][iCent];
        h_bkg = nominal->bkg->h_trk_dphi[iSpc][iPtZ][maxNPtchBins+1][iCent];
        relVardPhi[iPtZ][1][iCent] = (TH1D*) h_sub->Clone (Form ("h_nonClosureRelVardPhi_comb_lt4_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()));
        absVardPhi[iPtZ][1][iCent] = (TH1D*) h_sub->Clone (Form ("h_nonClosureAbsVardPhi_comb_lt4_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()));

        avg_rel_err = 1;
        avg_abs_err = 0;
        avg_s2b = 0;
        avg_sub = 0;
        avg_bkg = 0;
        num_points = 0;
        for (int iX = 1; iX <= h_sub->GetNbinsX (); iX++) {
          if (h_sub->GetBinCenter (iX) > 3.*pi/4.)
            continue; // only average < 3pi/4 points
          avg_sub += h_sub->GetBinContent (iX);
          avg_bkg += h_bkg->GetBinContent (iX);
          num_points++;
        }
        assert (num_points > 0 && avg_bkg != 0);
        avg_s2b = avg_sub / avg_bkg;
        assert (avg_s2b > 0);

        avg_sub = avg_sub/num_points;
        avg_bkg = avg_bkg/num_points;

        if (iCent == 0) avg_rel_err = 1.02; // constant 2% error in pp
        else avg_rel_err = 1.+0.5*fabs (f_closure_sigToBkg_pTch->Eval (avg_s2b)-1.);
        avg_abs_err = (avg_rel_err-1.) * fabs (avg_sub);

        for (int iX = 1; iX <= relVardPhi[iPtZ][1][iCent]->GetNbinsX (); iX++) {
          if (h_sub->GetBinCenter (iX) < 3.*pi/4.) {
            relVardPhi[iPtZ][1][iCent]->SetBinContent (iX, avg_rel_err);
            relVardPhi[iPtZ][1][iCent]->SetBinError (iX, 0);
            absVardPhi[iPtZ][1][iCent]->SetBinContent (iX, avg_abs_err);
            absVardPhi[iPtZ][1][iCent]->SetBinError (iX, 0);
          }
          else {
            const float rel_err = (iCent == 0 ? 1.02 : 1.+0.5*fabs (f_closure_sigToBkg_pTch->Eval (h_s2b->GetBinContent (iX))-1.));
            relVardPhi[iPtZ][1][iCent]->SetBinContent (iX, rel_err);
            relVardPhi[iPtZ][1][iCent]->SetBinError (iX, 0);
            absVardPhi[iPtZ][1][iCent]->SetBinContent (iX, (rel_err-1.) * h_sub->GetBinContent (iX));
            absVardPhi[iPtZ][1][iCent]->SetBinError (iX, 0);
          }
        }
      } // end loop over iCent
    } // end iSpc = 2 scope
  } // end loop over iPtZ

  f_pTch->Close ();
  f_xhZ->Close ();
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void NonClosureVariation :: ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar) {
  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        if (upVar) {
          if (iCent == 0) {
            a->h_trk_pt_ptz[iSpc][iPtZ][iCent]->Multiply (relVarPt[iSpc][iPtZ][iCent]);
            //TH1D* h = a->h_trk_pt_ptz[iSpc][iPtZ][iCent];
            //for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            //  assert (h->GetBinCenter (iX) == relVarPt[iSpc][iPtZ][iCent]->GetBinCenter (iX));
            //  h->SetBinContent (iX, h->GetBinContent (iX) * relVarPt[iSpc][iPtZ][iCent]->GetBinContent (iX));
            //}

            a->h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Multiply (relVarX[iSpc][iPtZ][iCent]);
            //h = a->h_trk_xhz_ptz[iSpc][iPtZ][iCent];
            //for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            //  assert (h->GetBinCenter (iX) == relVarX[iSpc][iPtZ][iCent]->GetBinCenter (iX));
            //  h->SetBinContent (iX, h->GetBinContent (iX) * relVarX[iSpc][iPtZ][iCent]->GetBinContent (iX));
            //}
          }
          //a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]->Multiply (relVarPt[iSpc][iPtZ][iCent]);
          TH1D* h = a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent];
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            assert (h->GetBinCenter (iX) == relVarPt[iSpc][iPtZ][iCent]->GetBinCenter (iX));
            h->SetBinContent (iX, h->GetBinContent (iX) * relVarPt[iSpc][iPtZ][iCent]->GetBinContent (iX));
          }

          //a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]->Multiply (relVarX[iSpc][iPtZ][iCent]);
          h = a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent];
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            assert (h->GetBinCenter (iX) == relVarX[iSpc][iPtZ][iCent]->GetBinCenter (iX));
            h->SetBinContent (iX, h->GetBinContent (iX) * relVarX[iSpc][iPtZ][iCent]->GetBinContent (iX));
          }
          //a->h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->Multiply (relVarPt[iSpc][iPtZ][iCent]);
          //a->h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->Multiply (relVarX[iSpc][iPtZ][iCent]);
        }
        else {
          if (iCent == 0) {
            //a->h_trk_pt_ptz[iSpc][iPtZ][iCent]->Divide (relVarPt[iSpc][iPtZ][iCent]);
            TH1D* h = a->h_trk_pt_ptz[iSpc][iPtZ][iCent];
            for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
              assert (h->GetBinCenter (iX) == relVarPt[iSpc][iPtZ][iCent]->GetBinCenter (iX));
              h->SetBinContent (iX, h->GetBinContent (iX) / relVarPt[iSpc][iPtZ][iCent]->GetBinContent (iX));
            }

            //a->h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Divide (relVarX[iSpc][iPtZ][iCent]);
            h = a->h_trk_xhz_ptz[iSpc][iPtZ][iCent];
            for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
              assert (h->GetBinCenter (iX) == relVarX[iSpc][iPtZ][iCent]->GetBinCenter (iX));
              h->SetBinContent (iX, h->GetBinContent (iX) / relVarX[iSpc][iPtZ][iCent]->GetBinContent (iX));
            }
          }

          a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]->Divide (relVarPt[iSpc][iPtZ][iCent]);
          TH1D* h = a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent];
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            assert (h->GetBinCenter (iX) == relVarPt[iSpc][iPtZ][iCent]->GetBinCenter (iX));
            h->SetBinContent (iX, h->GetBinContent (iX) / relVarPt[iSpc][iPtZ][iCent]->GetBinContent (iX));
          }

          //a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]->Divide (relVarX[iSpc][iPtZ][iCent]);
          h = a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent];
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            assert (h->GetBinCenter (iX) == relVarX[iSpc][iPtZ][iCent]->GetBinCenter (iX));
            h->SetBinContent (iX, h->GetBinContent (iX) / relVarX[iSpc][iPtZ][iCent]->GetBinContent (iX));
          }
          //a->h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->Divide (relVarPt[iSpc][iPtZ][iCent]);
          //a->h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->Divide (relVarX[iSpc][iPtZ][iCent]);
        }
      } // end loop over iCent
    } // end loop over iSpc

    {
      const short iSpc = 2;
      for (short iCent : {0, numCentBins}) {
        TH1D* h = a->h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins][iCent];
        for (int iX = 1; iX <= h->GetNbinsX (); iX++)
          h->SetBinContent (iX, h->GetBinContent (iX) + (upVar ? 1 : -1) * absVardPhi[iPtZ][0][iCent]->GetBinContent (iX));

        h = a->h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins+1][iCent];
        for (int iX = 1; iX <= h->GetNbinsX (); iX++)
          h->SetBinContent (iX, h->GetBinContent (iX) + (upVar ? 1 : -1) * absVardPhi[iPtZ][1][iCent]->GetBinContent (iX));
      } // end loop over iCent
    } // end iSpc = 2 scope
  } // end loop over iPtZ

  return;
}


#endif
