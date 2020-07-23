#ifndef __SystematicFits_h__
#define __SystematicFits_h__

#include "Params.h"
#include "PhysicsAnalysis.h"

#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>
#include <algorithm>
#include <map>

using namespace std;

class SystematicFits {

  private:
  string name;

  TF1*** relVarPt = Get2DArray <TF1*> (nPtZBins, numCentBins); // iPtZ, iCent
  TF1*** relVarX  = Get2DArray <TF1*> (nPtZBins, numCentBins); // iPtZ, iCent
  TF1**** relVardPhi = Get3DArray <TF1*> (nPtZBins, 2, 2); // iPtZ, iPtch (gt4 or lt4), iCent (pp or 0-30% PbPb)

  public:

  SystematicFits (const char* _name = "systematics") {
    name = _name;
  }

  ~SystematicFits () {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        SaferDelete (&relVarPt[iPtZ][iCent]);
        SaferDelete (&relVarX[iPtZ][iCent]);
      } // end loop over iCent
      SaferDelete (&relVardPhi[iPtZ][0][0]);
      SaferDelete (&relVardPhi[iPtZ][0][1]);
      SaferDelete (&relVardPhi[iPtZ][1][0]);
      SaferDelete (&relVardPhi[iPtZ][1][1]);
    } // end loop over iPtZ
    Delete2DArray (&relVarPt, nPtZBins, numCentBins);
    Delete2DArray (&relVarX, nPtZBins, numCentBins);
    Delete3DArray (&relVardPhi, nPtZBins, 2, 2);
  }

  string Name ()              { return name; }
  void SetName (string _name) { name = _name; }

  void GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* var);
  void ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar = true);

  //virtual void AddVariations (); // variations add linearly
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void SystematicFits :: GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* var) {
  TH1D* hn = nullptr, *hd = nullptr;
  TF1* f = nullptr;

  cout << "Calculating variation fits on total yields." << endl;

  //for (short iSpc = 2; iSpc < 3; iSpc++) {
  {
    const short iSpc = 2;
    const char* spc = "comb";
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        hn = (TH1D*) var->h_trk_pt_ptz[iSpc][iPtZ][iCent]->Clone ("temp");
        hd = nominal->h_trk_pt_ptz[iSpc][iPtZ][iCent];
        assert (hn->GetNbinsX () == hd->GetNbinsX ());
        for (int iX = 1; iX <= hn->GetNbinsX (); iX++) {
          const float yn = hn->GetBinContent (iX);
          const float yne = hn->GetBinError (iX);
          const float yd = hd->GetBinContent (iX);
          const float yde = hd->GetBinError (iX);
          hn->SetBinContent (iX, yn/yd);
          hn->SetBinError (iX, fabs(yn/yd) * sqrt (fabs (pow (yne/yn, 2) + pow (yde/yd, 2) - 2.*yne*yne/(yn*yd))));
        }
        f = new TF1 (Form ("f_relVarPt_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()), "[0]+[1]*log(x)+[2]*log(x)*log(x)", pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
        f->SetParameter (0, 1);
        f->SetParameter (1, 0);
        f->SetParameter (2, 0);
        hn->Fit (f, "RN0Q");
        delete hn;
        relVarPt[iPtZ][iCent] = f;

        hn = (TH1D*) var->h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Clone ("temp");
        hd = nominal->h_trk_xhz_ptz[iSpc][iPtZ][iCent];
        assert (hn->GetNbinsX () == hd->GetNbinsX ());
        for (int iX = 1; iX <= hn->GetNbinsX (); iX++) {
          const float yn = hn->GetBinContent (iX);
          const float yne = hn->GetBinError (iX);
          const float yd = hd->GetBinContent (iX);
          const float yde = hd->GetBinError (iX);
          hn->SetBinContent (iX, yn/yd);
          hn->SetBinError (iX, fabs(yn/yd) * sqrt (fabs (pow (yne/yn, 2) + pow (yde/yd, 2) - 2.*yne*yne/(yn*yd))));
        }
        f = new TF1 (Form ("f_relVarX_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()), "[0]+[1]*log(x)+[2]*log(x)*log(x)", xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
        f->SetParameter (0, 1);
        f->SetParameter (1, 0);
        f->SetParameter (2, 0);
        hn->Fit (f, "RN0Q");
        delete hn;
        relVarX[iPtZ][iCent] = f;
      } // end loop over iCent


      for (short iCent : {0, 1}) {
        for (short iPtch : {0, 1}) {

          hn = (TH1D*) var->h_trk_dphi[iSpc][iPtZ][maxNPtchBins+iPtch][numCentBins*iCent]->Clone ("temp");
          hd = nominal->h_trk_dphi[iSpc][iPtZ][maxNPtchBins+iPtch][numCentBins*iCent];
          assert (hn->GetNbinsX () == hd->GetNbinsX ());
          for (int iX = 1; iX <= hn->GetNbinsX (); iX++) {
            const float yn = hn->GetBinContent (iX);
            const float yne = hn->GetBinError (iX);
            const float yd = hd->GetBinContent (iX);
            const float yde = hd->GetBinError (iX);
            hn->SetBinContent (iX, yn/yd);
            hn->SetBinError (iX, fabs(yn/yd) * sqrt (fabs (pow (yne/yn, 2) + pow (yde/yd, 2) - 2.*yne*yne/(yn*yd))));
          }
          f = new TF1 (Form ("f_relVardPhi_%s_iPtZ%i_iCent%i_%s", (iPtch == 0 ? "gt4" : "lt4"), iPtZ, iCent, name.c_str ()), "[0]"/*+[1]*sin(x-[2])+[3]*sin(2*(x-[4]))+[5]*sin(3*(x-[6]))"*/, 0, pi);
          f->SetParameter (0, 1);
          //f->SetParameter (1, 0);
          //f->SetParameter (2, 0);
          //f->SetParameter (3, 0);
          //f->SetParameter (4, 0);
          //f->SetParameter (5, 0);
          //f->SetParameter (6, 0);
          hn->Fit (f, "RN0Q");
          delete hn;
          relVardPhi[iPtZ][iPtch][iCent] = f;

        } // end loop over iPtch
      } // end loop over iCent

    } // end loop over iPtZ
  } // end iSpc == 2 scope
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void SystematicFits :: ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar) {
  {
    const short iSpc = 2;
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

      for (short iCent = 0; iCent < numCentBins; iCent++) {

        for (TH1D* h : {a->h_trk_pt_ptz[iSpc][iPtZ][iCent], a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]}) {
          if (!h) continue;
          TF1* f = relVarPt[iPtZ][iCent];
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            float factor = 1;
            if (upVar)  factor = f->Eval (h->GetBinCenter (iX));
            else        factor = 2 - f->Eval (h->GetBinCenter (iX));
            h->SetBinContent (iX, h->GetBinContent (iX) * factor);
            h->SetBinError (iX, h->GetBinError (iX) * factor);
          }
        }

        for (TH1D* h : {a->h_trk_xhz_ptz[iSpc][iPtZ][iCent], a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]}) {
          if (!h) continue;
          TF1* f = relVarX[iPtZ][iCent];
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            float factor = 1;
            if (upVar)  factor = f->Eval (h->GetBinCenter (iX));
            else        factor = 2 - f->Eval (h->GetBinCenter (iX));
            h->SetBinContent (iX, h->GetBinContent (iX) * factor);
            h->SetBinError (iX, h->GetBinError (iX) * factor);
          }
        }
      } // end loop over iCent

      for (short iPtch : {0, 1}) {
        for (short iCent : {0, 1}) {

          TH1D* h = a->h_trk_dphi[iSpc][iPtZ][maxNPtchBins+iPtch][numCentBins*iCent];
          if (!h) continue;
          TF1* f = relVardPhi[iPtZ][iPtch][iCent];
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            float factor = 1;
            if (upVar)  factor = f->Eval (h->GetBinCenter (iX));
            else        factor = 2 - f->Eval (h->GetBinCenter (iX));
            h->SetBinContent (iX, h->GetBinContent (iX) * factor);
            h->SetBinError (iX, h->GetBinError (iX) * factor);
          }


          h = a->h_trk_dphi_sub[iSpc][iPtZ][maxNPtchBins+iPtch][numCentBins*iCent];
          if (!h) continue;

          f = relVardPhi[iPtZ][iPtch][iCent];
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            float factor = 1;
            if (upVar)  factor = f->Eval (h->GetBinCenter (iX));
            else        factor = 2 - f->Eval (h->GetBinCenter (iX));

            h->SetBinContent (iX, h->GetBinContent (iX) * factor);
            h->SetBinError (iX, h->GetBinError (iX) * factor);
          }

        } // end loop over iCent
      } // end loop over iPtch
    } // end loop over iPtZ
  } // end iSpc == 2 scope

  return;
}


#endif
