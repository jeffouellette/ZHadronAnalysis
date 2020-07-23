#ifndef __ReweightingVariation_h__
#define __ReweightingVariation_h__

#include "Params.h"
#include "PhysicsAnalysis.h"

#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;

class ReweightingVariation {

  protected:
  string name;

  public:
  TH1D**** relVarPt     = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D**** relVarX      = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D***** relVardPhi  = Get4DArray <TH1D*> (3, nPtZBins, maxNPtchBins+2, numCentBins); // iSpc, iPtZ, iCent
  //TH1D***** absVardPhi  = Get4DArray <TH1D*> (3, nPtZBins, maxNPtchBins+2, numCentBins+1); // iSpc, iPtZ, iPtch, iCent

  //ReweightingVariation (PhysicsAnalysis* nom, const char* _name = "systematics", const char* _desc = "systematic") : Systematic (nom, _name, _desc) { }
  ReweightingVariation (const char* _name = "reweightingVariation") {
    name = _name;
  }

  ~ReweightingVariation () {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          SaferDelete (&relVarPt[iSpc][iPtZ][iCent]);
          SaferDelete (&relVarX[iSpc][iPtZ][iCent]);
          for (short iPtch = maxNPtchBins; iPtch < maxNPtchBins+2; iPtch++) {
            SaferDelete (&relVardPhi[iSpc][iPtZ][iPtch][iCent]);
            //SaferDelete (&absVardPhi[iSpc][iPtZ][iPtch][iCent]);
          }
        } // end loop over iCent
      } // end loop over iPtZ
    } // end loop over iSpc
    Delete3DArray (&relVarPt, 3, nPtZBins, numCentBins);
    Delete3DArray (&relVarX, 3, nPtZBins, numCentBins);
    Delete4DArray (&relVardPhi, 3, nPtZBins, maxNPtchBins+2, numCentBins);
    //Delete4DArray (&absVardPhi, 3, nPtZBins, maxNPtchBins+2, numCentBins);
  }

  string Name ()              { return name; }
  void SetName (string _name) { name = _name; }

  void ScaleRelativeVariations (const float scaleFactor = 1.);

  virtual void GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* var);
  virtual void ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar = true);
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Rescales the relative variation weights by some constant factor
////////////////////////////////////////////////////////////////////////////////////////////////
void ReweightingVariation :: ScaleRelativeVariations (const float scaleFactor) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
      for (int iPhi = 0; iPhi <= numPhiBins; iPhi++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          for (TH1D* h : {relVarPt[iSpc][iPtZ][iCent], relVarX[iSpc][iPtZ][iCent]}) {
            if (!h)
              continue;
            for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
              float relVar = h->GetBinContent (ix);
              relVar = (relVar-1.)*scaleFactor + 1.;
              h->SetBinContent (ix, relVar);
            }
          } // end loop over h
        } // end loop over iCent
      } // end loop over iPhi
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void ReweightingVariation :: GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* var) {
  TH1D* h_nom = nullptr, *h_var = nullptr;

  cout << "Calculating reweighting variation on total yields." << endl;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        h_nom = (TH1D*) nominal->h_trk_pt_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_relVarPt_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        relVarPt[iSpc][iPtZ][iCent] = h_nom;
        h_nom->Divide (var->h_trk_pt_ptz[iSpc][iPtZ][iCent]);

        h_nom = (TH1D*) nominal->h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_relVarX_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        relVarX[iSpc][iPtZ][iCent] = h_nom;
        h_nom->Divide (var->h_trk_xhz_ptz[iSpc][iPtZ][iCent]);


        for (short iPtch = maxNPtchBins; iPtch < maxNPtchBins+2; iPtch++) {
          h_nom = (TH1D*) nominal->h_trk_dphi[iSpc][iPtZ][iPtch][iCent];
          h_var = (TH1D*) var->h_trk_dphi[iSpc][iPtZ][iPtch][iCent];
          if(!h_nom || !h_var) continue;

          float nom_avg_y = 0.;
          float var_avg_y = 0.;
          int num_points = 0;
          for (int iX = 1; iX <= h_nom->GetNbinsX (); iX++) {
            if (h_nom->GetBinCenter (iX) > 3.*pi/4.) continue;
            nom_avg_y += h_nom->GetBinContent (iX);
            var_avg_y += h_var->GetBinContent (iX);
            num_points++;
          }
          assert (num_points > 0);
          nom_avg_y = nom_avg_y / num_points;
          var_avg_y = var_avg_y / num_points;

          //cout << "iSpc = " << iSpc << ", iCent = " << iCent << ", iPtch = " << iPtch << "ivar_avg_y var_avg_y << ", " << nom_avg_y << endl;

          //absVardPhi[iSpc][iPtZ][iPtch][iCent] = (TH1D*) h_nom->Clone (Form ("h_absVardPhi_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ()));
          relVardPhi[iSpc][iPtZ][iPtch][iCent] = (TH1D*) h_nom->Clone (Form ("h_relVardPhi_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ()));

          for (int iX = 1; iX <= h_nom->GetNbinsX (); iX++) {
            if (h_nom->GetBinCenter (iX) > 3.*pi/4.) {
              assert (h_var->GetBinContent (iX) != 0);
              //absVardPhi[iSpc][iPtZ][iPtch][iCent]->SetBinContent (iX, fabs (((h_nom->GetBinContent (iX) / h_var->GetBinContent (iX)) - 1.) * (h_nom->GetBinContent (iX))));
              relVardPhi[iSpc][iPtZ][iPtch][iCent]->SetBinContent (iX, fabs (h_nom->GetBinContent (iX) / h_var->GetBinContent (iX)));
            }
            else {
              //absVardPhi[iSpc][iPtZ][iPtch][iCent]->SetBinContent (iX, fabs (nom_avg_y - var_avg_y));
              relVardPhi[iSpc][iPtZ][iPtch][iCent]->SetBinContent (iX, fabs (nom_avg_y/var_avg_y));
            }
          } // end loop over iX
        } // end loop over iPtch
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void ReweightingVariation :: ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        if (upVar) {
          a->h_trk_pt_ptz[iSpc][iPtZ][iCent]->Multiply (relVarPt[iSpc][iPtZ][iCent]);
          a->h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Multiply (relVarX[iSpc][iPtZ][iCent]);
        }
        else {
          a->h_trk_pt_ptz[iSpc][iPtZ][iCent]->Divide (relVarPt[iSpc][iPtZ][iCent]);
          a->h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Divide (relVarX[iSpc][iPtZ][iCent]);
        }

        for (short iPtch = maxNPtchBins; iPtch < maxNPtchBins+2; iPtch++) {
          TH1D* h = a->h_trk_dphi[iSpc][iPtZ][iPtch][iCent];
          if (!h) continue;

          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            h->SetBinContent (iX, h->GetBinContent (iX) * pow (relVardPhi[iSpc][iPtZ][iPtch][iCent]->GetBinContent (iX), upVar ? 1. : -1.));
            //h->SetBinContent (iX, h->GetBinContent (iX) + (upVar ? 1. : -1.) * absVardPhi[iSpc][iPtZ][iPtch][iCent]->GetBinContent (iX));
          } // end loop over iX
        } // end loop over iPtch
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc

  return;
}


#endif
