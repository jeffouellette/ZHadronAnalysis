#ifndef __TrackIDSystematic_h__
#define __TrackIDSystematic_h__

#include "Params.h"
#include "PhysicsAnalysis.h"

#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>
#include <algorithm>
#include <map>

using namespace std;

class TrackIDSystematic : public Systematic {

  public:

  float**** relVar = Get4DArray <float> (3, nPtZBins, numPhiBins+1, numCentBins); // iSpc, iPtZ, iPhi (or integrated at numPhiBins), iCent
  TrackIDSystematic (FullAnalysis* nom, const char* _name = "systematics", const char* _desc = "systematic") : Systematic (nom, _name, _desc) { }


  void GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* var);
  void ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar = true);

  //virtual void AddVariations (); // variations add linearly
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void TrackIDSystematic :: GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* var) {
  TH1D* hnom = nullptr, *hvar = nullptr, *hrat = nullptr, *eff1 = nullptr, *eff2 = nullptr, *effrat = nullptr, *purrat = nullptr;

  cout << "Calculating track ID variation on total yields." << endl;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

      //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
      //  for (short iCent = 0; iCent < numCentBins; iCent++) {
      //    //cout << "Getting relative variation for iSpc=" << iSpc << ", iPtZ=" << iPtZ << ", iPhi=" << iPhi << ", iCent=" << iCent << endl;

      //    hnom = nominal->h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent];
      //    hvar = var->h_trk_pt_dphi_raw[iSpc][iPtZ][iPhi][iCent];

      //    hrat = (TH1D*) hvar->Clone ("rattemp");
      //    hrat->Divide (hnom);
      //    for (int ix = 1; ix <= hrat->GetNbinsX (); ix++)
      //      hrat->SetBinError (ix, sqrt ((hrat->GetBinContent (ix)) * (1-hrat->GetBinContent (ix)) / hnom->GetBinContent (ix)));

      //    nominal->LoadTrackingEfficiencies (true);
      //    var->LoadTrackingEfficiencies (true);
      //    nominal->LoadTrackingPurities (true);
      //    var->LoadTrackingPurities (true);

      //    eff1 = (TH1D*) nominal->h2_num_trk_effs[iCent]->ProjectionY ("1");
      //    effrat = (TH1D*) nominal->h2_den_trk_effs[iCent]->ProjectionY ("2");
      //    eff1->Divide (effrat);
      //    delete effrat;
      //    eff2 = (TH1D*) var->h2_num_trk_effs[iCent]->ProjectionY ("3");
      //    effrat = (TH1D*) var->h2_den_trk_effs[iCent]->ProjectionY ("4");
      //    eff2->Divide (effrat);
      //    delete effrat;

      //    effrat = (TH1D*) eff2->Clone ("efftemp");
      //    effrat->Divide (eff1);
      //    delete eff1, eff2;

      //    eff1 = (TH1D*) nominal->h2_num_trk_purs[iCent]->ProjectionY ("1");
      //    eff1->Divide ((TH1D*) nominal->h2_den_trk_purs[iCent]->ProjectionY ("2"));
      //    eff2 = (TH1D*) var->h2_num_trk_purs[iCent]->ProjectionY ("3");
      //    eff2->Divide ((TH1D*) var->h2_den_trk_purs[iCent]->ProjectionY ("4"));

      //    purrat = (TH1D*) eff2->Clone ("purtemp");
      //    purrat->Divide (eff1);
      //    delete eff1, eff2;

      //    const int bin1 = effrat->FindFixBin (hrat->GetBinCenter (1));
      //    const int bin2 = purrat->FindFixBin (hrat->GetBinCenter (1));
      //    for (int iy = 1; iy <= hrat->GetNbinsX (); iy++)
      //      hrat->SetBinContent (iy, hrat->GetBinContent (iy) * purrat->GetBinContent (iy+bin2-1) / effrat->GetBinContent (iy+bin1-1));

      //    TF1* fit = new TF1 ("fit", "[0]", pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]);
      //    fit->SetParameter (0, 1);
      //    hrat->Fit (fit, "RQN0");
      //    relVar[2][iPtZ][iPhi][iCent] = fit->GetParameter (0);
      //    relVar[0][iPtZ][iPhi][iCent] = fit->GetParameter (0);
      //    relVar[1][iPtZ][iPhi][iCent] = fit->GetParameter (0);
      //    delete hrat, effrat, purrat, fit;
      //  } // end loop over iCent
      //} // end loop over iPhi

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        hnom = (TH1D*) nominal->h_trk_pt_dphi_raw[iSpc][iPtZ][1][iCent]->Clone ("hnom");
        hvar = (TH1D*) var->h_trk_pt_dphi_raw[iSpc][iPtZ][1][iCent]->Clone ("hvar");
        hnom->Add (nominal->h_trk_pt_dphi_raw[iSpc][iPtZ][2][iCent]);
        hvar->Add (var->h_trk_pt_dphi_raw[iSpc][iPtZ][2][iCent]);

        hrat = (TH1D*) hvar->Clone ("rattemp");
        hrat->Divide (hnom);
        for (int ix = 1; ix <= hrat->GetNbinsX (); ix++)
          hrat->SetBinError (ix, sqrt ((hrat->GetBinContent (ix)) * (1-hrat->GetBinContent (ix)) / hnom->GetBinContent (ix)));

        nominal->LoadTrackingEfficiencies (true);
        var->LoadTrackingEfficiencies (true);
        nominal->LoadTrackingPurities (true);
        var->LoadTrackingPurities (true);

        eff1 = (TH1D*) nominal->h2_num_trk_effs[iCent]->ProjectionY ("1");
        effrat = (TH1D*) nominal->h2_den_trk_effs[iCent]->ProjectionY ("2");
        eff1->Divide (effrat);
        SaferDelete (&effrat);
        eff2 = (TH1D*) var->h2_num_trk_effs[iCent]->ProjectionY ("3");
        effrat = (TH1D*) var->h2_den_trk_effs[iCent]->ProjectionY ("4");
        eff2->Divide (effrat);
        SaferDelete (&effrat);

        effrat = (TH1D*) eff2->Clone ("efftemp");
        effrat->Divide (eff1);
        SaferDelete (&eff1);
        SaferDelete (&eff2);

        eff1 = (TH1D*) nominal->h2_num_trk_purs[iCent]->ProjectionY ("1");
        eff1->Divide ((TH1D*) nominal->h2_den_trk_purs[iCent]->ProjectionY ("2"));
        eff2 = (TH1D*) var->h2_num_trk_purs[iCent]->ProjectionY ("3");
        eff2->Divide ((TH1D*) var->h2_den_trk_purs[iCent]->ProjectionY ("4"));

        purrat = (TH1D*) eff2->Clone ("purtemp");
        purrat->Divide (eff1);
        SaferDelete (&eff1);
        SaferDelete (&eff2);

        const int bin1 = effrat->FindFixBin (hrat->GetBinCenter (1));
        const int bin2 = purrat->FindFixBin (hrat->GetBinCenter (1));
        for (int iy = 1; iy <= hrat->GetNbinsX (); iy++)
          hrat->SetBinContent (iy, hrat->GetBinContent (iy) * purrat->GetBinContent (iy+bin2-1) / effrat->GetBinContent (iy+bin1-1));

        TF1* fit = new TF1 ("fit", "[0]", pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
        fit->SetParameter (0, 1);
        hrat->Fit (fit, "RQN0");
        relVar[2][iPtZ][numPhiBins][iCent] = fit->GetParameter (0);
        relVar[0][iPtZ][numPhiBins][iCent] = fit->GetParameter (0);
        relVar[1][iPtZ][numPhiBins][iCent] = fit->GetParameter (0);

        SaferDelete (&hrat);
        SaferDelete (&effrat);
        SaferDelete (&purrat);
        SaferDelete (&fit);
        SaferDelete (&hnom);
        SaferDelete (&hvar);
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc
}



////////////////////////////////////////////////////////////////////////////////////////////////
// Only applies to pT dependent histograms, i.e. will need to subtract background still.
////////////////////////////////////////////////////////////////////////////////////////////////
void TrackIDSystematic :: ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {

      // Hadron yield systematics, signal & signal+bkg levels
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        //  h_trk_pt_dphi[iSpc][iPtZ][iPhi][iCent]->Scale (upVar ? relVar[iSpc][iPtZ][iPhi][iCent] : 1/relVar[iSpc][iPtZ][iPhi][iCent]);
        //  h_trk_xhz_dphi[iSpc][iPtZ][iPhi][iCent]->Scale (upVar ? relVar[iSpc][iPtZ][iPhi][iCent] : 1/relVar[iSpc][iPtZ][iPhi][iCent]);
        //} // end loop over iPhi

        a->h_trk_pt_ptz[iSpc][iPtZ][iCent]->Scale (upVar ? relVar[iSpc][iPtZ][numPhiBins][iCent] : 1/relVar[iSpc][iPtZ][numPhiBins][iCent]);
        a->h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Scale (upVar ? relVar[iSpc][iPtZ][numPhiBins][iCent] : 1/relVar[iSpc][iPtZ][numPhiBins][iCent]);

        for (short iPtch = maxNPtchBins; iPtch < maxNPtchBins+2; iPtch++) {
          TH1D* h = a->h_trk_dphi[iSpc][iPtZ][iPtch][iCent];
          if (!h) continue;

          float avg_y = 0; 
          short num_points = 0;
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            if (h->GetBinCenter (iX) > 3.*pi/4.) continue;
            avg_y += h->GetBinContent (iX);
            num_points++;
          }
          assert (num_points > 0);
          avg_y = avg_y / num_points;

          const float avg_abs_err = avg_y * fabs (relVar[iSpc][iPtZ][numPhiBins][iCent] - 1.);

          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            if (h->GetBinCenter (iX) < 3.*pi/4.) {
              h->SetBinContent (iX, h->GetBinContent (iX) + (upVar ? 1. : -1) * avg_abs_err);
              h->SetBinError (iX, h->GetBinError (iX) + (upVar ? 1. : -1) * avg_abs_err);
            }
            else {
              h->SetBinContent (iX, h->GetBinContent (iX) * pow (relVar[iSpc][iPtZ][numPhiBins][iCent], upVar ? 1. : -1.));
              h->SetBinError (iX, h->GetBinError (iX) * pow (relVar[iSpc][iPtZ][numPhiBins][iCent], upVar ? 1. : -1.));
            }
          } // end loop over iX
        } // end loop over iPtch
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc
} 


#endif
