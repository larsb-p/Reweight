//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightXSecMEC

\brief    Model-independent tweak dials for reweighting MEC events

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\edits    Lars Bathe-Peters <lars.bathe-peters@physics.ox.ac.uk>
          University of Oxford

\created  Sep 11, 2019

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_MEC_H_
#define _G_REWEIGHT_NU_XSEC_MEC_H_

#include <map>
#include <memory>
#include <string>
#include <gsl/gsl_sf_legendre.h> // Include the GSL header for Legendre polynomials

#include <TGraph.h>

// GENIE includes
#include "Framework/Interaction/InteractionType.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "RwCalculators/GReWeightModel.h"
#include "RwFramework/GSyst.h"

namespace genie {

namespace rew   {

 class GReWeightXSecMEC : public GReWeightModel {
 public:

   GReWeightXSecMEC();
   GReWeightXSecMEC(std::string model, std::string type);
  ~GReWeightXSecMEC();

   // implement the GReWeightI interface
   bool   AppliesTo      (const EventRecord & event) const;
   bool   IsHandled      (GSyst_t syst) const;
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord& event);

 private:

   void Init(void);
   double CalcWeightNorm(const EventRecord& event);
   double CalcWeightAngularDist(const EventRecord& event);
   double CalcWeightAngularDistLegendre(const EventRecord& event);
   double CalcWeightPNDelta(const EventRecord& event);
   double CalcWeightXSecShape(const EventRecord& event);
   double CalcWeightXSecShape_Empirical(const EventRecord& event);
   double CalcWeightXSecShape_Martini(const EventRecord& event);
   double CalcWeightEnergyDependence(const EventRecord& event);

   /// Helper function for CalcWeightXSecShape
   double GetXSecIntegral(const XSecAlgorithmI* xsec_alg,
   const Interaction* interaction);

   /// Get total xsec, searching all loaded tunes for splines before
   /// falling back to numerical integration
   double GetXSecFromSplineOrIntegral(const XSecAlgorithmI* xsec_alg,
   const Interaction* interaction);

   /// Helper function for CalcWeightEnergyDependence
   double CalcWeight2p2hEnergyDependence(const EventRecord& event);

   /// Build energy-dependence ratio graphs from all three MEC models
   void BuildEnergyDepRatioGraphs(const EventRecord& event);

   /// Helper function for CalcWeightAngularDistLegendre. Returns
   /// max(0, 1 + k * F(theta)), where F(theta) = sum_l coeff_l * P_l(costheta)
   /// is the raw Legendre-polynomial shape and k is the per-throw
   /// amplitude scale. Since every P_l (l >= 1) integrates to zero over the
   /// isotropic base distribution in costheta, <F(theta)> = 0 and therefore
   /// <weight> = 1 for any k, preserving the overall normalization. The
   /// max(0, ...) is a hard backstop against negative weights in case the
   /// discrete theta-scan used to pick k (see ComputeLegendreAmplitudeK)
   /// slightly underestimates the true minimum of F(theta)
   double CalcWeightDecayAngMECLegendre(double theta_rad, double twk_dial, double twk_dial2, double twk_dial3, double twk_dial4, double twk_dial5, double twk_dial6, double k);

   /// Given the per-order Legendre coefficients (twk_dial...twk_dial6) for
   /// one systematic throw, find the largest amplitude k for which
   /// w(theta) = 1 + k*F(theta) stays non-negative over the full angular
   /// domain [0, pi]
   double ComputeLegendreAmplitudeK(double twk_dial, double twk_dial2, double twk_dial3, double twk_dial4, double twk_dial5, double twk_dial6);

   /// Simple struct containing tweak dial information for the
   /// normalization of one MEC interaction type (CC, NC, EM)
   struct NormMapEntry {
     NormMapEntry() : fTwkDial(0.), fNormDef(0.), fNormCurr(0.) {}
     NormMapEntry(double twk_dial, double def, double curr)
       : fTwkDial( twk_dial ), fNormDef( def ), fNormCurr( curr ) {}
     double fTwkDial;
     double fNormDef;
     double fNormCurr;
   };

   /// Map linking MEC interaction modes to normalization tweak dial values
   std::map<InteractionType_t, NormMapEntry> fNormMap;

   /// Lookup table linking MEC GSyst_t tweak dial enum labels to
   /// interaction modes
   static std::map<GSyst_t, InteractionType_t> fGSystToIntTypeMap;

   /// Tweak dial value for adjusting the nucleon cluster decay angular
   /// distribution
   double fDecayAngTwkDial;

   /// Another tweak dial value for adjusting the nucleon cluster decay
   /// angular distribution
   double fDecayAng2TwkDial;

   /// Tweak dial value for adjusting the nucleon cluster decay angular
   /// distribution for DecayAngMECLegendre parameter
   double fDecayAngLegendreTwkDial;
   double fDecayAngLegendre2TwkDial;
   double fDecayAngLegendre3TwkDial;
   double fDecayAngLegendre4TwkDial;
   double fDecayAngLegendre5TwkDial;
   double fDecayAngLegendre6TwkDial;

   /// Per-throw amplitude scale k applied to the DecayAngMECLegendre shape.
   /// Cached (along with the coefficient values it was computed from) so
   /// that it is only recomputed when the coefficients change, rather than
   /// on every event
   double fDecayAngLegendreK;
   double fDecayAngLegendreKCachedDials[6];
   bool   fDecayAngLegendreKCacheValid;

   /// Tweak dial value for adjusting the fraction of CC events that
   /// involve an initial pn pair
   double fFracPN_CCTwkDial;

   /// Tweak dial value for adjusting the fraction of CC events that
   /// involve an internal Delta
   double fFracDelta_CCTwkDial;

   /// CCMEC cross section model used to generate the events (untweaked)
   XSecAlgorithmI* fXSecAlgCCDef;

   /// Alternate CCMEC cross section models
   XSecAlgorithmI* fXSecAlgCCAlt_Nieves;
   XSecAlgorithmI* fXSecAlgCCAlt_SuSAv2;
   XSecAlgorithmI* fXSecAlgCCAlt_Empirical;
   XSecAlgorithmI* fXSecAlgCCAlt_Martini;

   /// Integrator used by the CalcWeightXSecShape function
   const XSecIntegratorI* fXSecIntegrator;

   /// Tweak dials that interpolates the shape of the CCMEC differential 
   // cross section between models
   double fCCXSecShapeTwkDial;
   double fCCXSecShapeEmpiricalTwkDial;
   double fCCXSecShapeMartiniTwkDial;

   /// Tweak dial value for adjusting the energy dependence of the CCMEC
   //cross section
   double fEnergyDependenceTwkDial;

   /// Energy-dependence ratio graphs (upper/lower uncertainty envelopes)
   bool fEnergyDepRatioInitialized;
   std::unique_ptr<TGraph> fEnergyDepUpperGraph;
   std::unique_ptr<TGraph> fEnergyDepLowerGraph;
};

} // rew   namespace
} // genie namespace

#endif // _G_REWEIGHT_NU_XSEC_MEC_H_
