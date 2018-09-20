#ifndef UNIGENQA_H
#define UNIGENQA_H 1

#include <iostream>
#include <fstream>
#include <vector>
#include <TMath.h>
#include <Rtypes.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TProfile.h>
#include <TChain.h>
#include <TString.h>

#include "URun.h"
#include "UEvent.h"
#include "UParticle.h"

using namespace std;

namespace qa {
/**
 * Axes definition
 */

enum EMomentumAxis {
  kPT = 0, kETA, kPHI, kYM, kEcm, kElab, kPcm, kPlab, kMcm, kMlab, kA, kZ, kMpdg, kMcm_Ecm, kMlab_Elab, kAxes
};

struct TMomentumAxis {
  Int_t id;
  TString name;
  TString displayName;
  Int_t nBins;
  Double_t min;
  Double_t max;
} gMomentumAxes[kAxes] = {
    /* traverse momentum */
    {.id = kPT, .name = "Pt", .displayName="p_{T} (GeV/#it{c})", .nBins = 1000, .min=0.0, .max=5.},
    /* pseudorapidity */
    {.id = kETA, .name = "Eta", .displayName="\\eta", .nBins = 1000, .min=-10., .max=12.},
    /* azimuthal angle */
    {.id = kPHI, .name = "Phi", .displayName="\\varphi (rad)", .nBins = 1000, .min=-3.15, .max=3.15},
    /* rapidity in CoM frame */
    {.id = kYM, .name = "Y", .displayName="#it{y}", .nBins = 1000, .min=-5.0, .max=5.0},
    /* energy in CoM frame */
    {.id = kEcm, .name = "Ecm", .displayName="E_{CM} (GeV)", .nBins = 1000, .min=0., .max=1000.0},
    /* energy in lab frame */
    {.id = kElab, .name = "Elab", .displayName="E_{LS} (GeV)", .nBins = 1000, .min=0., .max=2400.0},
    /* full momentum in CoM frame */
    {.id = kPcm, .name = "Pcm", .displayName="p_{CM} (GeV/#it{c})", .nBins = 1000, .min=0., .max=500.},
    /* full momentum in lab frame */
    {.id = kPlab, .name = "Plab", .displayName="p_{LS} (GeV/#it{c})", .nBins = 1000, .min=0., .max=2400.},
    /* mass in CoM frame */
    {.id = kMcm, .name = "Mcm", .displayName="M_{CMS} (GeV/#it{c^{2}})", .nBins = 200, .min=0., .max=200},
    /* mass in lab frame */
    {.id = kMlab, .name = "Mlab", .displayName="M_{LS} (GeV/#it{c^{2}})", .nBins = 200, .min=0., .max=200},
    /* mass number */
    {.id = kA, .name = "A", .displayName="Mass number", .nBins = 200, .min=0., .max=200},
    /* number of protons */
    {.id = kZ, .name = "Z", .displayName="Charge number", .nBins = 82, .min=0., .max=82},
    /* mass from PDG code */
    {.id = kMpdg, .name = "Mpdg", .displayName="M_{PDG}", .nBins = 200, .min=0., .max=200},
    /* E/M (CMS) */
    {.id = kMcm_Ecm, .name = "Mcm_Ecm", .displayName="M/E (CMS)", .nBins = 100, .min=0., .max=1.},
    /* E/M (LS) */
    {.id = kMlab_Elab, .name = "Mlab_Elab", .displayName="M/E (LS)", .nBins = 100, .min=0., .max=1.},
};

/**
 * Particles definition
 */
enum EParticles {
  kALLSPECIES = 0,
  kLEPTONS,
  kFRAGMENTS,
  kPHOTON,
  kPROTON,
  kPROTONBAR,
  kPIPLUS,
  kPIMINUS,
  kKPLUS,
  kKMINUS,
  kLAMBDA,
  kLAMBDABAR,
  kParticles
};

const struct TParticle {
  Int_t id;
  Int_t pdg;
  Double_t mass;
  Int_t charge;
  std::string name;
  std::string displayName;
} gParticles[kParticles] = {
    {kALLSPECIES, 9999, 9999, 9999, "_all_species", "all species"},
    {kLEPTONS, 99999, 99999, 99999, "_leptons", "leptons"},
    {kFRAGMENTS, 999999, 999999, 999999, "_fragments", "fragments"},
    {kPHOTON, 22, 0., 0, "_photons", "photons"},
    {kPROTON, 2212, 0.938, 1, "_p", "p"},
    {kPROTONBAR, -2212, 0.938, -1, "_pbar", "\\bar{p}"},
    {kPIPLUS, 211, 0.138, 1, "_piplus", "\\pi^{+}"},
    {kPIMINUS, -211, 0.138, -1, "_piminus", "\\pi^{-}"},
    {kKPLUS, 321, 0.494, 1, "_kplus", "K^{+}"},
    {kKMINUS, -321, 0.494, -1, "_kminus", "K^{-}"},
    {kLAMBDA, 3122, 1.116, 1, "_lambda", "\\Lambda"},
    {kLAMBDABAR, -3122, 1.116, -1, "_lambda_bar", "\\bar{\\Lambda}"}
};

/**
 * PSD Groups definition
 */
enum EPSDGroups {
  kEFull = 0, kPSDAll, kPSD1, kPSD2, kPSD3, kPSDGroups
};
const struct TPSDGroup {
  Int_t id;
  std::string name;
  std::string displayName;
  std::vector<double> theta;
} gPSDGroups[kPSDGroups] = {
//				{.id = kEFull,  .name = "Efull", 	.displayName = "Full E", 		.theta = {0, TMath::Pi ()}},
    {.id = kEFull, .name = "Efull", .displayName = "Full E", .theta = {0, 3.15}},
    {.id = kPSDAll, .name = "PSDAll", .displayName = "Full PSD", .theta = {0, 0.074860}},
    {.id = kPSD1, .name = "PSD1", .displayName = "PSD1", .theta = {0., 0.024995}},
    {.id = kPSD2, .name = "PSD2", .displayName = "PSD2", .theta = {0.024995, 0.049958}},
    {.id = kPSD3, .name = "PSD3", .displayName = "PSD3", .theta = {0.049958, 0.074860}}
};

class UnigenQA {

 public:

  UnigenQA();

  ~UnigenQA();

  void Init(TString fileName, TString treeName);

  void Init(TChain *chain);

  void Init_Histograms();

  void FillEventInfo();

  void FillTracks();

  void Write_Histograms(const TString filePath);

  void SetReferenceChain(TChain *fReferenceChain) { UnigenQA::fReferenceChain = fReferenceChain; }

  void Run(Int_t nEvents = 1e9);

 private:

  vector<vector<int>> gPSDCorr = {
      {kPSD1, kPSD2},
      {kPSD2, kPSD3},
      {kPSD1, kPSD3}
  };

  vector<vector<int>> gTH2Axes = {
      {kETA, kPT},
      {kPHI, kPT},
      {kPHI, kETA},
      {kYM, kPT},
      {kPHI, kYM},
      {kPcm, kEcm},
      {kPlab, kElab},
      {kPlab, kPcm},
      {kElab, kEcm},
      {kPcm, kMcm_Ecm},
      {kPlab, kMlab_Elab},
      {kA, kEcm},
      {kA, kElab},
      {kA, kPcm},
      {kA, kPlab},
      {kMcm, kEcm},
      {kMlab, kElab},
      {kMcm, kMlab},
      {kMcm, kMpdg},
      {kMlab, kMpdg},
  };

  const TMomentumAxis
      gMultiplicity = {.id = kAxes, .name = "Mult", .displayName = "Multiplicity", .nBins = 250, .min = 0, .max = 250};

  Bool_t fDeleteChainOnExit{kFALSE};
  TChain *fChain;

  TChain *fReferenceChain{};
  UEvent *event_;
  double fPSDMax, fSnn{-999.}, fPcm{-999.}, fPlab{-999.}, fElab{-999.}, fEkin{-999.}, fBeta{-999.}, fA{-999.},
      fZ{-999.};
  vector<vector<int>> fPidGroups = {{0, 1099999999}, {37, 1099999999}};
  vector<TString> fPidGroupNames = {"all species", "hadrons"};
  double fPSDGroupEnergy[kPSDGroups][2];
  int fAmax, fZmax;

  /* PSD correlation histograms */
  TH1D *hPSDGroupEnergy[kPSDGroups][2];
  TH2D *hPSDMultCorr[kPSDGroups][2];
  TH2D *hBPSDCorr[kPSDGroups][2];
  TH2D *hPSDGroupsCorr[3][2];

  TH1D *hM;
  TH1D *hB;
  TH1D *hPsi;
  TH2D *hMBcorr;
  TH2D *h2BAmax;
  TH2D *h2BZmax;

  TH1D *hPdg;
  TH1D *hA;
  TH1D *hZ;
  TH2D *h2ZA;
  TH1D *hTrackMomentum[kAxes][kParticles];
  vector<TH2D **> hTrackMomentumCorr;
  TH2D *h2ElabA;
  TH2D *h2EcmA;
  TH2D *h2MlabA;
  TH2D *h2McmA;
  TH2D *h2PcmMcm_Ecm;
  TH2D *h2PlabMlab_Elab;

  // Flow
  TProfile *pVn_pT[2][kParticles];
  TProfile *pVn_Y[2][kParticles];

  // Yields
  TProfile *hYields[kParticles];
};

}
#endif
