#include <TLorentzVector.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include "UnigenQA.h"

using namespace TMath;
using namespace std;
using namespace qa;

UnigenQA::UnigenQA()
{
}

UnigenQA::~UnigenQA()
{
}

void UnigenQA::Init(TString filePath, TString treeName)
{
  fChain = new TChain(treeName);
  fChain->Add(filePath);
  event_ = new UEvent;
  fChain->SetBranchAddress("event", &event_);

  // Get run description
  TString fileName;

  if(!filePath.Contains("*"))
    fileName = filePath;
  else
  {
    TString dirName = filePath.Remove(filePath.Last('/') + 1, 10000);
    TSystemDirectory dir(dirName, dirName);
    TList *files = dir.GetListOfFiles();
    if(files)
    {
      TSystemFile *sysFile;
      TIter next(files);
      while((sysFile = (TSystemFile *)next()) && !fileName.EndsWith(".root"))
      {
        fileName = dirName + sysFile->GetName();
      }
    }
  }

  cout << "fileName = " << fileName << endl;

  TFile *file = new TFile(fileName, "open");
  URun *run = (URun *)file->Get("run");
  if(NULL == run)
  {
    cout << "No run description in input file." << endl;
  }
  run->Dump();
  if(TMath::Abs(run->GetPTarg()) > 0.001)
    cout << "Input data is in CM frame" << endl;
  double mProton = 0.938272029;
  fA = run->GetAProj() > run->GetATarg() ? run->GetAProj() : run->GetATarg();
  fZ = run->GetZProj() > run->GetZTarg() ? run->GetZProj() : run->GetZTarg();
  fSnn = run->GetNNSqrtS();
  fPcm = 0.5 * sqrt(fSnn * fSnn - 4. * mProton * mProton);
  fElab = 2. * fPcm * fPcm / mProton + mProton;
  fPlab = sqrt(fElab * fElab - mProton * mProton);
  fEkin = fSnn * fSnn / 1.87 - 1.87;
  fBeta = 2. * fPcm / fSnn;
  cout << "A = " << fA << endl;
  cout << "Z = " << fZ << endl;
  cout << "Snn = " << fSnn << " AGeV" << endl;
  cout << "Pcm = " << fPcm << " AGeV" << endl;
  cout << "Plab = " << fPlab << " AGeV" << endl;
  cout << "Elab = " << fElab << " AGeV" << endl;
  cout << "Ekin = " << fEkin << " AGeV" << endl;
  cout << "beta = " << fBeta << endl;
}

TChain *UnigenQA::MakeChain(TString filename, TString treename)
{
  TChain *chain = new TChain(treename);
  ifstream in;
  in.open(filename);
  string line;
  cout << "Adding files to chain:\n";
  while((in >> line).good())
  {
    if(!line.empty())
    {
      chain->AddFile(line.data());
      cout << line << endl;
    }
  }
  return chain;
}

void UnigenQA::Run(Long64_t nEvents)
{
  Long64_t nEntries = nEvents < fChain->GetEntries() ? nEvents : fChain->GetEntries();
  Long64_t outputStep = nEntries / 10;
  if(outputStep == 0)
    outputStep = 1;
  std::cout << "Entries = " << nEntries << std::endl;

  for(int i = 0; i < nEntries; i++)
  {
    if((i + 1) % outputStep == 0)
      std::cout << i + 1 << "/" << nEntries << "\r" << std::flush;
    fChain->GetEntry(i);
    FillTracks();     // !!! particle loop goes before the event loop (energy summ is calculated in the former)
    FillEventInfo();  // !!! particle loop goes before the event loop (energy summ is calculated in the former)
  }
}

void UnigenQA::Init_Histograms()
{
  gMomentumAxes[kEcm].max = round(fSnn * 2);
  gMomentumAxes[kPcm].max = round(fSnn * 2);
  gMomentumAxes[kElab].max = round(fElab * 3);
  gMomentumAxes[kPlab].max = round(fElab * 3);
  gMomentumAxes[kPzLab].max = round(fElab * 3);
  gMomentumAxes[kMcm].max = fA + 2;
  // gMomentumAxes[kMlab].max = fA;
  gMomentumAxes[kA].max = fA + 2;
  gMomentumAxes[kA].nBins = fA + 2;
  gMomentumAxes[kZ].max = fZ + 2;
  gMomentumAxes[kZ].nBins = fZ + 2;

  if(fReferenceChain == nullptr)
    fReferenceChain = fChain;
  else
  {
    cout << "Using reference chain..." << endl;
  }
  fPSDMax = round(fElab*fA/100)*110;
  Double_t fMmax = round(fReferenceChain->GetMaximum("fNpa")/100.)*110;

  cout << "fPSDMax = " << fPSDMax << endl;
  cout << "fMmax = " << fMmax << endl;

  TString name, title;

  hM = new TH1D("hM", "Track multiplicity;Multiplicity;nEvents", Nint(fMmax), 0, fMmax);
  hB = new TH1D("hB", "Impact parameter; B (fm);nEvents", 200, 0, 20);
  hPsi = new TH1D("hPsi", "#Psi_{RP};#Psi_{RP} (rad);Nevents", 100, -3.15, 3.15);

  hMBcorr = new TH2D("hMBcorr", "M : B;Multiplicity;B (fm)", Nint(fMmax), 0, fMmax, 200, 0, 20);
  h2BAmax = new TH2D("h2BAmax", "B : A_{max};B (fm);A_{max}", 200, 0, 20, 200, 0, 200);
  h2BZmax = new TH2D("h2BZmax", "B : Z_{max};B (fm);Z_{max}", 200, 0, 20, 82, 0, 82);

  /* PSD histogram initialization */
  fPSDGroupEnergy.resize(kPSDGroups, vector<double>(fPidGroups.size()));
  hPSDGroupEnergy.resize(kPSDGroups);
  hPSDGroupsCorr.resize(kPSDGroups);
  hPSDMultCorr.resize(kPSDGroups);
  hBPSDCorr.resize(kPSDGroups);
  Int_t nbins = 500;
  for(uint pidGroup = 0; pidGroup < fPidGroups.size(); pidGroup++)
  {
    for(uint psdCorr = 0; psdCorr < gPSDCorr.size(); psdCorr++)
    {
      auto group1 = gPSDGroups[gPSDCorr[psdCorr][0]];
      auto group2 = gPSDGroups[gPSDCorr[psdCorr][1]];
      
      name = "hPSDGroupsCorr_" + group1.name + "_" + group2.name + "_" + fPidGroupNames.at(pidGroup);
      hPSDGroupsCorr.at(psdCorr).push_back(new TH2D(name, name, nbins, 0, fPSDMax, nbins, 0, fPSDMax));
      hPSDGroupsCorr.at(psdCorr).at(pidGroup)->SetXTitle(group1.displayName.c_str());
      hPSDGroupsCorr.at(psdCorr).at(pidGroup)->SetYTitle(group2.displayName.c_str());
    }
    for(auto group : gPSDGroups)
    {
      hPSDGroupEnergy.resize(kPSDGroups);
      name = "hEnergy_" + group.name + "_" + fPidGroupNames.at(pidGroup);
      title = name + ";E (GeV);nEvents";
      hPSDGroupEnergy.at(group.id).push_back(new TH1D(name, title, nbins, 0, fPSDMax));
      hPSDMultCorr.resize(kPSDGroups);
      name = "hPSDMultCorr_" + group.name + "_" + fPidGroupNames.at(pidGroup);
      title = name + ";Multiplicity;E (GeV)";
      hPSDMultCorr.at(group.id).push_back(new TH2D(name, title, Nint(fMmax), 0, Nint(fMmax), nbins, 0, fPSDMax));
      hBPSDCorr.resize(kPSDGroups);
      name = "hBPSDCorr_" + group.name + "_" + fPidGroupNames.at(pidGroup);
      title = name + ";E (GeV);B (fm)";
      hBPSDCorr.at(group.id).push_back(new TH2D(name, title, nbins, 0, fPSDMax, 200, 0, 20));
    }
  }

  for(Int_t iPart = 0; iPart < kParticles; ++iPart)
  {
    auto particle = gParticles[iPart];

    for(Int_t iMom = 0; iMom < kAxes; ++iMom)
    {
      auto axis = gMomentumAxes[iMom];
      name = "hTrack" + axis.name + particle.name;
      title = axis.name + " (" + particle.displayName + ")";
      hTrackMomentum[iMom][iPart] = new TH1D(name, name, axis.nBins, axis.min, axis.max);
      hTrackMomentum[iMom][iPart]->SetXTitle(axis.displayName);
      hTrackMomentum[iMom][iPart]->SetYTitle("Counts");
    }

    for(uint iCorr = 0; iCorr < gTH2Axes.size(); ++iCorr)
    {
      vector<int> axes = gTH2Axes[iCorr];

      auto xAxis = gMomentumAxes[axes[0]];
      auto yAxis = gMomentumAxes[axes[1]];

      name = "hTrack" + xAxis.name + yAxis.name + particle.name;
      title = yAxis.name + " : " + xAxis.name + " (" + particle.displayName + ")";
      hTrackMomentumCorr.push_back(new TH2D *[kParticles]);
      hTrackMomentumCorr[iCorr][iPart] =
          new TH2D(name, title, xAxis.nBins, xAxis.min, xAxis.max, yAxis.nBins, yAxis.min, yAxis.max);
      hTrackMomentumCorr[iCorr][iPart]->SetXTitle(xAxis.displayName);
      hTrackMomentumCorr[iCorr][iPart]->SetYTitle(yAxis.displayName);
    }

    name = "h2Yield" + particle.name;
    title = "h2Yield_" + particle.displayName + ";B (fm);#it{y};#LTdN/d#it{y}#GT";
    h2Yields[iPart] = new TProfile2D(name, title, 100, 0., 20., (gMomentumAxes[kYM].max-gMomentumAxes[kYM].min)/dy, gMomentumAxes[kYM].min, gMomentumAxes[kYM].max);

    for(Int_t iHarm = 0; iHarm < 2; ++iHarm)
    {
      name = Form("hv%d%s_pT", iHarm + 1, particle.name.c_str());
      title = Form("v_{%d}^{%s} (p_{T});p_{T} (GeV/#it{c});v%d", iHarm + 1, particle.displayName.c_str(), iHarm + 1);
      pVn_pT[iHarm][iPart] = new TProfile(name, name, 50, gMomentumAxes[kPT].min, gMomentumAxes[kPT].max);
      name = Form("hv%d%s_Y", iHarm + 1, particle.name.c_str());
      title = Form("v_{%d}^{%s} (#it{y});#it{y};v%d", iHarm + 1, particle.name.c_str(), iHarm + 1);
      pVn_Y[iHarm][iPart] = new TProfile(name, name, 50, gMomentumAxes[kYM].min, gMomentumAxes[kYM].max);
    }
  }

  hPdg = new TH1D("hPdg", "PDG code;PDG code;nCounts", 7000, -3500, 3500);

  cout << "Initialization finished" << endl;
}

void UnigenQA::FillEventInfo()
{
  double M = event_->GetNpa();
  double B = event_->GetB();
  double psiRP = event_->GetPhi();

  if(fPSDGroupEnergy[0][0] > fPSDMax)
    cout << "EventId = " << event_->GetEventNr() << ", full energy = " << fPSDGroupEnergy[0][0] << endl;

  for(uint pidGroup = 0; pidGroup < fPidGroups.size(); pidGroup++)
  {
    for(auto group : gPSDGroups)
    {
      hPSDGroupEnergy.at(group.id).at(pidGroup)->Fill(fPSDGroupEnergy.at(group.id).at(pidGroup));
      hPSDMultCorr.at(group.id).at(pidGroup)->Fill(M, fPSDGroupEnergy.at(group.id).at(pidGroup));
      hBPSDCorr.at(group.id).at(pidGroup)->Fill(fPSDGroupEnergy.at(group.id).at(pidGroup), B);
    }

    for(uint psdCorr = 0; psdCorr < gPSDCorr.size(); ++psdCorr)
    {
      auto group1 = gPSDCorr[psdCorr][0];
      auto group2 = gPSDCorr[psdCorr][1];
      hPSDGroupsCorr.at(psdCorr).at(pidGroup)->Fill(fPSDGroupEnergy.at(group1).at(pidGroup),
                                                    fPSDGroupEnergy.at(group2).at(pidGroup));
    }
  }

  hM->Fill(M);
  hB->Fill(B);
  hPsi->Fill(psiRP);
  hMBcorr->Fill(M, B);
  h2BAmax->Fill(B, fAmax);
  h2BZmax->Fill(B, fZmax);
}

void UnigenQA::FillTracks()
{
  UParticle *track;
  double psiRP = event_->GetPhi();
  Int_t nTracks = event_->GetNpa();
  TLorentzVector momentum;
  TH2I yield("yield","",kParticles+1, 0, kParticles+1, h2Yields[0]->GetNbinsY(), h2Yields[0]->GetYaxis()->GetXmin(), h2Yields[0]->GetYaxis()->GetXmax());
  Int_t pdg, A, Z, flowSign;
  double y, theta, flow, Elab, Ecm, Pcm, Plab, PzLab, Mcm, Mlab;

  fAmax = 0;
  fZmax = 0;
  for(auto group : gPSDGroups)
  {
    for(uint pidGroup = 0; pidGroup < fPidGroups.size(); pidGroup++)
    {
      fPSDGroupEnergy.at(group.id).at(pidGroup) = 0.;
    }
  }

  for(int i = 0; i < nTracks; i++)
  {
    track = event_->GetParticle(i);
    pdg = track->GetPdg();
    if(pdg / 1000000000 != 0)
    {
      A = abs(pdg % 10000 / 10);
      Z = abs(pdg % 10000000 / 10000);
    }
    else if(pdg == 22)
    {
      A = 0;
      Z = 0;
    }
    else
    {
      A = 1;
      Z = 1;
    }
    if(A > fAmax)
      fAmax = A;
    if(Z > fZmax)
      fZmax = Z;
    momentum = track->GetMomentum();
    y = momentum.Rapidity();
    Ecm = momentum.E();
    Pcm = momentum.P();
    Mcm = sqrt(fabs(Ecm * Ecm - Pcm * Pcm));
    momentum.Boost(0., 0., fBeta);
    theta = momentum.Theta();
    Elab = momentum.E();
    Plab = momentum.P();
    PzLab = momentum.Pz();
    Mlab = sqrt(fabs(Elab * Elab - Plab * Plab));

    if(abs(pdg) < 3500)
      hPdg->Fill(pdg);

    double mom[kAxes];
    mom[kPT] = momentum.Pt();
    mom[kETA] = momentum.PseudoRapidity();
    mom[kPHI] = momentum.Phi();
    mom[kYM] = y;
    mom[kEcm] = Ecm;
    mom[kElab] = Elab;
    mom[kPcm] = Pcm;
    mom[kPlab] = Plab;
    mom[kPzLab] = PzLab;
    mom[kMcm] = Mcm;
    // mom [kMlab] = Mlab;
    mom[kA] = A;
    mom[kZ] = Z;
    // mom [kMpdg] = A * 0.931;
    // mom [kMcm_Ecm] = Mcm / Ecm;
    // mom [kMlab_Elab] = Mlab / Elab;

    //				if (mom [kPlab] > gMomentumAxes[kPlab].max || mom [kElab] > gMomentumAxes[kElab].max) cout << i << "\t" <<
    //mom [kA] << "\t" << mom [kPcm] << "\t" << mom [kEcm] << endl;

    for(Int_t iMom = 0; iMom < kAxes; ++iMom)
    {
      hTrackMomentum[iMom][kALLSPECIES]->Fill(mom[iMom]);
      //						if (mom[iMom] > gMomentumAxes[iMom].max) cout << i << "\t" <<
      //gMomentumAxes[iMom].name << " = " << mom[iMom] << endl;
    }

    for(uint j = 0; j < gTH2Axes.size(); ++j)
    {
      vector<int> axes = gTH2Axes[j];
      hTrackMomentumCorr[j][kALLSPECIES]->Fill(mom[axes[0]], mom[axes[1]]);
    }

    for(Int_t iHarm = 0; iHarm < 2; ++iHarm)
    {
      if(iHarm % 2 == 0 && y < 0.)
        flowSign = -1.;
      else
        flowSign = 1.;
      flow = cos((iHarm + 1) * (mom[kPHI] - psiRP));
      pVn_Y[iHarm][kALLSPECIES]->Fill(mom[kYM], flow);
      if (fabs(y) < 1.)
        pVn_pT[iHarm][kALLSPECIES]->Fill(mom[kPT], flowSign * flow);
    }
    yield.Fill(kALLSPECIES, mom[kYM]);

    for(auto group : gPSDGroups)
    {
      for(uint pidGroup = 0; pidGroup < fPidGroups.size(); pidGroup++)
      {
        if(theta >= group.theta[0] && theta < group.theta[1])
          if(abs(pdg) > fPidGroups.at(pidGroup)[0] && abs(pdg) < fPidGroups.at(pidGroup)[1])
            fPSDGroupEnergy.at(group.id).at(pidGroup) += Elab;
      }
    }

    for(Int_t iPart = 0; iPart < kParticles; ++iPart)
    {
      if(iPart == kLEPTONS && abs(pdg) < 38 && pdg != 22)
        pdg = 99999;
      if(iPart == kFRAGMENTS && pdg / 100000000 != 0)
        pdg = 999999;

      if(gParticles[iPart].pdg == pdg)
      {
        yield.Fill(iPart, mom[kYM]);

        for(Int_t iMom = 0; iMom < kAxes; ++iMom)
        {
          hTrackMomentum[iMom][iPart]->Fill(mom[iMom]);
        }

        for(uint iCorr = 0; iCorr < gTH2Axes.size(); ++iCorr)
        {
          vector<int> axes = gTH2Axes[iCorr];
          hTrackMomentumCorr[iCorr][iPart]->Fill(mom[axes[0]], mom[axes[1]]);
        }

        for(Int_t iHarm = 0; iHarm < 2; ++iHarm)
        {
          if(iHarm % 2 == 0 && y < 0.)
            flowSign = -1.;
          else
            flowSign = 1.;
          flow = cos((iHarm + 1) * (mom[kPHI] - psiRP));
          pVn_Y[iHarm][iPart]->Fill(mom[kYM], flow);
          if (fabs(y) < 1.)
            pVn_pT[iHarm][iPart]->Fill(mom[kPT], flowSign * flow);
        }
      }
      if(iPart == kLEPTONS || iPart == kFRAGMENTS)
        pdg = track->GetPdg();
    }
  }

  for(Int_t iPart = 0; iPart < kParticles; iPart++)
    for(int iYbin=1; iYbin<=yield.GetNbinsY(); iYbin++)
      h2Yields[iPart]->Fill(event_->GetB(), yield.GetYaxis()->GetBinCenter(iYbin), yield.GetBinContent(iPart+1,iYbin));
}

void UnigenQA::Write_Histograms(const TString filename)
{
  auto outputFile = new TFile(filename, "RECREATE");
  if(outputFile->IsOpen())
    std::cout << "File '" << filename << "' is opened successfully" << std::endl;
  TDirectory *outputDir;

  hM->Write();
  hB->Write();
  hPsi->Write();
  hMBcorr->Write();
  h2BAmax->Write();
  h2BZmax->Write();
  hPdg->Write();

  for(auto psdGroup : gPSDGroups)
  {
    if(fPidGroupNames.size() > 1)
    {
      outputDir = outputFile->mkdir(psdGroup.name.c_str());
      outputDir->cd();
    }
    for(auto hist : hPSDGroupEnergy[psdGroup.id])
      hist->Write();
  }
  for(auto psdGroup : gPSDGroups)
  {
    if(fPidGroupNames.size() > 1)
    {
      outputDir = outputFile->mkdir(Form("%s_Mult", psdGroup.name.c_str()));
      outputDir->cd();
    }
    for(auto hist : hPSDMultCorr[psdGroup.id])
      hist->Write();
  }
  for(auto psdGroup : gPSDGroups)
  {
    if(fPidGroupNames.size() > 1)
    {
      outputDir = outputFile->mkdir(Form("%s_B", psdGroup.name.c_str()));
      outputDir->cd();
    }
    for(auto hist : hBPSDCorr[psdGroup.id])
      hist->Write();
  }

  for(uint psdCorr = 0; psdCorr < gPSDCorr.size(); psdCorr++)
  {
    auto group1 = gPSDGroups[gPSDCorr[psdCorr][0]];
    auto group2 = gPSDGroups[gPSDCorr[psdCorr][1]];
    if(fPidGroupNames.size() > 1)
    {
      outputDir = outputFile->mkdir((group1.name + group2.name).c_str());
      outputDir->cd();
    }
    for(auto hist : hPSDGroupsCorr.at(psdCorr))
      hist->Write();
  }

  for(auto axis : gMomentumAxes)
  {
    outputDir = outputFile->mkdir(axis.name);
    outputDir->cd();
    for(auto hist : hTrackMomentum[axis.id])
      hist->Write();
  }

  for(uint iCorr = 0; iCorr < gTH2Axes.size(); ++iCorr)
  {
    vector<int> axes = gTH2Axes[iCorr];
    auto xAxis = gMomentumAxes[axes[0]];
    auto yAxis = gMomentumAxes[axes[1]];
    outputDir = outputFile->mkdir(xAxis.name + yAxis.name);
    outputDir->cd();
    for(uint iPart = 0; iPart < kParticles; iPart++)
      hTrackMomentumCorr[iCorr][iPart]->Write();
  }

  for(Int_t iHarm = 0; iHarm < 2; ++iHarm)
  {
    outputDir = outputFile->mkdir(Form("V%i_pT", iHarm + 1));
    outputDir->cd();
    for(auto hist : pVn_pT[iHarm])
      hist->Write();
    outputDir = outputFile->mkdir(Form("V%i_Y", iHarm + 1));
    outputDir->cd();
    for(auto hist : pVn_Y[iHarm])
      hist->Write();
  }

  outputDir = outputFile->mkdir("Yields");
  outputDir->cd();
  for(auto hist : h2Yields)
    hist->Write();

  outputFile->Close();
}
