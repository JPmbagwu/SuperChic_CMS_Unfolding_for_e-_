#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TMath.h>
#include <TLegend.h>
#include <TChain.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

// Centralized cut definitions
struct Cuts {
    static constexpr float MaxRapidity = 2.4;        // Individual electron |eta| < 2.4
    static constexpr float SingleElePtMin = 0.0;     // Single electron pT > 0.0 GeV/c
    static constexpr float maxDielePt = 1.0;         // Dielectron pT < 1.0 GeV/c
    static constexpr float minMass = 7.0;            // Minimum dielectron mass (GeV/c²)
    static constexpr float maxMass = 10.0;           // Maximum dielectron mass (GeV/c²)
    static constexpr float maxHFpCut = 7.3;          // Max HF+ energy (GeV)
    static constexpr float maxHFmCut = 7.6;          // Max HF- energy (GeV)
    static constexpr float maxVertexZ = 20.0;        // Max vertex z (cm)
    static constexpr float electronMass = 0.000511;  // Electron mass (GeV)
};

// Function to create TChain from file or list
TChain* CreateChain(const std::string& inFile, const char* treeName) {
    TChain* chain = new TChain(treeName);
    if (inFile.find(".txt") != std::string::npos) {
        std::ifstream file(inFile);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file list " << inFile << std::endl;
            return chain;
        }
        std::string rootFile;
        while (std::getline(file, rootFile)) {
            if (!rootFile.empty()) {
                chain->Add(rootFile.c_str());
            }
        }
    } else {
        chain->Add(inFile.c_str());
    }
    return chain;
}

// Function to compute the angle between two vectors in the transverse plane
float GetAngle(float vPX1, float vPY1, float vPX2, float vPY2) {
    TVector2 DielectronVect(vPX1, vPY1);
    TVector2 ElectronVect(vPX2, vPY2);

    float DielectronMag = sqrt(DielectronVect.X() * DielectronVect.X() + DielectronVect.Y() * DielectronVect.Y());
    float ElectronMag = sqrt(ElectronVect.X() * ElectronVect.X() + ElectronVect.Y() * ElectronVect.Y());

    if (DielectronMag == 0 || ElectronMag == 0) return 0; // Handle zero magnitude case

    float cosAngle = (DielectronVect.X() * ElectronVect.X() + DielectronVect.Y() * ElectronVect.Y()) / (DielectronMag * ElectronMag);
    float sinAngle = (DielectronVect.X() * ElectronVect.Y() - DielectronVect.Y() * ElectronVect.X()) / (DielectronMag * ElectronMag);
    float tanAngle = atan2(sinAngle, cosAngle);
    return tanAngle;
}

void PlotDataAndSuperchicNorm() {
    // File and tree names
    struct Dataset {
        const char* filename;
        const char* treename;
        const char* label;
    };
    std::vector<Dataset> datasets = {
        {"/eos/user/j/jmbagwu/Data2018/merged/merged_ntuples_*.root", "ggHiNtuplizer/EventTree", "Data"},
        {"/eos/user/j/jmbagwu/SC2018/merged_ntuples_*.root", "ggHiNtuplizer/EventTree", "Superchic"}
    };

    // Histogram definitions
    int nBins = 100;
    float ptMin = 0, ptMax = 25.0;
    float dielectronPtMin = 0, dielectronPtMax = 2.0;
    float etaMin = -2.5, etaMax = 2.5;
    float phiMin = -TMath::Pi(), phiMax = TMath::Pi();
    float mMin = 0, mMax = 100.0;
    float deltaPhiMin = -TMath::Pi(), deltaPhiMax = TMath::Pi();
    float angleMin = -TMath::Pi(), angleMax = TMath::Pi();
    float yMin = -3.0, yMax = 3.0;  // rapidity range

    // Arrays to store histograms for Data and Superchic, with nocuts (0) and cuts (1)
    TH1F* hElePt[2][2];
    TH1F* hEleEta[2][2];
    TH1F* hElePhi[2][2];
    TH1F* hPosPt[2][2];
    TH1F* hPosEta[2][2];
    TH1F* hPosPhi[2][2];
    TH1F* hSumPt[2][2];
    TH1F* hSumEta[2][2];
    TH1F* hSumPhi[2][2];
    TH1F* hSumM[2][2];
    TH1F* hSumRapidity[2][2];
    TH1F* hDeltaPhi[2][2];
    TH1F* hAngle[2][2];
    TH2F* hElePtPhi[2][2];
    TH2F* hPosPtPhi[2][2];
    TH2F* hElePtEta[2][2];
    TH2F* hPosPtEta[2][2];
    TH2F* hDiePtM[2][2];
    TH2F* hSumRapidityMass[2][2];

    // Colors and styles for Data and Superchic
    int colors[2] = {kBlack, kRed};
    int markerStyles[2] = {20, 24};  // Filled circle for Data, open circle for Superchic

    // Initialize histograms
    for (int i = 0; i < 2; ++i) {
        for (int c = 0; c < 2; ++c) {
            std::string cutSuffix = (c == 1 ? "_cuts" : "_nocuts");
            std::string suffix = cutSuffix + "_" + std::string(datasets[i].label);
            hElePt[i][c]  = new TH1F(("hElePt" + suffix).c_str(),  ";Reco Level Electron p_{T} (GeV);Normalized Entries", nBins, ptMin, ptMax);
            hEleEta[i][c] = new TH1F(("hEleEta" + suffix).c_str(), ";Reco Level Electron #eta;Normalized Entries", nBins, etaMin, etaMax);
            hElePhi[i][c] = new TH1F(("hElePhi" + suffix).c_str(), ";Reco Level Electron #phi;Normalized Entries", nBins, phiMin, phiMax);
            hPosPt[i][c]  = new TH1F(("hPosPt" + suffix).c_str(),  ";Reco Level Positron p_{T} (GeV);Normalized Entries", nBins, ptMin, ptMax);
            hPosEta[i][c] = new TH1F(("hPosEta" + suffix).c_str(), ";Reco Level Positron #eta;Normalized Entries", nBins, etaMin, etaMax);
            hPosPhi[i][c] = new TH1F(("hPosPhi" + suffix).c_str(), ";Reco Level Positron #phi;Normalized Entries", nBins, phiMin, phiMax);
            hSumPt[i][c]  = new TH1F(("hSumPt" + suffix).c_str(),  ";Reco Level Dielectron p_{T} (GeV);Normalized Entries", nBins, dielectronPtMin, dielectronPtMax);
            hSumEta[i][c] = new TH1F(("hSumEta" + suffix).c_str(), ";Reco Level Dielectron #eta;Normalized Entries", nBins, etaMin, etaMax);
            hSumPhi[i][c] = new TH1F(("hSumPhi" + suffix).c_str(), ";Reco Level Dielectron #phi;Normalized Entries", nBins, phiMin, phiMax);
            hSumM[i][c]   = new TH1F(("hSumM" + suffix).c_str(),   ";Reco Level Dielectron Mass (GeV);Normalized Entries", nBins, mMin, mMax);
            hSumRapidity[i][c] = new TH1F(("hSumRapidity" + suffix).c_str(), ";Reco Level Dielectron Rapidity;Normalized Entries", nBins, yMin, yMax);
            hDeltaPhi[i][c] = new TH1F(("hDeltaPhi" + suffix).c_str(), ";Reco Level #Delta#phi = #phi_{ee} - #phi_{e-} (rad);Normalized Entries", nBins, deltaPhiMin, deltaPhiMax);
            hAngle[i][c]    = new TH1F(("hAngle" + suffix).c_str(),    ";Reco Level #Delta#phi = #phi_{ee} - #phi_{e-} (rad);Normalized Entries", nBins, angleMin, angleMax);
            hElePtPhi[i][c] = new TH2F(("hElePtPhi" + suffix).c_str(), ";Reco Level Electron #phi;Reco Level Electron p_{T} (GeV)", nBins, phiMin, phiMax, nBins, ptMin, ptMax);
            hPosPtPhi[i][c] = new TH2F(("hPosPtPhi" + suffix).c_str(), ";Reco Level Positron #phi;Reco Level Positron p_{T} (GeV)", nBins, phiMin, phiMax, nBins, ptMin, ptMax);
            hElePtEta[i][c] = new TH2F(("hElePtEta" + suffix).c_str(), ";Reco Level Electron #eta;Reco Level Electron p_{T} (GeV)", nBins, etaMin, etaMax, nBins, ptMin, ptMax);
            hPosPtEta[i][c] = new TH2F(("hPosPtEta" + suffix).c_str(), ";Reco Level Positron #eta;Reco Level Positron p_{T} (GeV)", nBins, etaMin, etaMax, nBins, ptMin, ptMax);
            hDiePtM[i][c]   = new TH2F(("hDiePtM" + suffix).c_str(),   ";Reco Level Dielectron Mass (GeV);Reco Level Dielectron p_{T} (GeV)", nBins, mMin, mMax, nBins, dielectronPtMin, dielectronPtMax);
            hSumRapidityMass[i][c] = new TH2F(("hSumRapidityMass" + suffix).c_str(), ";Reco Level Dielectron Rapidity;Dielectron Mass (GeV)", nBins, yMin, yMax, nBins, mMin, mMax);

            // Set histogram styles for 1D
            for (auto h : {hElePt[i][c], hEleEta[i][c], hElePhi[i][c], hPosPt[i][c], hPosEta[i][c], hPosPhi[i][c], hSumPt[i][c], hSumEta[i][c], hSumPhi[i][c], hSumM[i][c], hSumRapidity[i][c], hDeltaPhi[i][c], hAngle[i][c]}) {
                h->SetStats(0);
                h->SetLineColor(colors[i]);
                h->SetLineWidth(2);
                h->SetMarkerStyle(markerStyles[i]);
                h->SetMarkerColor(colors[i]);
            }

            // Set histogram styles for 2D
            for (auto h2 : {hElePtPhi[i][c], hPosPtPhi[i][c], hElePtEta[i][c], hPosPtEta[i][c], hDiePtM[i][c], hSumRapidityMass[i][c]}) {
                h2->SetStats(0);
                h2->SetMarkerStyle(markerStyles[i]);
                h2->SetMarkerColor(colors[i]);
            }
        }
    }

    // Process each dataset
    for (int i = 0; i < 2; ++i) {
        TChain *tree = CreateChain(datasets[i].filename, datasets[i].treename);
        if (tree->GetEntries() == 0) {
            std::cerr << "Error: No entries in tree for " << datasets[i].label << std::endl;
            delete tree;
            continue;
        }

        // Set branch addresses
        std::vector<float> *pt = nullptr, *eta = nullptr, *phi = nullptr, *trkvz = nullptr;
        std::vector<int> *charge = nullptr;
        std::vector<float> *CaloTower_e = nullptr, *CaloTower_eta = nullptr;
        int nEle = 0, nTrk = 0, nTower = 0;
        tree->SetBranchAddress("elePt", &pt);
        tree->SetBranchAddress("eleEta", &eta);
        tree->SetBranchAddress("elePhi", &phi);
        tree->SetBranchAddress("eleCharge", &charge);
        tree->SetBranchAddress("trkvz", &trkvz);
        tree->SetBranchAddress("nEle", &nEle);
        tree->SetBranchAddress("nTrk", &nTrk);
        tree->SetBranchAddress("CaloTower_e", &CaloTower_e);
        tree->SetBranchAddress("CaloTower_eta", &CaloTower_eta);
        tree->SetBranchAddress("nTower", &nTower);

        // Loop over events
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t j = 0; j < nEntries; ++j) {
            tree->GetEntry(j);

            if (!pt || !eta || !phi || !charge || !trkvz || pt->size() < 2 || eta->size() < 2 || phi->size() < 2 || charge->size() < 2 || trkvz->size() == 0) continue;
            if (nEle != 2 || nTrk != 2) continue;
            if (charge->at(0) * charge->at(1) != -1) continue; // Opposite charge requirement

            // Assign electron and positron based on charge
            int EleP = (charge->at(0) == 1) ? 0 : 1;
            int EleN = 1 - EleP;

            // Construct TLorentzVectors
            TLorentzVector ele, pos, sum;
            ele.SetPtEtaPhiM(pt->at(EleN), eta->at(EleN), phi->at(EleN), Cuts::electronMass);
            pos.SetPtEtaPhiM(pt->at(EleP), eta->at(EleP), phi->at(EleP), Cuts::electronMass);
            sum = ele + pos;

            // Compute DeltaPhi using the electron (e-, ele)
            float DielePx = sum.Pt() * cos(sum.Phi());
            float DielePy = sum.Pt() * sin(sum.Phi());
            float ElePx = ele.Pt() * cos(ele.Phi());
            float ElePy = ele.Pt() * sin(ele.Phi());
            float deltaPhi = GetAngle(DielePx, DielePy, ElePx, ElePy);
            float angle = deltaPhi; // Ensure angle is consistent with deltaPhi

            // Fill nocuts histograms
            hElePt[i][0]->Fill(pt->at(EleN));
            hEleEta[i][0]->Fill(eta->at(EleN));
            hElePhi[i][0]->Fill(phi->at(EleN));
            hPosPt[i][0]->Fill(pt->at(EleP));
            hPosEta[i][0]->Fill(eta->at(EleP));
            hPosPhi[i][0]->Fill(phi->at(EleP));
            hSumPt[i][0]->Fill(sum.Pt());
            hSumEta[i][0]->Fill(sum.Eta());
            hSumPhi[i][0]->Fill(sum.Phi());
            hSumM[i][0]->Fill(sum.M());
            hSumRapidity[i][0]->Fill(sum.Rapidity());
            hDeltaPhi[i][0]->Fill(deltaPhi);
            hAngle[i][0]->Fill(angle);
            hElePtPhi[i][0]->Fill(phi->at(EleN), pt->at(EleN));
            hPosPtPhi[i][0]->Fill(phi->at(EleP), pt->at(EleP));
            hElePtEta[i][0]->Fill(eta->at(EleN), pt->at(EleN));
            hPosPtEta[i][0]->Fill(eta->at(EleP), pt->at(EleP));
            hDiePtM[i][0]->Fill(sum.M(), sum.Pt());
            hSumRapidityMass[i][0]->Fill(sum.Rapidity(), sum.M());

            // Apply cuts
            if (pt->at(EleN) <= Cuts::SingleElePtMin || pt->at(EleP) <= Cuts::SingleElePtMin) continue;
            if (fabs(eta->at(EleN)) >= Cuts::MaxRapidity || fabs(eta->at(EleP)) >= Cuts::MaxRapidity) continue;
            if (sum.M() < Cuts::minMass || sum.M() > Cuts::maxMass) continue;
            if (sum.Pt() >= Cuts::maxDielePt) continue;
            if (fabs(trkvz->at(0)) > Cuts::maxVertexZ) continue;

            // Calculate max HF energies
            float maxHFp = 0, maxHFm = 0;
            if (CaloTower_eta && CaloTower_e) {
                for (int t = 0; t < nTower; t++) {
                    float tower_eta = CaloTower_eta->at(t);
                    float tower_E = CaloTower_e->at(t);
                    if (tower_eta > 2.9 && tower_eta < 5.2)
                        maxHFp = TMath::Max(maxHFp, tower_E);
                    if (tower_eta < -2.9 && tower_eta > -5.2)
                        maxHFm = TMath::Max(maxHFm, tower_E);
                }
            }
            if (maxHFp > Cuts::maxHFpCut || maxHFm > Cuts::maxHFmCut) continue;

            // Fill cuts histograms
            hElePt[i][1]->Fill(pt->at(EleN));
            hEleEta[i][1]->Fill(eta->at(EleN));
            hElePhi[i][1]->Fill(phi->at(EleN));
            hPosPt[i][1]->Fill(pt->at(EleP));
            hPosEta[i][1]->Fill(eta->at(EleP));
            hPosPhi[i][1]->Fill(phi->at(EleP));
            hSumPt[i][1]->Fill(sum.Pt());
            hSumEta[i][1]->Fill(sum.Eta());
            hSumPhi[i][1]->Fill(sum.Phi());
            hSumM[i][1]->Fill(sum.M());
            hSumRapidity[i][1]->Fill(sum.Rapidity());
            hDeltaPhi[i][1]->Fill(deltaPhi);
            hAngle[i][1]->Fill(angle);
            hElePtPhi[i][1]->Fill(phi->at(EleN), pt->at(EleN));
            hPosPtPhi[i][1]->Fill(phi->at(EleP), pt->at(EleP));
            hElePtEta[i][1]->Fill(eta->at(EleN), pt->at(EleN));
            hPosPtEta[i][1]->Fill(eta->at(EleP), pt->at(EleP));
            hDiePtM[i][1]->Fill(sum.M(), sum.Pt());
            hSumRapidityMass[i][1]->Fill(sum.Rapidity(), sum.M());
        }
        delete tree;
    }

    // Normalize 1D histograms
    for (int i = 0; i < 2; ++i) {
        for (int c = 0; c < 2; ++c) {
            for (auto h : {hElePt[i][c], hEleEta[i][c], hElePhi[i][c], hPosPt[i][c], hPosEta[i][c], hPosPhi[i][c], hSumPt[i][c], hSumEta[i][c], hSumPhi[i][c], hSumM[i][c], hSumRapidity[i][c], hDeltaPhi[i][c], hAngle[i][c]}) {
                if (h->Integral() > 0) {
                    h->Scale(1.0 / h->Integral());
                }
            }
        }
    }

    // Function to draw CMS and PbPb labels
    auto drawCMSLabels = []() {
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(62);
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.11, 0.91, "CMS");
        latex.SetTextFont(52);
        latex.SetTextSize(0.035);
        latex.DrawLatex(0.17, 0.91, "Work in Progress");
        latex.SetTextFont(42);
        latex.SetTextSize(0.035);
        latex.DrawLatex(0.65, 0.91, "PbPb #sqrt{S_{NN}} = 5.02 TeV");
    };

    // Function to draw legend in top-right corner
    auto drawLegend = [](TH1* h1, TH1* h2) {
        TLegend *leg = new TLegend(0.60, 0.75, 0.80, 0.85); // Top-right corner
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(h1, "Data 2018", "lep");
        leg->AddEntry(h2, "2018 Superchic + FSR", "lep");
        leg->Draw();
        return leg;
    };

    // Plot and save 1D histograms
    auto plotHist = [&](TH1F* h1, TH1F* h2, const char* outname, const char* title) {
        TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
        if (std::string(outname).find("dielectron_pt") != std::string::npos) {
            canvas->SetLogy(1);
        }
        double max1 = h1->GetMaximum();
        double max2 = h2->GetMaximum();
        h1->SetMaximum(std::max(max1, max2) * 1.2);
        h1->Draw("PE");
        h2->Draw("PE SAME");
        drawCMSLabels();
        auto leg = drawLegend(h1, h2);
        if (strlen(title) > 0) {
            TLatex *titleLatex = new TLatex();
            titleLatex->SetNDC();
            titleLatex->SetTextFont(42);
            titleLatex->SetTextSize(0.04);
            titleLatex->DrawLatex(0.15, 0.82, title);
        }
        canvas->SaveAs(outname);
        delete leg;
        delete canvas;
    };

    // Plot and save 2D histograms
    auto plot2DHist = [&](TH2F* h1, TH2F* h2, const char* outname, const char* title, bool withCuts) {
        TCanvas *canvas = new TCanvas("canvas2D", "canvas2D", 1600, 600);
        canvas->Divide(2, 1);
        canvas->cd(1);
        h1->Draw("COLZ");
        drawCMSLabels();
        TLatex *label1 = new TLatex();
        label1->SetNDC();
        label1->SetTextFont(42);
        label1->SetTextSize(0.04);
        label1->DrawLatex(0.15, 0.82, "Data 2018");
        label1->SetTextSize(0.03);
        if (withCuts) {
            label1->DrawLatex(0.15, 0.77, Form("|#eta^{e^{#pm}}| < %.1f", Cuts::MaxRapidity));
            label1->DrawLatex(0.15, 0.72, Form("%.1f < M_{ee} < %.1f GeV", Cuts::minMass, Cuts::maxMass));
            label1->DrawLatex(0.15, 0.67, Form("p_{T}^{ee} < %.1f GeV", Cuts::maxDielePt));
            label1->DrawLatex(0.15, 0.62, Form("|z_{vtx}| < %.1f cm", Cuts::maxVertexZ));
            label1->DrawLatex(0.15, 0.57, Form("HF+ < %.1f GeV, HF- < %.1f GeV", Cuts::maxHFpCut, Cuts::maxHFmCut));
        } else {
            label1->DrawLatex(0.15, 0.77, "No cuts applied");
        }
        canvas->cd(2);
        h2->Draw("COLZ");
        drawCMSLabels();
        TLatex *label2 = new TLatex();
        label2->SetNDC();
        label2->SetTextFont(42);
        label2->SetTextSize(0.04);
        label2->DrawLatex(0.15, 0.82, "2018 Superchic + FSR");
        label2->SetTextSize(0.03);
        if (withCuts) {
            label2->DrawLatex(0.15, 0.77, Form("|#eta^{e^{#pm}}| < %.1f", Cuts::MaxRapidity));
            label2->DrawLatex(0.15, 0.72, Form("%.1f < M_{ee} < %.1f GeV", Cuts::minMass, Cuts::maxMass));
            label2->DrawLatex(0.15, 0.67, Form("p_{T}^{ee} < %.1f GeV", Cuts::maxDielePt));
            label2->DrawLatex(0.15, 0.62, Form("|z_{vtx}| < %.1f cm", Cuts::maxVertexZ));
            label2->DrawLatex(0.15, 0.57, Form("HF+ < %.1f GeV, HF- < %.1f GeV", Cuts::maxHFpCut, Cuts::maxHFmCut));
        } else {
            label2->DrawLatex(0.15, 0.77, "No cuts applied");
        }
        if (strlen(title) > 0) {
            TLatex *titleLatex = new TLatex();
            titleLatex->SetNDC();
            titleLatex->SetTextFont(42);
            titleLatex->SetTextSize(0.04);
            titleLatex->DrawLatex(0.15, 0.76, title);
        }
        canvas->SaveAs(outname);
        delete label1;
        delete label2;
        delete canvas;
    };

    // Plot all 1D histograms for cuts
    plotHist(hElePt[0][1], hElePt[1][1], "electron_pt_comparison_cuts.pdf", "");
    plotHist(hEleEta[0][1], hEleEta[1][1], "electron_eta_comparison_cuts.pdf",
             Form("|#eta^{e}| < %.1f", Cuts::MaxRapidity));
    plotHist(hElePhi[0][1], hElePhi[1][1], "electron_phi_comparison_cuts.pdf", "");
    plotHist(hPosPt[0][1], hPosPt[1][1], "positron_pt_comparison_cuts.pdf", "");
    plotHist(hPosEta[0][1], hPosEta[1][1], "positron_eta_comparison_cuts.pdf",
             Form("|#eta^{e^{+}}| < %.1f", Cuts::MaxRapidity));
    plotHist(hPosPhi[0][1], hPosPhi[1][1], "positron_phi_comparison_cuts.pdf", "");
    plotHist(hSumPt[0][1], hSumPt[1][1], "dielectron_pt_comparison_cuts.pdf",
             Form("p_{T}^{ee} < %.1f GeV", Cuts::maxDielePt));
    plotHist(hSumEta[0][1], hSumEta[1][1], "dielectron_eta_comparison_cuts.pdf",
             Form("|#eta^{e^{#pm}}| < %.1f", Cuts::MaxRapidity));
    plotHist(hSumPhi[0][1], hSumPhi[1][1], "dielectron_phi_comparison_cuts.pdf", "");
    plotHist(hSumM[0][1], hSumM[1][1], "dielectron_mass_comparison_cuts.pdf",
             Form("%.1f < M_{ee} < %.1f GeV", Cuts::minMass, Cuts::maxMass));
    plotHist(hSumRapidity[0][1], hSumRapidity[1][1], "dielectron_rapidity_comparison_cuts.pdf",
             Form("|y^{ee}| < %.1f", Cuts::MaxRapidity));
    plotHist(hDeltaPhi[0][1], hDeltaPhi[1][1], "delta_phi_comparison_cuts.pdf", "");
    plotHist(hAngle[0][1], hAngle[1][1], "angle_between_vectors_comparison_cuts.pdf", "");

    // Plot all 1D histograms for nocuts
    plotHist(hElePt[0][0], hElePt[1][0], "electron_pt_comparison_nocuts.pdf", "No cuts");
    plotHist(hEleEta[0][0], hEleEta[1][0], "electron_eta_comparison_nocuts.pdf", "No cuts");
    plotHist(hElePhi[0][0], hElePhi[1][0], "electron_phi_comparison_nocuts.pdf", "No cuts");
    plotHist(hPosPt[0][0], hPosPt[1][0], "positron_pt_comparison_nocuts.pdf", "No cuts");
    plotHist(hPosEta[0][0], hPosEta[1][0], "positron_eta_comparison_nocuts.pdf", "No cuts");
    plotHist(hPosPhi[0][0], hPosPhi[1][0], "positron_phi_comparison_nocuts.pdf", "No cuts");
    plotHist(hSumPt[0][0], hSumPt[1][0], "dielectron_pt_comparison_nocuts.pdf", "No cuts");
    plotHist(hSumEta[0][0], hSumEta[1][0], "dielectron_eta_comparison_nocuts.pdf", "No cuts");
    plotHist(hSumPhi[0][0], hSumPhi[1][0], "dielectron_phi_comparison_nocuts.pdf", "No cuts");
    plotHist(hSumM[0][0], hSumM[1][0], "dielectron_mass_comparison_nocuts.pdf", "No cuts");
    plotHist(hSumRapidity[0][0], hSumRapidity[1][0], "dielectron_rapidity_comparison_nocuts.pdf", "No cuts");
    plotHist(hDeltaPhi[0][0], hDeltaPhi[1][0], "delta_phi_comparison_nocuts.pdf", "No cuts");
    plotHist(hAngle[0][0], hAngle[1][0], "angle_between_vectors_comparison_nocuts.pdf", "No cuts");

    // Plot all 2D histograms for cuts
    plot2DHist(hElePtPhi[0][1], hElePtPhi[1][1], "electron_pt_vs_phi_comparison_cuts.pdf", "", true);
    plot2DHist(hPosPtPhi[0][1], hPosPtPhi[1][1], "positron_pt_vs_phi_comparison_cuts.pdf", "", true);
    plot2DHist(hElePtEta[0][1], hElePtEta[1][1], "electron_pt_vs_eta_comparison_cuts.pdf", "", true);
    plot2DHist(hPosPtEta[0][1], hPosPtEta[1][1], "positron_pt_vs_eta_comparison_cuts.pdf", "", true);
    plot2DHist(hDiePtM[0][1], hDiePtM[1][1], "dielectron_pt_vs_mass_comparison_cuts.pdf", "", true);
    plot2DHist(hSumRapidityMass[0][1], hSumRapidityMass[1][1], "dielectron_rapidity_vs_mass_comparison_cuts.pdf", "", true);

    // Plot all 2D histograms for nocuts
    plot2DHist(hElePtPhi[0][0], hElePtPhi[1][0], "electron_pt_vs_phi_comparison_nocuts.pdf", "", false);
    plot2DHist(hPosPtPhi[0][0], hPosPtPhi[1][0], "positron_pt_vs_phi_comparison_nocuts.pdf", "", false);
    plot2DHist(hElePtEta[0][0], hElePtEta[1][0], "electron_pt_vs_eta_comparison_nocuts.pdf", "", false);
    plot2DHist(hPosPtEta[0][0], hPosPtEta[1][0], "positron_pt_vs_eta_comparison_nocuts.pdf", "", false);
    plot2DHist(hDiePtM[0][0], hDiePtM[1][0], "dielectron_pt_vs_mass_comparison_nocuts.pdf", "", false);
    plot2DHist(hSumRapidityMass[0][0], hSumRapidityMass[1][0], "dielectron_rapidity_vs_mass_comparison_nocuts.pdf", "", false);

    // Cleanup
    for (int i = 0; i < 2; ++i) {
        for (int c = 0; c < 2; ++c) {
            delete hElePt[i][c]; delete hEleEta[i][c]; delete hElePhi[i][c];
            delete hPosPt[i][c]; delete hPosEta[i][c]; delete hPosPhi[i][c];
            delete hSumPt[i][c]; delete hSumEta[i][c]; delete hSumPhi[i][c]; delete hSumM[i][c]; delete hSumRapidity[i][c];
            delete hDeltaPhi[i][c]; delete hAngle[i][c];
            delete hElePtPhi[i][c]; delete hPosPtPhi[i][c];
            delete hElePtEta[i][c]; delete hPosPtEta[i][c];
            delete hDiePtM[i][c]; delete hSumRapidityMass[i][c];
        }
    }
}
