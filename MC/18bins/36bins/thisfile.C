#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include <iostream>
#include <TLine.h>
#include <TPaveText.h>

void compareChi2_withNormalization() {
    gStyle->SetOptStat(0);  // Disable stats box globally
    gStyle->SetOptTitle(0);

    // Load ROOT file
    TFile* file = TFile::Open("FDC36enhanced_superchic2018_dielectrons_full.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open ROOT file!" << std::endl;
        return;
    }

    // Retrieve histograms
    TH1F* hDataReco = (TH1F*)file->Get("hDeltaPhiRaw");
    TH1F* hMCReco   = (TH1F*)file->Get("hrecoTrue");
    TH1D* hDataUnfold = (TH1D*)file->Get("hUnfoldDataOverAcceptance");
    TH1F* hMCGen    = (TH1F*)file->Get("hgenTrue");

    if (!hDataReco || !hMCReco || !hDataUnfold || !hMCGen) {
        std::cerr << "Failed to retrieve one or more histograms!" << std::endl;
        if (!hDataReco) std::cerr << "hDeltaPhiRaw not found!" << std::endl;
        if (!hMCReco) std::cerr << "hrecoTrue not found!" << std::endl;
        if (!hDataUnfold) std::cerr << "hUnfoldDataOverAcceptance not found!" << std::endl;
        if (!hMCGen) std::cerr << "hgenTrue not found!" << std::endl;
        file->Close();
        return;
    }

    // Lambda function to calculate chi2
    auto calculateChi2 = [](TH1* h1, TH1* h2, const std::string& name = "") {
        double chi2 = 0;
        int ndf = 0;
        
        std::cout << "\n=== " << name << " χ² Calculation ===" << std::endl;
        std::cout << "Bin\tData\t±err\tMC\t±err\tχ² contrib" << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
        
        for (int i = 1; i <= h1->GetNbinsX(); i++) {
            double val1 = h1->GetBinContent(i);
            double val2 = h2->GetBinContent(i);
            double err1 = h1->GetBinError(i);
            double err2 = h2->GetBinError(i);

            if (err1 == 0 && err2 == 0) continue;
            if (val1 == 0 && val2 == 0) continue;

            double combined_err = sqrt(err1*err1 + err2*err2);
            if (combined_err > 0) {
                double bin_chi2 = pow((val1 - val2) / combined_err, 2);
                chi2 += bin_chi2;
                ndf++;
                std::cout << i << "\t" << val1 << "\t" << err1 << "\t"
                          << val2 << "\t" << err2 << "\t" << bin_chi2 << std::endl;
            }
        }

        std::cout << "------------------------------------------------" << std::endl;
        std::cout << "Total χ² = " << chi2 << ", NDF = " << ndf << std::endl;
        return std::make_pair(chi2, ndf);
    };

    // ================= Reco-level Comparison =================
    TH1F* hDataNorm = (TH1F*)hDataReco->Clone("hDataNorm");
    TH1F* hMCNorm   = (TH1F*)hMCReco->Clone("hMCNorm");
    hDataNorm->SetStats(0);
    hMCNorm->SetStats(0);

    if (hDataNorm->Integral() > 0) hDataNorm->Scale(1.0 / hDataNorm->Integral());
    if (hMCNorm->Integral() > 0)   hMCNorm->Scale(1.0 / hMCNorm->Integral());

    auto [chi2_reco, ndf_reco] = calculateChi2(hDataNorm, hMCNorm, "Reco Level");
    double pvalue_reco = (ndf_reco > 0) ? TMath::Prob(chi2_reco, ndf_reco) : 0;

    int chi2_reco_int       = TMath::Nint(chi2_reco);
    int chi2_over_18_reco_int = TMath::Nint(chi2_reco / 18.0);
    int chi2_over_36_reco_int = TMath::Nint(chi2_reco / 36.0);

    // Reco plot
    TCanvas* c1 = new TCanvas("c1", "Reco-level Comparison", 1000, 600);
    c1->SetLeftMargin(0.12); c1->SetRightMargin(0.08); c1->SetTopMargin(0.08); c1->SetBottomMargin(0.12);

    hDataNorm->SetLineColor(kRed); hDataNorm->SetLineWidth(3);
    hDataNorm->SetMarkerColor(kRed); hDataNorm->SetMarkerStyle(20); hDataNorm->SetMarkerSize(1.5);
    hMCNorm->SetLineColor(kBlue); hMCNorm->SetLineWidth(3);
    hMCNorm->SetMarkerColor(kBlue); hMCNorm->SetMarkerStyle(24); hMCNorm->SetMarkerSize(1.5);

    double max_val_reco = TMath::Max(hDataNorm->GetMaximum(), hMCNorm->GetMaximum());
    hDataNorm->GetYaxis()->SetRangeUser(0, max_val_reco * 1.4);
    hDataNorm->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e-}");
    hDataNorm->GetYaxis()->SetTitle("Normalized Entries");
    hDataNorm->Draw("E");
    hMCNorm->Draw("E SAME");

    TLatex latex;
    latex.SetNDC(); latex.SetTextFont(62); latex.SetTextSize(0.045);
    latex.DrawLatex(0.12, 0.93, "CMS");
    latex.SetTextFont(52); latex.SetTextSize(0.035);
    latex.DrawLatex(0.18, 0.93, "Work in Progress");
    latex.SetTextFont(42); latex.SetTextSize(0.035);
    latex.DrawLatex(0.65, 0.94, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    latex.SetTextFont(62); latex.SetTextAlign(13); latex.SetTextSize(0.030);
    latex.DrawLatex(0.18, 0.85, "p_{T}_{ee} < 1.0 GeV");
    latex.DrawLatex(0.18, 0.80, "|#eta_{e}| < 2.4");
    latex.DrawLatex(0.18, 0.75, "7.0 < M_{ee} < 10.0 GeV");

    TPaveText* pt_reco = new TPaveText(0.65, 0.55, 0.85, 0.70, "NDC");
    pt_reco->SetFillColor(0); pt_reco->SetBorderSize(1);
    pt_reco->AddText(Form("#chi^{2}/ndf = %d/%d", chi2_reco_int, ndf_reco));
    pt_reco->AddText(Form("#chi^{2}/36 = %d", chi2_over_36_reco_int));
    pt_reco->Draw();

    TLegend* leg1 = new TLegend(0.65, 0.75, 0.85, 0.85);
    leg1->SetBorderSize(0); leg1->SetFillStyle(0);
    leg1->AddEntry(hDataNorm, "Reconstructed Data", "lep");
    leg1->AddEntry(hMCNorm, "Reconstructed SuperChic", "lep");
    leg1->Draw();

    c1->SaveAs("RecoLevelComparison_Normalized.pdf");
    c1->SaveAs("RecoLevelComparison_Normalized.png");

    std::cout << "Reco-level χ²/NDF: " << chi2_reco_int << "/" << ndf_reco
              << ", χ²/36 = " << chi2_over_36_reco_int << ", p = " << pvalue_reco << std::endl;

    // ================= Gen-level Comparison =================
    TH1D* hUnfoldNorm = (TH1D*)hDataUnfold->Clone("hUnfoldNorm");
    TH1F* hGenMCNorm = (TH1F*)hMCGen->Clone("hGenMCNorm");
    hUnfoldNorm->SetStats(0); hGenMCNorm->SetStats(0);

    if (hUnfoldNorm->Integral() > 0) hUnfoldNorm->Scale(1.0 / hUnfoldNorm->Integral());
    if (hGenMCNorm->Integral() > 0) hGenMCNorm->Scale(1.0 / hGenMCNorm->Integral());

    auto [chi2_gen, ndf_gen] = calculateChi2(hUnfoldNorm, hGenMCNorm, "Gen Level");
    double pvalue_gen = (ndf_gen > 0) ? TMath::Prob(chi2_gen, ndf_gen) : 0;

    int chi2_gen_int       = TMath::Nint(chi2_gen);
    int chi2_over_18_gen_int = TMath::Nint(chi2_gen / 18.0);
    int chi2_over_36_gen_int = TMath::Nint(chi2_gen / 36.0);

    TCanvas* c2 = new TCanvas("c2", "Gen-level Comparison", 1000, 600);
    c2->SetLeftMargin(0.12); c2->SetRightMargin(0.08); c2->SetTopMargin(0.08); c2->SetBottomMargin(0.12);

    hUnfoldNorm->SetLineColor(kRed); hUnfoldNorm->SetLineWidth(3);
    hUnfoldNorm->SetMarkerColor(kRed); hUnfoldNorm->SetMarkerStyle(20); hUnfoldNorm->SetMarkerSize(1.5);
    hGenMCNorm->SetLineColor(kBlue); hGenMCNorm->SetLineWidth(3);
    hGenMCNorm->SetMarkerColor(kBlue); hGenMCNorm->SetMarkerStyle(24); hGenMCNorm->SetMarkerSize(1.5);

    double max_val_gen = TMath::Max(hUnfoldNorm->GetMaximum(), hGenMCNorm->GetMaximum());
    hUnfoldNorm->GetYaxis()->SetRangeUser(0, max_val_gen * 1.4);
    hUnfoldNorm->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e-}");
    hUnfoldNorm->GetYaxis()->SetTitle("Normalized Entries");

    hUnfoldNorm->Draw("E");
    hGenMCNorm->Draw("E SAME");

    latex.SetTextFont(62); latex.SetTextSize(0.045);
    latex.DrawLatex(0.12, 0.96, "CMS");
    latex.SetTextFont(52); latex.SetTextSize(0.035);
    latex.DrawLatex(0.18, 0.96, "Work in Progress");
    latex.SetTextFont(42); latex.SetTextSize(0.035);
    latex.DrawLatex(0.65, 0.97, "PbPb #sqrt{s_{NN}} = 5.02 TeV");

    latex.SetTextFont(62); latex.SetTextAlign(13); latex.SetTextSize(0.030);
    latex.DrawLatex(0.18, 0.85, "p_{T}_{ee} < 1.0 GeV");
    latex.DrawLatex(0.18, 0.80, "|#eta_{e}| < 2.4");
    latex.DrawLatex(0.18, 0.75, "7.0 < M_{ee} < 10.0 GeV");

    TPaveText* pt_gen = new TPaveText(0.65, 0.55, 0.85, 0.70, "NDC");
    pt_gen->SetFillColor(0); pt_gen->SetBorderSize(1);
    pt_gen->AddText(Form("#chi^{2}/ndf = %d/%d", chi2_gen_int, ndf_gen));
    pt_gen->AddText(Form("#chi^{2}/36 = %d", chi2_over_36_gen_int));
    pt_gen->Draw();

    TLegend* leg2 = new TLegend(0.65, 0.75, 0.85, 0.85);
    leg2->SetBorderSize(0); leg2->SetFillStyle(0);
    leg2->AddEntry(hUnfoldNorm, "Unfolded Data", "lep");
    leg2->AddEntry(hGenMCNorm, "Generated-Level SuperChic", "lep");
    leg2->Draw();

    c2->SaveAs("GenLevelComparison_Normalized.pdf");
    c2->SaveAs("GenLevelComparison_Normalized.png");

    std::cout << "Gen-level χ²/NDF: " << chi2_gen_int << "/" << ndf_gen
              << ", χ²/36 = " << chi2_over_36_gen_int << ", p = " << pvalue_gen << std::endl;

    // Cleanup
    delete hDataNorm; delete hMCNorm;
    delete hUnfoldNorm; delete hGenMCNorm;
    delete c1; delete c2; delete leg1; delete leg2; delete pt_reco; delete pt_gen;

    file->Close();
    std::cout << "\nBoth plots saved." << std::endl;
}

