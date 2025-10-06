#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <iostream>

void plot_hDeltaPhiPurityCorrected_withUnfoldData2018() {
    // Open ROOT file
    TFile *f = TFile::Open("SSC36enhanced_superchic2018_dielectrons_full.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve histograms
    TH1F* hPurity = (TH1F*)f->Get("hDeltaPhiPurityCorrected;1");
    TH1F* hUnfold  = (TH1F*)f->Get("UnfoldData2018");
    if (!hPurity || !hUnfold) {
        std::cerr << "Could not find one of the histograms!" << std::endl;
        f->Close();
        return;
    }

    // --- Print integrals and statistical uncertainties ---
    double errPurity = 0.0;
    double errUnfold = 0.0;

    double integralPurity = hPurity->IntegralAndError(1, hPurity->GetNbinsX(), errPurity);
    double integralUnfold = hUnfold->IntegralAndError(1, hUnfold->GetNbinsX(), errUnfold);

    std::cout << "Integral of hDeltaPhiPurityCorrected: "
              << integralPurity << " ± " << errPurity << std::endl;
    std::cout << "Integral of UnfoldData2018: "
              << integralUnfold << " ± " << errUnfold << std::endl;

    // Style for hDeltaPhiPurityCorrected (blue line + black errors)
    hPurity->SetStats(0);
    hPurity->SetLineColor(kBlue);
    hPurity->SetLineWidth(2);
    hPurity->SetMarkerStyle(20);
    hPurity->SetMarkerColor(kBlue);
    hPurity->SetMarkerSize(1.0);
    hPurity->SetMinimum(0);
    hPurity->SetTitle(";#Delta#phi_{reco};Entries");
    hPurity->GetXaxis()->SetTitleFont(62);
    hPurity->GetXaxis()->SetTitleSize(0.04);
    hPurity->GetYaxis()->SetTitleFont(62);
    hPurity->GetYaxis()->SetTitleSize(0.04);

    // Style for UnfoldData2018 (red line + markers)
    hUnfold->SetLineColor(kRed);
    hUnfold->SetLineWidth(2);
    hUnfold->SetMarkerStyle(21);
    hUnfold->SetMarkerColor(kRed);
    hUnfold->SetMarkerSize(1.0);

    // Canvas
    TCanvas *c = new TCanvas("c", "hDeltaPhiPurityCorrected + UnfoldData2018", 800, 600);

    // Draw histograms
    hPurity->Draw("HIST P");     // Draw blue histogram with line + markers
    hPurity->SetLineColor(kBlack);
    hPurity->SetMarkerColor(kBlack);
    hPurity->Draw("E SAME");     // Add black error bars

    hUnfold->Draw("HIST P SAME"); // Draw UnfoldData2018 on same canvas

    // Add legend
    TLegend *legend = new TLegend(0.65, 0.45, 0.88, 0.62);
    legend->AddEntry(hPurity, "Purity Corrected Data", "lep");
    legend->AddEntry(hUnfold, "Unfolded 2018 Data", "lep");
    legend->Draw();

    // CMS labels
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(62);
    latex.SetTextSize(0.04);
    latex.SetTextAlign(13);
    latex.DrawLatex(0.10, 0.93, "CMS");

    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.16, 0.93, "#it{work in progress}");

    latex.SetTextSize(0.038);
    latex.SetTextAlign(33);
    latex.DrawLatex(0.89, 0.935, "PbPb 2018 #sqrt{#it{s}_{NN}} = 5.02 TeV");

    // --- (Optional) Draw integrals + uncertainties on the plot ---
    latex.SetTextSize(0.03);
    latex.SetTextAlign(11);
    latex.DrawLatex(0.15, 0.85, Form("Purity int: %.2f #pm %.2f", integralPurity, errPurity));
    latex.DrawLatex(0.15, 0.80, Form("Unfold int: %.2f #pm %.2f", integralUnfold, errUnfold));

    c->Update();
    c->SaveAs("hDeltaPhiPurityCorrected_withUnfoldData2018.pdf");

    f->Close();
}

