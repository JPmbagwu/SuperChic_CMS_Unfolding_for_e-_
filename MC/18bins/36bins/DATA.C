#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <iostream>

void plot_UnfoldData2018_withErrors() {
    // Open ROOT file
    TFile *f = TFile::Open("FDC36enhanced_superchic2018_dielectrons_full.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve histogram
    TH1D* h = (TH1D*)f->Get("UnfoldData2018;1");
    if (!h) {
        std::cerr << "Could not find UnfoldData2018!" << std::endl;
        f->Close();
        return;
    }

    // Remove stats box
    h->SetStats(0);

    // Style for histogram line & markers (blue)
    h->SetLineColor(kBlue);
    h->SetLineWidth(2);
    h->SetMarkerStyle(20);
    h->SetMarkerColor(kBlue);
    h->SetMarkerSize(1.0);

    // Set axis titles bold & font size, Y-axis labeled Entries
    h->SetTitle(";#Delta#phi_{gen};Entries");
    h->GetXaxis()->SetTitleFont(62);  // bold title
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetXaxis()->SetLabelFont(62);  // bold tick labels
    h->GetXaxis()->SetLabelSize(0.035);

    h->GetYaxis()->SetTitleFont(62);  // bold title
    h->GetYaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetLabelFont(62);  // bold tick labels
    h->GetYaxis()->SetLabelSize(0.035);
    h->GetYaxis()->SetTitle("Entries");

    // Start y-axis from 0
    h->SetMinimum(0);

    // Create canvas
    TCanvas *c = new TCanvas("c", "UnfoldData2018", 800, 600);

    // Draw line and markers first (blue)
    h->Draw("HIST P");

    // Draw error bars only in black
    h->SetLineColor(kBlack);
    h->SetMarkerColor(kBlack);
    h->Draw("E SAME");

    // Draw CMS + labels
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

    // Update canvas and save
    c->Update();
    c->SaveAs("UnfoldData2018_withErrors.pdf");

    f->Close();
}

