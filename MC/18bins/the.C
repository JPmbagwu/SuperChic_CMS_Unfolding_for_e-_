#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include <TF1.h>
#include <iostream>
#include <TLine.h>
#include <tuple>

std::tuple<TF1*, double, double, double, double>
fitCosineFunctionAndReturn(TH1* hist, const char* name) {
    double xmin = hist->GetXaxis()->GetXmin();
    double xmax = hist->GetXaxis()->GetXmax();

    TF1* fitFunc = new TF1(Form("fit_%s", name),
                          "[0]*(1 + [1]*cos(2*x) + [2]*cos(4*x))",
                          xmin, xmax);
    fitFunc->SetParameters(1.0, 0.1, 0.05);
    fitFunc->SetParNames("C", "a_{2}", "a_{4}");

    hist->Fit(fitFunc, "RQN");

    double a2 = fitFunc->GetParameter(1);
    double a2_err = fitFunc->GetParError(1);
    double a4 = fitFunc->GetParameter(2);
    double a4_err = fitFunc->GetParError(2);

    std::cout << "Fit results for " << name << ":" << std::endl;
    std::cout << "  a2 = " << a2 << " ± " << a2_err << std::endl;
    std::cout << "  a4 = " << a4 << " ± " << a4_err << std::endl;

    return std::make_tuple(fitFunc, a2, a2_err, a4, a4_err);
}

void compareChi2_withNormalization() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TFile* file = TFile::Open("SSC18enhanced_superchic2018_dielectrons_full.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open ROOT file!" << std::endl;
        return;
    }

    TH1F* hDataReco = (TH1F*)file->Get("hDeltaPhiRaw");
    TH1F* hMCReco = (TH1F*)file->Get("hreco");
    TH1D* hDataUnfold = (TH1D*)file->Get("hUnfoldDataOverAcceptance");
    TH1F* hMCGen = (TH1F*)file->Get("hgen");

    if (!hDataReco || !hMCReco || !hDataUnfold || !hMCGen) {
        std::cerr << "Failed to retrieve histograms!" << std::endl;
        file->Close();
        return;
    }

    auto calculateChi2 = [](TH1* h1, TH1* h2, const std::string& name = "") {
        double chi2 = 0;
        int ndf = 0;
        for (int i = 1; i <= h1->GetNbinsX(); i++) {
            double val1 = h1->GetBinContent(i);
            double val2 = h2->GetBinContent(i);
            double err1 = h1->GetBinError(i);
            double err2 = h2->GetBinError(i);
            if ((err1 == 0 && err2 == 0) || (val1 == 0 && val2 == 0)) continue;
            double comb = sqrt(err1*err1 + err2*err2);
            if (comb > 0) {
                chi2 += pow((val1 - val2)/comb, 2);
                ndf++;
            }
        }
        return std::make_pair(chi2, ndf);
    };

    // === Reco-level ===
    TH1F* hDataNorm = (TH1F*)hDataReco->Clone("hDataNorm");
    TH1F* hMCNorm = (TH1F*)hMCReco->Clone("hMCNorm");

    hDataNorm->Sumw2();
    hMCNorm->Sumw2();

    hDataNorm->SetStats(0); hMCNorm->SetStats(0);

    double dataInt = hDataNorm->Integral();
    double mcInt = hMCNorm->Integral();
    double binWidth = hDataNorm->GetBinWidth(1);

    if (dataInt > 0 && binWidth > 0) {
        hDataNorm->Scale(2 * TMath::Pi() / (dataInt * binWidth));
    }
    if (mcInt > 0 && binWidth > 0) {
        hMCNorm->Scale(2 * TMath::Pi() / (mcInt * binWidth));
    }

    TF1* fitReco;
    double a2_r, a2e_r, a4_r, a4e_r;
    std::tie(fitReco, a2_r, a2e_r, a4_r, a4e_r) = fitCosineFunctionAndReturn(hDataNorm, "reco_data");

    fitReco->SetLineColor(kRed);
    fitReco->SetLineWidth(2);
    fitReco->SetLineStyle(2);

    auto [chi2_reco, ndf_reco] = calculateChi2(hDataNorm, hMCNorm, "Reco Level");
    int chi2_reco_int = TMath::Nint(chi2_reco);
    int chi2_over_36_reco_int = TMath::Nint(chi2_reco / 36.0);

    TCanvas* c1 = new TCanvas("c1", "Reco-level Comparison", 1000, 600);
    c1->SetLeftMargin(0.12); c1->SetRightMargin(0.08);
    c1->SetTopMargin(0.08); c1->SetBottomMargin(0.12);

    hDataNorm->SetLineColor(kRed); hDataNorm->SetLineWidth(3);
    hDataNorm->SetMarkerColor(kRed); hDataNorm->SetMarkerStyle(20); hDataNorm->SetMarkerSize(1.5);
    hMCNorm->SetLineColor(kBlue); hMCNorm->SetLineWidth(3);
    hMCNorm->SetMarkerColor(kBlue); hMCNorm->SetMarkerStyle(24); hMCNorm->SetMarkerSize(1.5);

    hDataNorm->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e-}");
    hDataNorm->GetXaxis()->SetTitleSize(0.045); hDataNorm->GetXaxis()->SetTitleOffset(1.1);
    hDataNorm->GetXaxis()->SetLabelSize(0.04);
    hDataNorm->GetYaxis()->SetTitle("\\frac{2\\pi}{N} \\frac{dN}{d\\Delta\\phi}");
    hDataNorm->GetYaxis()->SetTitleSize(0.045); hDataNorm->GetYaxis()->SetTitleOffset(1.2);
    hDataNorm->GetYaxis()->SetLabelSize(0.04);

    double max_val_reco = TMath::Max(hDataNorm->GetMaximum(), hMCNorm->GetMaximum());
    hDataNorm->GetYaxis()->SetRangeUser(0, max_val_reco * 1.4);

    hDataNorm->Draw("E");
    hMCNorm->Draw("E SAME");
    fitReco->Draw("SAME");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(62); latex.SetTextSize(0.045); latex.DrawLatex(0.12, 0.93, "CMS");
    latex.SetTextFont(52); latex.SetTextSize(0.035); latex.DrawLatex(0.18, 0.93, "Work in Progress");
    latex.SetTextFont(42); latex.SetTextSize(0.035); latex.DrawLatex(0.65, 0.94, "PbPb #sqrt{s_{NN}} = 5.02 TeV");

    latex.SetTextFont(62); latex.SetTextAlign(13); latex.SetTextSize(0.030);
    latex.DrawLatex(0.18, 0.85, "p_{T}_{ee} < 1.0 GeV");
    latex.DrawLatex(0.18, 0.80, "|#eta_{e}| < 2.4");
    latex.DrawLatex(0.18, 0.75, "7.0 < M_{ee} < 10.0 GeV");

    TPaveText* pt_reco = new TPaveText(0.65, 0.50, 0.88, 0.72, "NDC");
    pt_reco->SetFillColor(0); pt_reco->SetBorderSize(1);
    pt_reco->AddText(Form("#chi^{2}/ndf = %d/%d", chi2_reco_int, ndf_reco));
    pt_reco->AddText(Form("#chi^{2}/36 = %d", chi2_over_36_reco_int));
    pt_reco->AddText(Form("a_{2} = %.3f #pm %.3f", a2_r, a2e_r));
    pt_reco->AddText(Form("a_{4} = %.3f #pm %.3f", a4_r, a4e_r));
    pt_reco->Draw();

    TLegend* leg1 = new TLegend(0.65, 0.75, 0.85, 0.85);
    leg1->SetBorderSize(0); leg1->SetFillStyle(0);
    leg1->AddEntry(hDataNorm, "Reconstructed Data", "lep");
    leg1->AddEntry(fitReco, "Fit: C(1+a_{2}cos2#Delta#phi+a_{4}cos4#Delta#phi)", "l");
    leg1->AddEntry(hMCNorm, "Reconstructed SuperChic", "lep");
    leg1->Draw();

    c1->SaveAs("RecoLevelComparison_Normalizedfit.pdf");
    c1->SaveAs("RecoLevelComparison_Normalizedfit.png");

    // === Gen-level ===
    TH1D* hUnfoldNorm = (TH1D*)hDataUnfold->Clone("hUnfoldNorm");
    TH1F* hGenMCNorm = (TH1F*)hMCGen->Clone("hGenMCNorm");

    hUnfoldNorm->Sumw2();
    hGenMCNorm->Sumw2();

    hUnfoldNorm->SetStats(0); hGenMCNorm->SetStats(0);

    double unfInt = hUnfoldNorm->Integral();
    double genInt = hGenMCNorm->Integral();
    double binWidthGen = hUnfoldNorm->GetBinWidth(1);

    if (unfInt > 0 && binWidthGen > 0) {
        hUnfoldNorm->Scale(2 * TMath::Pi() / (unfInt * binWidthGen));
    }
    if (genInt > 0 && binWidthGen > 0) {
        hGenMCNorm->Scale(2 * TMath::Pi() / (genInt * binWidthGen));
    }

    TF1* fitGen;
    double a2_g, a2e_g, a4_g, a4e_g;
    std::tie(fitGen, a2_g, a2e_g, a4_g, a4e_g) = fitCosineFunctionAndReturn(hUnfoldNorm, "gen_data");

    fitGen->SetLineColor(kRed);
    fitGen->SetLineWidth(2);
    fitGen->SetLineStyle(2);

    auto [chi2_gen, ndf_gen] = calculateChi2(hUnfoldNorm, hGenMCNorm, "Gen Level");
    int chi2_gen_int = TMath::Nint(chi2_gen);
    int chi2_over_36_gen_int = TMath::Nint(chi2_gen / 36.0);

    TCanvas* c2 = new TCanvas("c2", "Gen-level Comparison", 1000, 600);
    c2->SetLeftMargin(0.12); c2->SetRightMargin(0.08);
    c2->SetTopMargin(0.08); c2->SetBottomMargin(0.12);

    hUnfoldNorm->SetLineColor(kRed); hUnfoldNorm->SetLineWidth(3);
    hUnfoldNorm->SetMarkerColor(kRed); hUnfoldNorm->SetMarkerStyle(20); hUnfoldNorm->SetMarkerSize(1.5);
    hGenMCNorm->SetLineColor(kBlue); hGenMCNorm->SetLineWidth(3);
    hGenMCNorm->SetMarkerColor(kBlue); hGenMCNorm->SetMarkerStyle(24); hGenMCNorm->SetMarkerSize(1.5);

    hUnfoldNorm->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e-}");
    hUnfoldNorm->GetXaxis()->SetTitleSize(0.045); hUnfoldNorm->GetXaxis()->SetTitleOffset(1.1);
    hUnfoldNorm->GetXaxis()->SetLabelSize(0.04);
    hUnfoldNorm->GetYaxis()->SetTitle("\\frac{2\\pi}{N} \\frac{dN}{d\\Delta\\phi}");
    hUnfoldNorm->GetYaxis()->SetTitleSize(0.045); hUnfoldNorm->GetYaxis()->SetTitleOffset(1.2);
    hUnfoldNorm->GetYaxis()->SetLabelSize(0.04);

    double max_val_gen = TMath::Max(hUnfoldNorm->GetMaximum(), hGenMCNorm->GetMaximum());
    hUnfoldNorm->GetYaxis()->SetRangeUser(0, max_val_gen * 1.4);

    hUnfoldNorm->Draw("E");
    hGenMCNorm->Draw("E SAME");
    fitGen->Draw("SAME");

    latex.SetNDC();
    latex.SetTextFont(62); latex.SetTextSize(0.045); latex.DrawLatex(0.12, 0.96, "CMS");
    latex.SetTextFont(52); latex.SetTextSize(0.035); latex.DrawLatex(0.18, 0.96, "Work in Progress");
    latex.SetTextFont(42); latex.SetTextSize(0.035); latex.DrawLatex(0.65, 0.97, "PbPb #sqrt{s_{NN}} = 5.02 TeV");

    latex.SetTextFont(62); latex.SetTextAlign(13); latex.SetTextSize(0.030);
    latex.DrawLatex(0.18, 0.85, "p_{T}_{ee} < 1.0 GeV");
    latex.DrawLatex(0.18, 0.80, "|#eta_{e}| < 2.4");
    latex.DrawLatex(0.18, 0.75, "7.0 < M_{ee} < 10.0 GeV");

    TPaveText* pt_gen = new TPaveText(0.65, 0.50, 0.88, 0.72, "NDC");
    pt_gen->SetFillColor(0); pt_gen->SetBorderSize(1);
    pt_gen->AddText(Form("#chi^{2}/ndf = %d/%d", chi2_gen_int, ndf_gen));
    pt_gen->AddText(Form("#chi^{2}/36 = %d", chi2_over_36_gen_int));
    pt_gen->AddText(Form("a_{2} = %.3f #pm %.3f", a2_g, a2e_g));
    pt_gen->AddText(Form("a_{4} = %.3f #pm %.3f", a4_g, a4e_g));
    pt_gen->Draw();

    TLegend* leg2 = new TLegend(0.65, 0.75, 0.85, 0.85);
    leg2->SetBorderSize(0); leg2->SetFillStyle(0);
    leg2->AddEntry(hUnfoldNorm, "Unfolded Data", "lep");
    leg2->AddEntry(fitGen, "Fit: C(1+a_{2}cos2#Delta#phi+a_{4}cos4#Delta#phi)", "l");
    leg2->AddEntry(hGenMCNorm, "Generated-Level SuperChic", "lep");
    leg2->Draw();

    c2->SaveAs("GenLevelComparison_Normalizedfit.pdf");
    c2->SaveAs("GenLevelComparison_Normalizedfit.png");

    // === Double Ratio ===
    TH1F* hRecoRatio = (TH1F*)hMCReco->Clone("hRecoRatio");
    TH1F* hGenRatio = (TH1F*)hMCGen->Clone("hGenRatio");
    TH1F* hDoubleRatio = (TH1F*)hMCReco->Clone("hDoubleRatio");

    hRecoRatio->Sumw2();
    hGenRatio->Sumw2();
    hDoubleRatio->Sumw2();

    hRecoRatio->Divide(hDataReco, hMCReco, 1.0, 1.0, "B");
    hGenRatio->Divide(hDataUnfold, hMCGen, 1.0, 1.0, "B");
    hDoubleRatio->Divide(hRecoRatio, hGenRatio, 1.0, 1.0, "B");

    hDoubleRatio->SetStats(0);

    TCanvas* c3 = new TCanvas("c3", "Double Ratio", 1000, 600);
    c3->SetLeftMargin(0.12); c3->SetRightMargin(0.08);
    c3->SetTopMargin(0.08); c3->SetBottomMargin(0.12);

    hDoubleRatio->SetLineColor(kBlack); hDoubleRatio->SetLineWidth(3);
    hDoubleRatio->SetMarkerColor(kBlack); hDoubleRatio->SetMarkerStyle(20); hDoubleRatio->SetMarkerSize(1.5);

    hDoubleRatio->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e-}");
    hDoubleRatio->GetXaxis()->SetTitleSize(0.045); hDoubleRatio->GetXaxis()->SetTitleOffset(1.1);
    hDoubleRatio->GetXaxis()->SetLabelSize(0.04);
    hDoubleRatio->GetYaxis()->SetTitle("\\frac{MC Reco/ Data Reco}{UnfoldDataOverAcceptance / MC Gen}");
    hDoubleRatio->GetYaxis()->SetTitleSize(0.035); hDoubleRatio->GetYaxis()->SetTitleOffset(1.5);
    hDoubleRatio->GetYaxis()->SetLabelSize(0.03);

    hDoubleRatio->GetYaxis()->SetRangeUser(0, hDoubleRatio->GetMaximum() * 1.4);

    hDoubleRatio->Draw("E");

    TF1* fitDoubleRatio;
    double a2_dr, a2e_dr, a4_dr, a4e_dr;
    std::tie(fitDoubleRatio, a2_dr, a2e_dr, a4_dr, a4e_dr) = fitCosineFunctionAndReturn(hDoubleRatio, "double_ratio");

    fitDoubleRatio->SetLineColor(kRed);
    fitDoubleRatio->SetLineWidth(2);
    fitDoubleRatio->SetLineStyle(2);
    fitDoubleRatio->Draw("SAME");

    // Calculate Chi² between hDoubleRatio and fitDoubleRatio
    double chi2_dr = 0;
    int ndf_dr = 0;
    for (int i = 1; i <= hDoubleRatio->GetNbinsX(); i++) {
        double val_hist = hDoubleRatio->GetBinContent(i);
        double bin_center = hDoubleRatio->GetBinCenter(i);
        double val_fit = fitDoubleRatio->Eval(bin_center);
        double err_hist = hDoubleRatio->GetBinError(i);
        if (err_hist > 0 && val_hist != 0) { // Skip bins with zero error or content
            chi2_dr += pow((val_hist - val_fit) / err_hist, 2);
            ndf_dr++;
        }
    }
    ndf_dr -= 3; // Subtract number of fit parameters (C, a2, a4)
    int chi2_dr_int = TMath::Nint(chi2_dr);
    int chi2_over_36_dr_int = TMath::Nint(chi2_dr / 36.0);

//    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(62); latex.SetTextSize(0.045); latex.DrawLatex(0.12, 0.96, "CMS");
    latex.SetTextFont(52); latex.SetTextSize(0.035); latex.DrawLatex(0.18, 0.96, "Work in Progress");
    latex.SetTextFont(42); latex.SetTextSize(0.035); latex.DrawLatex(0.65, 0.97, "PbPb #sqrt{s_{NN}} = 5.02 TeV");

    latex.SetTextFont(62); latex.SetTextAlign(13); latex.SetTextSize(0.030);
    latex.DrawLatex(0.18, 0.85, "p_{T}_{ee} < 1.0 GeV");
    latex.DrawLatex(0.18, 0.80, "|#eta_{e}| < 2.4");
    latex.DrawLatex(0.18, 0.75, "7.0 < M_{ee} < 10.0 GeV");

    TPaveText* pt_dr = new TPaveText(0.65, 0.15, 0.88, 0.35, "NDC");
    pt_dr->SetFillColor(0); pt_dr->SetBorderSize(1);
    pt_dr->AddText(Form("#chi^{2}/ndf = %d/%d", chi2_dr_int, ndf_dr));
    pt_dr->AddText(Form("#chi^{2}/36 = %d", chi2_over_36_dr_int));
    pt_dr->AddText(Form("a_{2} = %.3f #pm %.3f", a2_dr, a2e_dr));
    pt_dr->AddText(Form("a_{4} = %.3f #pm %.3f", a4_dr, a4e_dr));
    pt_dr->Draw();

    TLegend* leg3 = new TLegend(0.65, 0.75, 0.85, 0.85);
    leg3->SetBorderSize(0); leg3->SetFillStyle(0);
    leg3->AddEntry(hDoubleRatio, "Double Ratio", "lep");
    leg3->AddEntry(fitDoubleRatio, "Fit: C(1+a_{2}cos2#Delta#phi+a_{4}cos4#Delta#phi)", "l");
    leg3->Draw();

    c3->SaveAs("DoubleRatioComparison_WithChi2.pdf");
    c3->SaveAs("DoubleRatioComparison_WithChi2.png");

    // Cleanup
    delete hDataNorm; delete hMCNorm;
    delete hUnfoldNorm; delete hGenMCNorm;
    delete hRecoRatio; delete hGenRatio; delete hDoubleRatio;
    delete c1; delete c2; delete c3;
    delete leg1; delete leg2; delete leg3;
    delete pt_reco; delete pt_gen; delete pt_dr;
    delete fitReco; delete fitGen; delete fitDoubleRatio;

    file->Close();

    std::cout << "\n✅ Done! Plots show (2π/N)·dN/dφ with cosine fits and double ratio." << std::endl;
}
