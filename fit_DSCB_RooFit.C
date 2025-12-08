#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooCrystalBall.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooFit.h"
#include "TSystem.h"
#include <TROOT.h>

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>

using namespace RooFit;

void fit_DSCB_RooFit()
{
    gStyle->SetOptStat(1); // no stats box
    gROOT->SetBatch(kTRUE);

    // ----------------------------------------------------------
    // 1. Load your signal histogram
    // ----------------------------------------------------------
    TFile *f = TFile::Open("/eos/user/h/hsiaoche/Signal/CstarToGJ_M1000_f0p1_13TeV_NANOAOD_v2/CstarToGJ.root");
    if (!f) { std::cerr << "File not found!" << std::endl; return; }

    TH1F *h = (TH1F*)f->Get("hM_reco_selected");
    if (!h) { std::cerr << "Histogram not found!" << std::endl; return; }

    // rms and mean for initial guesses
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();
    double peak = h->GetBinCenter(h->GetMaximumBin());

    // ----------------------------------------------------------
    // 2. Define RooFit observable (mass)
    // ----------------------------------------------------------
    RooRealVar x("x", "M_{#gamma+jet} [GeV]", xmin, xmax);

    // Convert TH1F to RooDataHist
    RooDataHist data("data", "dataset", x, h);

    // ----------------------------------------------------------
    // 3. Define DSCB parameters
    // ----------------------------------------------------------
    RooRealVar x0("x0", "mean", peak, peak-100, peak+100);

    RooRealVar sigmaL("sigmaL", "sigmaL", 50, 1, 200);
    RooRealVar sigmaR("sigmaR", "sigmaR", 50, 1, 200);

    RooRealVar alphaL("alphaL", "alphaL", 1.5, 0.1, 5.0);
    RooRealVar nL("nL", "nL", 5, 0.1, 20);

    RooRealVar alphaR("alphaR", "alphaR", 1.5, 0.1, 5.0);
    RooRealVar nR("nR", "nR", 5, 0.1, 20);

    // ----------------------------------------------------------
    // 4. Build the Double-Sided Crystal Ball PDF
    // ----------------------------------------------------------
    RooCrystalBall dscb("dscb", "Double-Sided Crystal Ball",
            x, x0, sigmaL, sigmaR,
            alphaL, nL, alphaR, nR);

    // ----------------------------------------------------------
    // 5. Fit the model to the data
    // ----------------------------------------------------------
    x.setRange("fitRange", 400, 2000);
    dscb.fitTo(data, SumW2Error(kTRUE), PrintLevel(-1), Range("fitRange"));



    // ----------------------------------------------------------
    // 6. Plot the fit result
    // ----------------------------------------------------------
    TCanvas *c = new TCanvas("c", "DSCB Fit", 800, 600);
    RooPlot *frame = x.frame(Bins(60), Title("DSCB fit to C* signal"));
    data.plotOn(frame, Name("data"),
            MarkerStyle(20),
            MarkerSize(0.8),
            MarkerColor(kBlack),
            LineColor(kBlack));

    dscb.plotOn(frame, Name("dscb"),
            LineColor(kRed),
            LineWidth(2));
    frame->GetXaxis()->SetTitle("M_{#gamma+jet} [GeV]");
    frame->GetYaxis()->SetTitle("Events");
    frame->Draw();


    double chi2 = frame->chiSquare("dscb", "data");
    std::cout << "\nChi2/Ndof = " << chi2 << std::endl;


    TLegend *leg = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(frame->getObject(0), "Signal (MC)", "lep"); // data points
    leg->AddEntry(frame->getObject(1), "DSCB fit",    "l");   // fit curve
    leg->Draw();

    TPaveText *pt = new TPaveText(0.6, 0.45, 0.88, 0.73, "NDC");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.03);

    pt->AddText(Form("Mean = %.1f #pm %.1f",  x0.getVal(),      x0.getError()));
    pt->AddText(Form("#sigma_{L} = %.1f #pm %.1f", sigmaL.getVal(), sigmaL.getError()));
    pt->AddText(Form("#sigma_{R} = %.1f #pm %.1f", sigmaR.getVal(), sigmaR.getError()));
    pt->AddText(Form("#alpha_{L} = %.2f #pm %.2f", alphaL.getVal(), alphaL.getError()));
    pt->AddText(Form("n_{L} = %.2f #pm %.2f",      nL.getVal(),     nL.getError()));
    pt->AddText(Form("#alpha_{R} = %.2f #pm %.2f", alphaR.getVal(), alphaR.getError()));
    pt->AddText(Form("n_{R} = %.2f #pm %.2f",      nR.getVal(),     nR.getError()));
    pt->Draw();

    TPaveText *pt2 = new TPaveText(0.6, 0.38, 0.88, 0.45, "NDC");
    pt2->SetBorderSize(0);
    pt2->SetFillColor(0);
    pt2->SetTextAlign(12);
    pt2->SetTextSize(0.03);
    pt2->AddText(Form("#chi^{2}/N_{dof} = %.6f", chi2));
    pt2->Draw();

    c->SaveAs("DSCB_fit_RooFit.png");


    // ----------------------------------------------------------
    // 7. Print parameters
    // ----------------------------------------------------------
    std::cout << "\n===== DSCB Fit Parameters =====\n";
    x0.Print();
    sigmaL.Print();
    sigmaR.Print();
    alphaL.Print();
    nL.Print();
    alphaR.Print();
    nR.Print();



    // ----------------------------------------------------------
    // 8. Save everything into a RooWorkspace for Combine
    // ----------------------------------------------------------
    RooWorkspace ws("ws", "workspace");
    ws.import(x);
    ws.import(dscb);
    ws.import(data);

    ws.writeToFile("signal_DSCB_workspace.root");

    std::cout << "\nWorkspace saved to: signal_DSCB_workspace.root\n";
    f->Close();
    }
