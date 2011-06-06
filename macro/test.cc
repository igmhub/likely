// Created 6-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// ROOT macros to analyze the output of the likelytest program.

void slicestats(TH2F *hist) {
    // Lookup the 2D histogram's x binning info
    Int_t nx = hist->GetNbinsX();
    Double_t xmin = hist->GetXaxis()->GetXmin();
    Double_t xmax = hist->GetXaxis()->GetXmax();
    // Create a new histogram of x-slice mean and rms values
    TH1F *mean = new TH1F(Form("%s_m",hist->GetName()),"x-slice mean values",nx,xmin,xmax);
    TH1F *rms = new TH1F(Form("%s_r",hist->GetName()),"x-slice rms values",nx,xmin,xmax);
    // Loop over x bins
    TH1D *slice;
    for(Int_t bin = 1; bin <= nx; bin++) {
        slice = hist->ProjectionY(Form("%s_slice",hist->GetName()),bin,bin);
        mean->SetBinContent(bin, slice->GetMean());
        mean->SetBinError(bin, slice->GetMeanError());
        rms->SetBinContent(bin, slice->GetRMS());
        rms->SetBinError(bin, slice->GetRMSError());
    }
}

void test() {
    // Initialize graphics options
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(1);
}