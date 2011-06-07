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

// Compares the performance of methods based on wether they use gradient info.
void grad(int npar=100) {
    TH2F *frame = new TH2F("frame","",1,0,15,1,0.5,4.5);
    frame->SetXTitle("log10(Accuracy)");
    frame->SetYTitle("log10(Cost)");
    frame->Draw();
    lt->SetMarkerColor(kRed);
    lt->Draw("ncall:accuracy","method<=5","same");
    lt->SetMarkerColor(kGreen);
    lt->Draw("ncall:accuracy","method>5","same");
    lt->SetMarkerColor(kBlue);
    lt->Draw(Form("log10(pow(10,ncall)+2*%d*pow(10,ngrad)):accuracy",npar),"method>5","same");
}

void test() {
    // Initialize graphics options
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(1);
    TCanvas *c = new TCanvas("test","likelytest analysis",800,600);
    c->SetGridx();
    c->SetGridy();
}