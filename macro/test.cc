// Created 6-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// ROOT macros to analyze the output of the likelytest program.

/*
Run benchmarks from ../build using:
./likelytest --npar 100 --ntrial 50 --rho 0.5 --alpha -0.1 > ! test100r.dat
./likelytest --npar 100 --ntrial 50 --rho 0.5 --alpha -1 > ! test100ra.dat
./likelytest --npar 10 --ntrial 50 --rho 0.5 --alpha -0.1 > ! test10r.dat
./likelytest --npar 10 --ntrial 50 --rho 0.5 --alpha -1 > ! test10ra.dat
text2root -i test100r.dat -o '../macro/test100r.root[lt]'
*/

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

// Overall comparison of method costs.
void cmp(int npar=100) {
    // Draw a frame with axis labels.
    TH2F *frame = new TH2F("frame","",1,0,14,1,1,5);
    frame->SetXTitle("log10(Accuracy)");
    frame->SetYTitle("log10(Cost)");
    frame->Draw();
    // Add a legend in the bottom right corner.
    TLegend *legend = new TLegend(0.7,0.12,0.88,0.37);
    legend->SetTextSize(0.016);
    legend->SetFillColor(kWhite);
    // Plot results in red for methods that do not use a gradient calculator.
    lt->SetMarkerColor(2);
    lt->SetMarkerStyle(2);
    lt->Draw("ncall:accuracy>>h1","method==1","same");
    legend->AddEntry(h1,"gsl::simplex2","P");
    lt->SetMarkerStyle(4);
    lt->Draw("ncall:accuracy>>h2","method==2","same");
    legend->AddEntry(h2,"gsl::simplex2rand","P");
    lt->SetMarkerColor(6);
    lt->SetMarkerStyle(5);
    lt->Draw("ncall:accuracy>>h3","method==3","same");
    legend->AddEntry(h3,"mn2::simplex","P");
    lt->SetMarkerStyle(25);
    lt->Draw("ncall:accuracy>>h4","method==4","same");
    legend->AddEntry(h4,"mn2::vmetric","P");
    lt->SetMarkerStyle(26);
    lt->Draw("ncall:accuracy>>h5","method==5","same");
    legend->AddEntry(h5,"mn2::vmetric_fast","P");
    lt->SetMarkerColor(kGreen);
    lt->SetMarkerStyle(27);
    lt->Draw("ncall:accuracy>>h6","method==6","same");
    legend->AddEntry(h6,"mc::walkabout","P");
    // Plot results in blue for methods that use a gradient calculator, with a combined
    // cost equal to ncall + 2*npar*ngrad.
    TString what;
    lt->SetMarkerColor(4);
    lt->SetMarkerStyle(2);
    what = Form("log10(pow(10,ncall)+2*%d*pow(10,ngrad)):accuracy>>h%d",npar,11);
    lt->Draw(what,"method==11","same");
    legend->AddEntry(h11,"gsl::conjugate_fr","P");
    lt->SetMarkerStyle(4);
    what = Form("log10(pow(10,ncall)+2*%d*pow(10,ngrad)):accuracy>>h%d",npar,12);
    lt->Draw(what,"method==12","same");
    legend->AddEntry(h12,"gsl::conjugate_pr","P");
    lt->SetMarkerStyle(5);
    what = Form("log10(pow(10,ncall)+2*%d*pow(10,ngrad)):accuracy>>h%d",npar,13);
    lt->Draw(what,"method==13","same");
    legend->AddEntry(h13,"gsl::vector_bfgs2","P");
    lt->SetMarkerStyle(25);
    what = Form("log10(pow(10,ncall)+2*%d*pow(10,ngrad)):accuracy>>h%d",npar,14);
    lt->Draw(what,"method==14","same");
    legend->AddEntry(h14,"gsl::steepest_descent","P");
    lt->SetMarkerColor(7);
    lt->SetMarkerStyle(26);
    what = Form("log10(pow(10,ncall)+2*%d*pow(10,ngrad)):accuracy>>h%d",npar,15);
    lt->Draw(what,"method==15","same");
    legend->AddEntry(h15,"mn2::vmetric_grad","P");
    lt->SetMarkerStyle(27);
    what = Form("log10(pow(10,ncall)+2*%d*pow(10,ngrad)):accuracy>>h%d",npar,16);
    lt->Draw(what,"method==16","same");
    legend->AddEntry(h16,"mn2::vmetric_grad_fast","P");
    legend->Draw();
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