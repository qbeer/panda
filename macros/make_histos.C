void make_histos(TString prefix)
{

    std::cout << prefix + "_output_ana.root" << std::endl;

    TFile *in = TFile::Open(prefix + "_output_ana.root");

    if (!in)
    {
        std::cout << "Couldn't open file!" << std::endl;
        return;
    }

    TH1F *hjpsim_all = new TH1F("hjpsim_all", "J/#psi mass (all)", 200, 0, 4.5);
    TH1F *hpsim_all = new TH1F("hpsim_all", "#psi(2S) mass (all)", 200, 0, 5);

    TH1F *hjpsim_fit_all = new TH1F("hjpsim_fit_all", "J/#psi fit mass (all)", 200, 0, 4.5);
    TH1F *hpsim_fit_all = new TH1F("hpsim_fit_all", "#psi(2S) fit mass (all)", 200, 0, 5);

    TH1F *hjpsim_fit_cut_all = new TH1F("hjpsim_fit_cut_all", "J/#psi fit with cut mass (all)", 200, 0, 4.5);
    TH1F *hpsim_fit_cut_all = new TH1F("hpsim_fit_cut_all", "#psi(2S) fit with cut mass (all)", 200, 0, 5);

    TH1F *fitchiq_all = new TH1F("fitchiq_all", "fit chi square all", 200, 0, 4.5);
    TH1F *fitchiq_flagged = new TH1F("fitchiq_flagged", "#fit chi square flagged", 200, 0, 4.5);

    std::cout << "setting up reader" << std::endl;
    TTreeReader npbarp("npbarp", in);

    TTreeReaderValue<float> jpsiMass(npbarp, "pbarp_d0d0m");
    TTreeReaderValue<float> psiMass(npbarp, "pbarp_d0m");

    TTreeReaderValue<float> jpsiFitMass(npbarp, "pbarpfit_d0d0m");
    TTreeReaderValue<float> psiFitMass(npbarp, "pbarpfit_d0m");

    TTreeReaderValue<float> fitchiq(npbarp, "fitchiq");
    TTreeReaderValue<float> fitflag(npbarp, "fitflag");

    std::cout << "setup reader values" << std::endl;

    while (npbarp.Next())
    {
        hjpsim_all->Fill(*jpsiMass);
        hjpsim_fit_all->Fill(*jpsiFitMass);

        hpsim_all->Fill(*psiMass);
        hpsim_fit_all->Fill(*psiFitMass);

        fitchiq_all->Fill(*fitchiq);

        if (*fitflag == 1 && *fitchiq < 100 && *fitchiq > 0)
        {
            hjpsim_fit_cut_all->Fill(*jpsiFitMass);
            hpsim_fit_cut_all->Fill(*psiFitMass);
        }

        if(*fitflag == 1) {
            fitchiq_flagged->Fill(*fitchiq);
        }
    }

    TCanvas *c1 = new TCanvas("c1", "mycanv1", 0, 0, 700, 1400);
    c1->Divide(1, 4);

    c1->cd(1);
    hjpsim_all->Draw();
    hjpsim_fit_all->SetFillColor(42);
    hjpsim_fit_all->Draw("same");

    c1->cd(2);
    hpsim_all->Draw();
    hpsim_fit_all->SetFillColor(42);
    hpsim_fit_all->Draw("same");

    c1->cd(3);
    hjpsim_fit_cut_all->Draw();

    c1->cd(4);
    hpsim_fit_cut_all->Draw();

    TCanvas *c2 = new TCanvas("c2", "mycanv2", 0, 0, 400, 400);
    TCanvas *c3 = new TCanvas("c3", "mycanv3", 0, 0, 400, 400);    
    
    c2->cd();
    fitchiq_all->Draw();

    c3->cd();
    fitchiq_flagged->Draw();

    TCanvas *c4 = new TCanvas("c4", "mycanv4", 0, 0, 400, 400);
    c4->cd();

    RooRealVar mass_peak("m0","m0",90, 60,120);
    RooRealVar width("width","width",10, 0,100);
    RooRealVar a0("a0","",0,1000);
    RooRealVar x("x","",40,140);

    x.setBins(1000);

    RooPolynomial pol5("pol5","",x,RooArgList(a0));//construct a polynomial

    RooRealVar coe_bw("coe_bw","",1,0,1);
    RooRealVar coe_pol("coe_pol","",.0,0,1);

    RooBreitWigner bw("bw","Breit-Wigner function fit",x,mass_peak,width);//construct a Breit-Wigner function
    
    RooDataHist *mass_data=new RooDataHist("JPsiMass","J/Psi mass",x,hjpsim_fit_cut_all);//convert data to roodatahist
   
    RooPlot *plot_frame=x.frame(RooFit::Title("signal invariant mass of muon+ muon- and photon"),RooFit::Bins(1000));
    
    RooAddPdf bwandpol("bwandpol","",bw,pol5,coe_bw);//construct function include background and peak
    bwandpol.fitTo(*mass_data);//fit the curve
    
    mass_data->plotOn(plot_frame);
    bwandpol.plotOn(plot_frame);
    
    plot_frame->Draw();

}
