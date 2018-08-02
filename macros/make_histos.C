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
    //TTreeReader npbarp("npbarp", in);
    TTree* ntp = (TTree*)in->Get("npbarp");

    // TTreeReaderValue<float> jpsiMass(npbarp, "pbarp_d0d0m");
    // TTreeReaderValue<float> psiMass(npbarp, "pbarp_d0m");

    // TTreeReaderValue<float> jpsiFitMass(npbarp, "pbarpfit_d0d0m");
    // TTreeReaderValue<float> psiFitMass(npbarp, "pbarpfit_d0m");

    // TTreeReaderValue<float> fitchiq(npbarp, "fitchiq");
    // TTreeReaderValue<float> fitflag(npbarp, "fitflag");

    // std::cout << "setup reader values" << std::endl;

    // while (npbarp.Next())
    // {
    //     hjpsim_all->Fill(*jpsiMass);
    //     hjpsim_fit_all->Fill(*jpsiFitMass);

    //     hpsim_all->Fill(*psiMass);
    //     hpsim_fit_all->Fill(*psiFitMass);

    //     fitchiq_all->Fill(*fitchiq);

    //     if (*fitflag == 1 && *fitchiq < 100 && *fitchiq > 0)
    //     {
    //         hjpsim_fit_cut_all->Fill(*jpsiFitMass);
    //         hpsim_fit_cut_all->Fill(*psiFitMass);
    //     }

    //     if(*fitflag == 1) {
    //         fitchiq_flagged->Fill(*fitchiq);
    //     }
    // }

    // TCanvas *c1 = new TCanvas("c1", "mycanv", 0, 0, 700, 1400);
    // c1->Divide(1, 4);

    // c1->cd(1);
    // hjpsim_all->Draw();
    // hjpsim_fit_all->SetFillColor(42);
    // hjpsim_fit_all->Draw("same");

    // c1->cd(2);
    // hpsim_all->Draw();
    // hpsim_fit_all->SetFillColor(42);
    // hpsim_fit_all->Draw("same");

    // c1->cd(3);
    // hjpsim_fit_cut_all->Draw();

    // c1->cd(4);
    // hpsim_fit_cut_all->Draw();

    TCanvas *c2 = new TCanvas("c2", "mycanv", 0, 0, 700, 700);
    c2->Divide(1, 2);

    c2->cd(1);
    ntp->Draw("fitchiq");
    //fitchiq_all->Draw();

    c2->cd(2);
    ntp->Draw("fitchiq","fitflag == 1");
    //fitchiq_flagged->Draw();

}

