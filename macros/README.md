# pbar-p system project

* I set up macros to run the simulation, digitization, reconstruction and particle identification before the analysis

* for the analysis I created a macro that combines muons to J/psi and J/spi with poins to Psi(2S) and the Psi(2S) with pions to the pbar - p system

* my final goal would be to anylize pbar - A systems with a heavy ion, preferably Xe (noble gas)

* my supervisor, György Wolf, wrote an article about the mass shift of deep charmonium states (J/psi, psi(3768), psi(3770)) have a mass shift that might be mesured in dileptons' spectra with PANDA

* his simulation code could be an event generator, inclusion in PANDAroot is expected in the future, so far the final state output is in development by him, so I'd be able to map that to a FairASCII file which I could feed the simulation software with

* I generated 10'000 pbar - p events in the meantime with PANDAroot and anylized that with my macro using QATools

```cpp
class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;
class RhoTuple;

void ana_ntp_qa_complete_aolar(int nevts = 5000)
{
  TString parAsciiFile = "all.par";
  TString prefix = "ppbar_psi2s2pi_jpsi2pi_2mu_5000";
  TString input = "psi2spipi_Jpsipipi_mumu.dec";
  TString output = "ana";
  TString friend1 = "pid";

  PndMasterRunAna *fRun = new PndMasterRunAna();
  fRun->SetInput(input);
  fRun->SetOutput(output);
  fRun->AddFriend(friend1);
  fRun->SetParamAsciiFile(parAsciiFile);
  fRun->Setup(prefix);

  gStyle->SetOptFit(1011);

  fRun->Init();

  TFile *out = TFile::Open(prefix + "_output_ana.root", "RECREATE");

  RhoTuple *npbarp = new RhoTuple("npbarp", "pbarp System");

  PndAnalysis *theAnalysis = new PndAnalysis();
  if (nevts == 0) nevts = theAnalysis->GetEntries();

  RhoCandList pbarp, muplus, muminus, piplus, piminus, jpsi, psi2s, jpsiFit, psi2sFit;

  double m0_jpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass(); // Get nominal PDG mass of the J/psi
  RhoMassParticleSelector *jpsiMassSel = new RhoMassParticleSelector("jpsi", m0_jpsi, 1.0);

  TString pidSelection = "PidAlgoEmcBayes;PidAlgoDrc;PidAlgoDisc;PidAlgoStt;PidAlgoMdtHardCuts";

  // changed it accordingly to sim_complete beam momentum, inital lorentz vector of pbarpSystem
  // (p_x, p_y, p_z, E/c) style
  TLorentzVector ini1(0, 0, 9.231552, 9.279084);
  TLorentzVector ini2(0,0,0,0.938);
  TLorentzVector ini(ini1 + ini2);

  PndRhoTupleQA qa(theAnalysis, ini.P());

  int i = 0;
  while (theAnalysis->GetEvent() && i++ < nevts)
  {
    
    if(i%100 == 0) std::cout << "EVENT: #" << i << std::endl;

    theAnalysis->FillList(muplus, "MuonAllPlus", pidSelection);
    theAnalysis->FillList(muminus, "MuonAllMinus", pidSelection);
    theAnalysis->FillList(piplus, "PionAllPlus", pidSelection);
    theAnalysis->FillList(piminus, "PionAllMinus", pidSelection);

    jpsi.Combine(muplus, muminus);
    jpsi.SetType(443); // checked - j/psi

    for (int j = 0; j < jpsi.GetLength(); ++j)
    {
      
      RhoKinVtxFitter vtxfitter(jpsi[j]);
      vtxfitter.Fit();
      RhoCandidate* fitvtx_jpsi = jpsi[j]->GetFit();
      
      RhoKinFitter mfitter(fitvtx_jpsi);
			mfitter.AddMassConstraint(m0_jpsi);
			mfitter.Fit();
      fitvtx_jpsi = fitvtx_jpsi->GetFit();
      
      if (fitvtx_jpsi) jpsiFit.Add(fitvtx_jpsi);

    }

    psi2s.Combine(jpsiFit, piplus, piminus);
    psi2s.SetType(100443); // checked - psi(2s)

    for (int j = 0; j < psi2s.GetLength(); ++j)
    {

      RhoKinVtxFitter vtxfitter(psi2s[j]); 
      vtxfitter.Fit();
      RhoCandidate* fitted_psi2s = psi2s[j]->GetFit();
      
      if(fitted_psi2s) psi2sFit.Add(fitted_psi2s);

    }

    pbarp.Combine(psi2sFit, piplus, piminus);
    pbarp.SetType(88888); // checked - pbarpSystem

    for (int j = 0; j < pbarp.GetLength(); ++j){
      
      qa.qaComp("pbarp_", pbarp[j], npbarp);

      RhoDecayTreeFitter dtf(pbarp[j],ini);
      dtf.setMassConstraint(pbarp[j]->Daughter(0)->Daughter(0));
      
      int check = dtf.Fit();

      npbarp->Column("fitflag",(Float_t)check, -999.9f);
      npbarp->Column("fitchiq",(Float_t)dtf.chiSquare(), -999.9f);
      qa.qaComp("pbarpfit_", pbarp[j]->GetFit(), npbarp);
      
      npbarp->DumpData();
    }

    psi2sFit.Cleanup();
    jpsiFit.Cleanup();

  } // end while loop

  out->cd();

  npbarp->GetInternalTree()->Write();

  out->Save();

}

```

* basically the reaction chain is the following:
        
        pbar -p --> pi+, pi-
                    Psi(2S)  --> pi+, pi-
                                 J/Psi    --> mu+, mu-
                            
* so in th end, I have 6 charged particles, the charmonium states decay ON SPOT, so I have to find a vertex in order to analyze the data correctly, this is the task of the `RhoDecayTreeFitter`, in order to acquire more precise data I set up `RhoKinVtxFitter` for the Psi(2S) and pbar - p system reconstruction as well: `jpsiFit, psi2sFit`
  
* for historgram fitting I wrote a macro using `RooFit` tools and my analysis output:

```cpp
void make_histos_fit(TString prefix)
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

        if (*fitflag == 1)
        {
            hjpsim_fit_cut_all->Fill(*jpsiFitMass);
            hpsim_fit_cut_all->Fill(*psiFitMass);
            fitchiq_flagged->Fill(*fitchiq);
        }
    }

    TCanvas *c1 = new TCanvas("c1", "mycanv1", 0, 0, 700, 700);
    c1->Divide(1, 2);

    c1->cd(1);
    hjpsim_all->Draw();
    hjpsim_fit_all->SetFillColor(42);
    hjpsim_fit_all->Draw("same");

    c1->cd(2);
    hpsim_all->Draw();
    hpsim_fit_all->SetFillColor(42);
    hpsim_fit_all->Draw("same");

    TCanvas *c1_5 = new TCanvas("c1.5", "mycanv1.5", 0,0, 700,700);
    c1_5->Divide(1,2);

    c1_5->cd(1);
    hjpsim_fit_cut_all->Draw();

    c1_5->cd(2);
    hpsim_fit_cut_all->Draw();

    TCanvas *c2 = new TCanvas("c2", "mycanv2", 0, 0, 400, 400);
    TCanvas *c3 = new TCanvas("c3", "mycanv3", 0, 0, 400, 400);

    c2->cd();
    fitchiq_all->Draw();

    c3->cd();
    fitchiq_flagged->Draw();

    TCanvas *c4 = new TCanvas("c4", "mycanv4", 0, 0, 400, 400);
    c4->cd();

    RooRealVar mass_peak("m0", "m0", 3.096, 3.09, 3.1);
    RooRealVar width("width", "width", 0.021, 0, 0.055);
    RooRealVar x("x", "", 0, 15);

    x.setBins(1500);

    RooRealVar coe_bw("coe_bw", "", 4, 3.097, 0.01);

    RooBreitWigner bw("bw", "Breit-Wigner function fit", x, mass_peak, width);

    RooDataHist *mass_data = new RooDataHist("JPsiMass", "J/Psi mass", x, hjpsim_fit_all);

    RooPlot *plot_frame = x.frame(RooFit::Title("J/Psi mass"), RooFit::Bins(1500));

    RooAddPdf bwPlot("bw", "", bw, coe_bw);
    RooFitResult *r = bwPlot.fitTo(*mass_data);
    //r->Print(); // for some reason is null, even with Save()

    mass_data->plotOn(plot_frame);
    bwPlot.plotOn(plot_frame);

    plot_frame->Draw();

    auto legend = new TLegend(0.1, 0.7, 0.48, 0.9);
    legend->SetHeader("J/Psi mass fit", "C"); // option "C" allows to center the header
    legend->AddEntry(hjpsim_fit_all, "J/Psi mass distro - DecayTree fitted", "f");
    legend->AddEntry("bwPlot", "Breit Wigner", "l");
    legend->Draw("same");
}
```

* this creates several canvases which can be saved interactively
  
* for fitting the J/Psi spectrum I use a Breit-Wigner function with E = 3.096 GeV for the J/Psi mass and a very thin width

* I couldn't get the fit parameters yet