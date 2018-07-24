class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;

// *** routine to only keep PID matched candidates in list
int SelectTruePid(PndAnalysis *ana, RhoCandList &l)
{
	int removed = 0;
	
	for (int ii=l.GetLength()-1;ii>=0;--ii)
	{
		if ( !(ana->McTruthMatch(l[ii])) )
		{
			l.Remove(l[ii]);
			removed++;
		}
	}
	
	return removed;
}


void tut_ana_pid(int nevts = 0, TString prefix = "signal")
{
 	// *** some variables
	int i=0,j=0, k=0, l=0;
	gStyle->SetOptFit(1011);
	
	// *** the output file for FairRunAna
	TString OutFile="out_dummy.root";  
					
	// *** the files coming from the simulation
	TString inPidFile  = prefix+"_pid.root";    // this file contains the PndPidCandidates and McTruth
	TString inParFile  = prefix+"_par.root";
	
	// *** PID table with selection thresholds; can be modified by the user
	TString pidParFile = TString(gSystem->Getenv("VMCWORKDIR"))+"/macro/params/all.par";	
	
	// *** initialization
	FairLogger::GetLogger()->SetLogToFile(kFALSE);
	FairRunAna* fRun = new FairRunAna();
	FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
	fRun->SetInputFile(inPidFile);
	
	// *** setup parameter database 	
	FairParRootFileIo* parIO = new FairParRootFileIo();
	parIO->open(inParFile);
	FairParAsciiFileIo* parIOPid = new FairParAsciiFileIo();
	parIOPid->open(pidParFile.Data(),"in");
	
	rtdb->setFirstInput(parIO);
	rtdb->setSecondInput(parIOPid);
	rtdb->setOutput(parIO);  
	
	fRun->SetOutputFile(OutFile);
	fRun->Init(); 
	
	// *** create an output file for all histograms
	TFile *out = TFile::Open(prefix+"_ana_pid.root","RECREATE");
	
	// *** create some histograms
	TH1F *hjpsim_all = new TH1F("hjpsim_all","J/#psi mass (all)",200,0,4.5);
	TH1F *hpsim_all  = new TH1F("hpsim_all","#psi(2S) mass (all)",200,0,5);
	
	TH1F *hjpsim_lpid = new TH1F("hjpsim_lpid","J/#psi mass (loose pid)",200,0,4.5);
	TH1F *hpsim_lpid  = new TH1F("hpsim_lpid","#psi(2S) mass (loose pid)",200,0,5);
	
	TH1F *hjpsim_tpid = new TH1F("hjpsim_tpid","J/#psi mass (tight pid)",200,0,4.5);
	TH1F *hpsim_tpid  = new TH1F("hpsim_tpid","#psi(2S) mass (tight pid)",200,0,5);
	
	TH1F *hjpsim_trpid = new TH1F("hjpsim_trpid","J/#psi mass (true pid)",200,0,4.5);
	TH1F *hpsim_trpid  = new TH1F("hpsim_trpid","#psi(2S) mass (true pid)",200,0,5);
	
	
	//
	// Now the analysis stuff comes...
	//
	
	
	// *** the data reader object
	PndAnalysis* theAnalysis = new PndAnalysis();
	if (nevts==0) nevts= theAnalysis->GetEntries();
	
	// *** RhoCandLists for the analysis
	RhoCandList muplus, muminus, piplus, piminus, jpsi, psi2s;
	
	// *** Mass selector for the jpsi cands
	RhoMassParticleSelector *jpsiMassSel=new RhoMassParticleSelector("jpsi",3.096,1.0);
	
	// ***
	// the event loop
	// ***
	while (theAnalysis->GetEvent() && i++<nevts)
	{
		if ((i%100)==0) cout<<"evt " << i << endl;
				
		// ***
		// *** NO PID combinatorics
		// ***
		
		// *** Select with no PID info ('All'); type and mass are set 		
		theAnalysis->FillList(muplus,  "MuonAllPlus");
		theAnalysis->FillList(muminus, "MuonAllMinus");
		theAnalysis->FillList(piplus,  "PionAllPlus");
		theAnalysis->FillList(piminus, "PionAllMinus");
		
		// *** combinatorics for J/psi -> mu+ mu-
		jpsi.Combine(muplus, muminus);
		
		for (j=0;j<jpsi.GetLength();++j) hjpsim_all->Fill( jpsi[j]->M() );
		
		// *** some rough mass selection
		jpsi.Select(jpsiMassSel);
		
		// *** combinatorics for psi(2S) -> J/psi pi+ pi-
		psi2s.Combine(jpsi, piplus, piminus);
		for (j=0;j<psi2s.GetLength();++j) hpsim_all->Fill( psi2s[j]->M() );
		
		
		// ***
		// *** TRUE PID combinatorics
		// ***
		
		// *** do MC truth match for PID type
		SelectTruePid(theAnalysis, muplus);
		SelectTruePid(theAnalysis, muminus);
		SelectTruePid(theAnalysis, piplus);
		SelectTruePid(theAnalysis, piminus);
				
		// *** all combinatorics again with true PID
		jpsi.Combine(muplus, muminus);
		for (j=0;j<jpsi.GetLength();++j) hjpsim_trpid->Fill( jpsi[j]->M() );
		jpsi.Select(jpsiMassSel);
		
		psi2s.Combine(jpsi, piplus, piminus);
		for (j=0;j<psi2s.GetLength();++j) hpsim_trpid->Fill( psi2s[j]->M() );
		
		
		// ***
		// *** LOOSE PID combinatorics
		// ***
		
		// *** and again with PidAlgoMvd;PidAlgoStt;PidAlgoDrc and loose selection
		theAnalysis->FillList(muplus,  "MuonLoosePlus",  "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
		theAnalysis->FillList(muminus, "MuonLooseMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
		theAnalysis->FillList(piplus,  "PionLoosePlus",  "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
		theAnalysis->FillList(piminus, "PionLooseMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
		
		jpsi.Combine(muplus, muminus);
		for (j=0;j<jpsi.GetLength();++j) hjpsim_lpid->Fill( jpsi[j]->M() );
		jpsi.Select(jpsiMassSel);
		
		psi2s.Combine(jpsi, piplus, piminus);
		for (j=0;j<psi2s.GetLength();++j) hpsim_lpid->Fill( psi2s[j]->M() );
		
		
		// ***
		// *** TIGHT PID combinatorics
		// ***
		
		// *** and again with PidAlgoMvd;PidAlgoStt and tight selection
		theAnalysis->FillList(muplus,  "MuonTightPlus",  "PidAlgoMdtHardCuts");
		theAnalysis->FillList(muminus, "MuonTightMinus", "PidAlgoMdtHardCuts");
		theAnalysis->FillList(piplus,  "PionLoosePlus",  "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
		theAnalysis->FillList(piminus, "PionLooseMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
		
		jpsi.Combine(muplus, muminus);
		for (j=0;j<jpsi.GetLength();++j) hjpsim_tpid->Fill( jpsi[j]->M() );
		jpsi.Select(jpsiMassSel);
		
		psi2s.Combine(jpsi, piplus, piminus);
		for (j=0;j<psi2s.GetLength();++j) hpsim_tpid->Fill( psi2s[j]->M() );
		
	}
	
	// *** write out all the histos
	out->cd();
	
	hjpsim_all->Write();
	hpsim_all->Write();
	hjpsim_lpid->Write();
	hpsim_lpid->Write();
	hjpsim_tpid->Write();
	hpsim_tpid->Write();
	hjpsim_trpid->Write();
	hpsim_trpid->Write();
			
	out->Save();
	
}
