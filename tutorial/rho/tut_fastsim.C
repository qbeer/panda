// *******
// Tutorial Macro for running fast simulation 
// *******

// The parameters are
// -------------------
// USAGE:\n";
// tut_fastsim.C+( [nevt], [prefix], [decfile], [mom], [res] )
//    [nevt]     : number of events; default = 1000
//    [pref]     : prefix for output
//    [decfile]  : decfile; 'DPM'/'FTF' uses DPM/FTF generator instead
//    [mom]      : pbar momentum
//    [res]      : resonance (ignored when running DPM); default = 'pbarpSystem0'


void tut_fastsim(Int_t nEvents = 1000, TString prefix = "signal", TString Decfile="pp_jpsi2pi_jpsi_mumu.dec", Float_t Mom=6.232, TString Resonance="pbarpSystem" )
{
	// if Mom<0, interprete as -E_cm -> compute mom
	double mp = 0.938272;
	if (Mom<0)
	{
		double X = (Mom*Mom-2*mp*mp)/(2*mp);
		Mom = sqrt(X*X-mp*mp);
	}
	
	// Allow shortcut for resonance
	if (Resonance=="pbp")  Resonance = "pbarpSystem";
	if (Resonance=="pbp0") Resonance = "pbarpSystem0";

	// Prevent generator from throwing a lot of warnings
	TLorentzVector fIni(0,0,Mom,mp+sqrt(Mom*Mom+mp*mp));
	TDatabasePDG::Instance()->AddParticle("pbarpSystem","pbarpSystem",fIni.M(),kFALSE,0.1,0, "",88888,0);
	TDatabasePDG::Instance()->AddParticle("pbarpSystem0","pbarpSystem0",fIni.M(),kFALSE,0.1,0, "",88880,0);
		
	//----- Switches for Simulation Options ------------------------------
	Bool_t enableSplitoff    = true;  // create e.-m. and hadronic split offs
	Bool_t mergeNeutrals     = true;  // merge neutrals (for merged pi0s)
	Bool_t electronBrems     = true;  // bremsstrahlung loss for electrons 
	
	//----- Presist simulation output ------------------------------
	Bool_t persist           = true;  

	//-----General settings-----------------------------------------------
	TString BaseDir =  gSystem->Getenv("VMCWORKDIR");
	TString splitpars = BaseDir+"/fsim/splitpars.dat";
	gRandom->SetSeed();

	//-----User Settings:-------------------------------------------------
	TString OutputFile = prefix+"_fast.root";
	
	gDebug             = 0;

	// choose your event generator
	Bool_t UseEvtGenDirect  = kTRUE;
	Bool_t UseDpm           = kFALSE;

	// use DPM generator; default: inelastic @ pbarmom = mom
	if (Decfile.BeginsWith("DPM"))
	{
		UseEvtGenDirect = kFALSE;
		UseDpm          = kTRUE;
	}

	// Start a stop watch
	TStopwatch timer;
	timer.Start();

	// Create the Simulation run manager
	// --------------------------------
	FairRunSim *fRun = new FairRunSim();
	fRun->SetOutputFile(OutputFile.Data());
	fRun->SetWriteRunInfoFile(kFALSE);

	FairLogger::GetLogger()->SetLogToFile(kFALSE);

	// -------------------------------
	// Create and Set Event Generator
	// -------------------------------
	//FairFilteredPrimaryGenerator* primGen = new FairFilteredPrimaryGenerator();
	FairPrimaryGenerator* primGen = new FairPrimaryGenerator();
	fRun->SetGenerator(primGen);
	fRun->SetName("TGeant3");

	if(UseDpm)
	{
		int mode = 0;
		if (Decfile=="DPM1") mode = 1;
		if (Decfile=="DPM2") mode = 2;
		
		PndDpmDirect *Dpm= new PndDpmDirect(Mom,mode);  // 0 = inelastic, 1 = inelastic & elastic, 2 = elastic
		Dpm->SetUnstable(111);   // pi0
		Dpm->SetUnstable(310);   // K_S0
		Dpm->SetUnstable(3122);  // Lambda
		Dpm->SetUnstable(-3122); // anti-Lambda
		Dpm->SetUnstable(221);   // eta
		primGen->AddGenerator(Dpm);
	}

	if(UseEvtGenDirect)
	{
		PndEvtGenDirect *EvtGen = new PndEvtGenDirect(Resonance, Decfile.Data(), Mom);
		EvtGen->SetStoreTree(kTRUE);
		primGen->AddGenerator(EvtGen);
	}

	// ------------- switch off the transport of particles
	primGen->DoTracking(kFALSE);

	//---------------------Create and Set the Field(s)----------
	PndMultiField *fField= new PndMultiField("AUTO");
	fRun->SetField(fField);

	
	// ***********************************
	// Setup the Fast Simulation Task
	// ***********************************

	PndFastSim* fastSim = new PndFastSim(persist);
		
	// increasing verbosity increases the amount of console output (mainly for debugging)
	fastSim->SetVerbosity(0);

	// enable the merging of neutrals if they have similar direction
	//-----------------------------
	fastSim->MergeNeutralClusters(mergeNeutrals);

	// enable bremsstahlung loss for electrons
	//-----------------------------
	fastSim->EnableElectronBremsstrahlung(electronBrems);

	//enable the producting of parametrized neutral (hadronic) split offs
	// generate electro-magnetic / hadronic split offs in the EMC? switch off when running w/o EMC

	if (enableSplitoff)	fastSim->EnableSplitoffs(splitpars.Data());

	fastSim->SetUseFlatCov(true);
	
	// -----------------------------------------------------------------------------------
	// Tracking: Set up in parts of theta coverage. All modelled by PndFsmSimpleTracker.
	// Mind: Numbers on resolution (pRes,thtRes,phiRes) and efficiency are guessed
	// -----------------------------------------------------------------------------------
	
	// - (Full Panda Tracking: STT MVD GEM FTS)

	fastSim->AddDetector("ScSttAlone",  "thtMin=145.  thtMax=159.5 ptmin=0.1 pmin=0.0 pRes=0.04 thtRes=0.006 phiRes=0.007 efficiency=0.25");
	fastSim->AddDetector("ScSttMvd",    "thtMin=20.9  thtMax=145.  ptmin=0.1 pmin=0.0 pRes=0.02 thtRes=0.001 phiRes=0.001 efficiency=0.85");
	fastSim->AddDetector("ScSttMvdGem", "thtMin=7.8   thtMax=20.9  ptmin=0.1 pmin=0.0 pRes=0.02 thtRes=0.001 phiRes=0.001 efficiency=0.85");
	fastSim->AddDetector("ScMvdGem",    "thtMin=5.    thtMax=7.8   ptmin=0.1 pmin=0.0 pRes=0.03 thtRes=0.001 phiRes=0.001 efficiency=0.60");

	// Fwd spectrometer enabled -> use Fwd tracking system
	fastSim->AddDetector("ScFts",       "thtMin=0.    thtMax=5.    ptmin=0.0 pmin=0.5 pRes=0.05  thtRes=0.002 phiRes=0.002 efficiency=0.80");

	// -----------------------------------------------------------------------------------
	// Vertexing
	// -----------------------------------------------------------------------------------

	// MVD and GEM are enabled -> better vertexing in central region
	fastSim->AddDetector("ScVtxMvd",   "thtMin=5. thtMax=145. ptmin=0.1 vtxRes=0.005 efficiency=1."); // efficiency=1: all tracks found in trackers will get a vertex information
	fastSim->AddDetector("ScVtxNoMvd", "thtMin=0. thtMax=5.   ptmin=0.0 vtxRes=0.05  efficiency=1."); // efficiency=1: all tracks found in trackers will get a vertex information

	// -----------------------------------------------------------------------------------
	// EM Calorimeters w/ default parameters
	// (don't have to be set, just to list the available parameters
	// -----------------------------------------------------------------------------------

	fastSim->AddDetector("EmcFwCap", "thtMin=10.0 thtMax=22.0 Emin=0.01 dist=2.5");
	fastSim->AddDetector("EmcBwCap", "thtMin=142.0 thtMax=160.0 Emin=0.01 dist=0.7");

	// EmcBarrel also allows to set phiMin and phiMax and can be added multiple times as EmcBarrel1, EmcBarrel2, etc.
	// Should be made constistent with EmcPidBarrel below
	fastSim->AddDetector("EmcBarrel","thtMin=22.0 thtMax=142.0 Emin=0.01 barrelRadius=0.5");

	// Fwd spectrometer enabled -> use Fwd EMC
	fastSim->AddDetector("EmcFS",    "thtMin=0.05 thtMax=10.0 aPar=0.013 bPar=0.0283 Emin=0.01 dist=8.2");

	// -----------------------------------------------------------------------------------
	// PID
	// -----------------------------------------------------------------------------------

	// PID detectors being always in: STT, MUO Barrel, EMC FwdCap, EMC BwdCap
	//Note: A dEdX parametrization from 2008
	fastSim->AddDetector("SttPid","thtMin=7.8 thtMax=159.5 ptmin=0.1 dEdxRes=0.15 efficiency=1.");
	fastSim->AddDetector("ScMdtPidBarrel", "thtMin=10.0 thtMax=130.0 pmin=0.5 efficiency=0.95 misId=0.01");
	fastSim->AddDetector("ScEmcPidFwCap",  "thtMin=10.0  thtMax=22.0  ptmin=0.0 pmin=0.0 efficiency=1.0");
	fastSim->AddDetector("ScEmcPidBwCap",  "thtMin=142.0 thtMax=160.0  ptmin=0.0 pmin=0.0 efficiency=1.0");

	// MVD and GEM are enabled -> MVD PID available
	//Note: A Bethe-Bloch-Landau-Gauss Prametrization from 2008
	fastSim->AddDetector("MvdPid","thtMin=5.  thtMax=133.6 ptmin=0.1  dEdxResMulti=1. efficiency=1.");

	// EMC Barrel enable -> EMC barrel PID available
	fastSim->AddDetector("ScEmcPidBarrel", "thtMin=22.0  thtMax=142.0 ptmin=0.2 pmin=0.0 efficiency=1.0");

	// Barrel DIRC enabled
	fastSim->AddDetector("DrcBarrel","thtMin=22.0 thtMax=140.0 dthtc=0.01 nPhotMin=5 effNPhotons=0.075");

	// Disc DIRC enabled
	fastSim->AddDetector("DrcDisc","thtMin=5.0 thtMax=22.0 dthtc=0.01 nPhotMin=5 effNPhotons=0.075");

	// Fwd spectrometer enabled -> use RICH, FwdMUO and EMC FS
	fastSim->AddDetector("ScEmcPidFS",     "thtMin=0.5   thtMax=10.0  ptmin=0.0 pmin=0.5 efficiency=1.0");
	fastSim->AddDetector("Rich","angleXMax=5.0 angleYMax=10.0 dthtc=0.01 nPhotMin=5 effNPhotons=0.075");
	fastSim->AddDetector("ScMdtPidForward","thtMin=0.0  thtMax=10.0  pmin=0.5 efficiency=0.95 misId=0.01");


	fRun->AddTask(fastSim);
		
	//-------------------------  Initialize the RUN  -----------------
	fRun->Init();
	
	//-------------------------  Run the Simulation  -----------------
	fRun->Run(nEvents);
	
	
	//------------------------Print some info and exit----------------
	timer.Stop();
	Double_t rtime = timer.RealTime();
	Double_t ctime = timer.CpuTime();
	printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);
}

