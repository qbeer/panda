// **********************************************************************************************
// Macro for running fast simulation + Software Trigger + simple analysis (PndSimpleCombinerTask)
// **********************************************************************************************

// The parameters are
// -------------------
// quickfsimana.C( <pref>, <decfile>, <mom>, <decay>, [nevt], [res], [parms], [runST], [runnum], [mode] )
//    <pref>     : output file names prefix
//    <decfile>  : EvtGen decfile; DPM/FTF/BOX uses DPM/FTF generator (inelastic mode) or box generator instead
//    <mom>      : EvtGen, DPM, FTF: pbar momentum (negative values are interpreted as -E_cm); BOX generator: maximum particle momentum
//    [decay]    : the decay pattern to be reconstructed, e.g. 'phi -> K+ K-; D_s+ -> phi pi+'; emtpy string: only fast sim will be run
//    [nevt]     : number of events; default = 1000
//    [res]      : initial resonance or particle type for BOX generator (ignored when running DPM); default = 'pbarpSystem0'
//    [parms]    : parameters for the analysis, e.g. 'mwin=0.4:mwin(phi)=0.1:emin=0.1:pmin=0.1:qamc'; 'qapart' runs PndParticleQATask
//    [runST]    : if 'true' runs Software Trigger (default: false)
//    [runnum]   : integer run number (default: 0)
//    [mode]     : arbitrary mode number (default: 0)
// -------------------

void    getRange(TString par, double &min, double &max);
TString getInitialResonance(TString &fEvtGenFile);

void quickfsimana(TString Prefix="", TString Decfile="", Float_t Mom=0., TString anadecay="", Int_t nEvents = 1000, TString anaparms="", bool runST=false, int run=0 , int runmode=0)
{
	if (Prefix=="" || Decfile=="" || Mom==0. ) 
	{
		cout << "USAGE:\n";
		cout << "quickfsimana.C( <pref>, <decfile>, <mom>, <decay>, [nevt], [parms], [runST], [runnum], [mode] )\n\n";
		cout << "   <pref>     : output file names prefix\n";
		cout << "   <decfile>  : EvtGen decfile 'xxx.dec' or 'xxx.dec:iniRes'; DPM/FTF/BOX uses DPM/FTF generator or BOX generator instead\n";
		cout << "                DPM settings: DPM = inelastic only, DPM1 = inel. + elastic, DPM2 = elastic only\n";
		cout << "                FTF settings: FTF = inelastic only, FTF1 = inel. + elastic\n";
		cout << "                BOX settings: optional ranges 'p/tht/phi[min,max]' separated with colon; single number = fixed value; example: 'BOX:p[1,5]:tht[45]:phi[90,210]'\n";
		cout << "   <mom>      : EvtGen, DPM, FTF: pbar momentum (negative values are interpreted as -E_cm); BOX generator w/o special settings: maximum particle momentum\n";
		cout << "   [decay]    : the decay pattern to be reconstructed, e.g. 'phi -> K+ K-; D_s+ -> phi pi+'; '': only fast sim w/o reco will be run\n";
		cout << "   [nevt]     : number of events; default = 1000\n";
		//cout << "   [res]      : initial resonance or particle type for BOX generator (ignored when running DPM); default = 'pbarpSystem0'\n";
		cout << "   [parms]    : parameters for the analysis, e.g. 'mwin=0.4:mwin(phi)=0.1:emin=0.1:pmin=0.1:qamc'; 'qapart' runs PndParticleQATask: 'persist' saves PndPidCandidates\n";
		cout << "   [runST]    : if 'true' runs Software Trigger (default: false)\n";
		cout << "   [runnum]   : integer run number (default: 0)\n";
		cout << "   [mode]     : arbitrary mode number (default: 0)\n\n";
		cout << "Example 1 - Do reco for EvtGen events : root -l -b -q 'quickfsimana.C(\"jpsi2pi\", \"pp_jpsi2pi_jpsi_mumu.dec\", 6.23, \"J/psi -> mu+ mu-; pbp -> J/psi pi+ pi-\", 1000, \"fit4c:mwin=0.6\")'\n";
		cout << "Example 2 - Particle QA for BOX gen   : root -l -b -q 'quickfsimana.C(\"single_kplus\", \"BOX:type(321,1)\", 8.0, \"\", 1000, \"qapart\")'\n";	
		cout << "Example 3 - Run fast sim only for DPM : root -l -b -q 'quickfsimana.C(\"bkg\", \"DPM\", 6.23, \"\", 1000)'\n\n";	
		return;
	}
	
	// persist fast sim output?
	bool persist = (anadecay == "" && anaparms == "") || anaparms.Contains("persist");
	
	// do some reconstruction ?
	bool doreco  = (anadecay != "" || anaparms.Contains("nevt"));
	
	// do particle QA?
	bool partQA  = (anaparms.Contains("qapart"));
	bool mc      = !(anaparms.Contains("!mc")) && !(anaparms.Contains("qamc"));
	bool neut    = !(anaparms.Contains("!neut"));
	bool chrg    = !(anaparms.Contains("!chrg"));
		
	// for submission to queue all blanks in decay string were replaced by '§'; now we replace again the other way around
	anadecay.ReplaceAll("§"," ");
	
	// if Mom<0, interprete as -E_cm
	double mp = 0.938272;
	if (Mom<0)
	{
		double X = (Mom*Mom-2*mp*mp)/(2*mp);
		Mom = sqrt(X*X-mp*mp);
	}
	
	// Allow shortcut for resonance pbarpSystem
	anadecay.ReplaceAll("pbp", "pbarpSystem");

	// Prevent generator from throwing a lot of warnings
	TLorentzVector fIni(0,0,Mom,mp+sqrt(Mom*Mom+mp*mp));
	TDatabasePDG *pdg=TDatabasePDG::Instance();
	pdg->AddParticle("pbarpSystem","pbarpSystem",fIni.M(),kFALSE,0.1,0, "",88888,0);
	pdg->AddParticle("pbarpSystem0","pbarpSystem0",fIni.M(),kFALSE,0.1,0, "",88880,0);

	//-----Evaluate Detector Setup ---------------------------------------
	bool SwMvdGem  = true;  // Enable MVD and GEM for central tracking in addition to STT
	bool SwEmcBar  = true;  // Enable EMC barrel for calorimetry (neutral detection and PID component)
	bool SwDrc     = true;  // Enable Barrel DIRC for PID
	bool SwDsc     = true;  // Enable Disc DIRC for PID
	bool SwFwdSpec = true;  // Enable complete Forward Spectrometer (= Fwd Spec. EMC, Fwd Tracking, RICH, Fwd MUO)
	
	//----- Switches for Simulation Options ------------------------------
	Bool_t enableSplitoff    = true;   // create e.-m. and hadronic split offs
	Bool_t mergeNeutrals     = true;   // merge neutrals (for merged pi0s)
	Bool_t electronBrems     = true;   // bremsstrahlung loss for electrons 
	Bool_t useEventFilter    = false;  // enable Fast Sim event filter. *** Needs configuration (see below) *** 
	Bool_t usePndEventFilter = false;  // enable Panda event filter.    *** Needs configuration (see below) *** 
	
	//-----General settings-----------------------------------------------
	TString BaseDir =  gSystem->Getenv("VMCWORKDIR");
	TString splitpars = BaseDir+"/fsim/splitpars.dat";
	gRandom->SetSeed();

	//-----User Settings:-------------------------------------------------
	TString  OutputFile     = TString::Format("%s_%d_ana.root",Prefix.Data(), run);
	if (persist) OutputFile = TString::Format("%s_%d_fsim.root",Prefix.Data(), run);
	
	gDebug             = 0;

	// choose your event generator
	Bool_t UseEvtGenDirect  = kTRUE;
	Bool_t UseFtf           = kFALSE;
	Bool_t UseDpm           = kFALSE;
	Bool_t UseBoxGenerator  = kFALSE;

	// use DPM generator; default: inelastic @ pbarmom = mom
	if (Decfile.BeginsWith("DPM") && !Decfile.EndsWith(".dec"))
	{
		UseEvtGenDirect = kFALSE;
		UseDpm 	      = kTRUE;
	}

	// use FTF generator; 
	if (Decfile.BeginsWith("FTF") && !Decfile.EndsWith(".dec"))
	{
		UseEvtGenDirect = kFALSE;
		UseFtf 	        = kTRUE;
	}

	// use BOX generator; defaults
	Double_t BoxMomMin  = 0.05;   // minimum momentum for box generator
	Double_t BoxMomMax  = Mom;    // maximum   "       "
	Double_t BoxThtMin  = 0. ;    // minimum theta for box generator
	Double_t BoxThtMax  = 180.;   // maximum   "       "
	Double_t BoxPhiMin  = 0. ;    // minimum phi for box generator
	Double_t BoxPhiMax  = 360.;   // maximum   "       "
	Bool_t   BoxCosTht  = false;  // isotropic in cos(theta) instead theta
  
	Int_t    BoxType    = 13;     // default particle muon
	Int_t    BoxMult    = 1;      // default particle multiplicity
	Double_t type=0,mult=0;       // ref. parameters for range function	

	if (Decfile.BeginsWith("BOX") && !Decfile.EndsWith(".dec"))
	{
		UseEvtGenDirect = kFALSE;
		UseBoxGenerator = kTRUE;
		Decfile.ToLower();
		
		if (Decfile!="box")
		{
			Decfile.ReplaceAll("box","");
			Decfile.ReplaceAll(" ","");
			Decfile += ":";
			
			while (Decfile.Contains(":"))
			{
				TString curpar = Decfile(0,Decfile.Index(":"));
				Decfile = Decfile(Decfile.Index(":")+1,1000);
				curpar.ReplaceAll("[","("); curpar.ReplaceAll("]",")"); 
				
				if (curpar.BeginsWith("type(")) {getRange(curpar,type,mult); BoxType = (Int_t)type; BoxMult = (Int_t)mult; }
				if (curpar.BeginsWith("p("))     getRange(curpar,BoxMomMin,BoxMomMax);
				if (curpar.BeginsWith("tht("))   getRange(curpar,BoxThtMin,BoxThtMax);
				if (curpar.BeginsWith("ctht(")) {getRange(curpar,BoxThtMin,BoxThtMax); BoxCosTht=true;}
				if (curpar.BeginsWith("phi("))   getRange(curpar,BoxPhiMin,BoxPhiMax);
			}
		}
		
		cout <<"BOX generator range: type["<<BoxType<<","<<BoxMult<<"]  p["<<BoxMomMin<<","<<BoxMomMax<<"]  tht["<<BoxThtMin<<","<<BoxThtMax<<"]"<<(BoxCosTht?"*":"")<<"  phi["<<BoxPhiMin<<","<<BoxPhiMax<<"]"<<endl;
	}


	// Start a stop watch
	TStopwatch timer;
	timer.Start();

	// Create the Simulation run manager
	// --------------------------------
	FairRunSim *fRun = new FairRunSim();
	fRun->SetOutputFile(OutputFile.Data());
	fRun->SetGenerateRunInfo(kFALSE);
	if (!persist) fRun->SetUserConfig(BaseDir+"/tutorials/analysis/g3ConfigNoMC.C");

	FairLogger::GetLogger()->SetLogToFile(kFALSE);


	// Create and Set Event Generator
	// -------------------------------
	FairFilteredPrimaryGenerator* primGen = new FairFilteredPrimaryGenerator();
	if (!usePndEventFilter) primGen->SetVerbose(0);
	//FairPrimaryGenerator* primGen = new FairPrimaryGenerator();
	fRun->SetGenerator(primGen);
	fRun->SetName("TGeant3");
	
	// Box Generator
	if(UseBoxGenerator)
	{  
		PndBoxGenerator* boxGen = new PndBoxGenerator(BoxType, BoxMult);
		boxGen->SetDebug(0);
		
		boxGen->SetPRange(BoxMomMin,BoxMomMax);      // GeV/c
		boxGen->SetPhiRange(BoxPhiMin, BoxPhiMax);   // Azimuth angle range [degree]
		boxGen->SetThetaRange(BoxThtMin, BoxThtMax); // Polar angle in lab system range [degree]
		
		if (BoxCosTht) boxGen->SetCosTheta();
		
		boxGen->SetXYZ(0., 0., 0.); //cm
		primGen->AddGenerator(boxGen);
	}
	
	// DPM Generator
	if(UseDpm)
	{
		int mode = 0;
		if (Decfile=="DPM1") mode = 1;
		if (Decfile=="DPM2") mode = 2;
		
		PndDpmDirect *Dpm= new PndDpmDirect(Mom,mode);  // 0 = inelastic, 1 = inelastic & elastic, 2 = elastic
		// since fastsim doesn't have a transport, let all long-living resonances decay by the generator
		Dpm->SetUnstable(111);   // pi0
		Dpm->SetUnstable(310);   // K_S0
		Dpm->SetUnstable(311);   // K0
		Dpm->SetUnstable(-311);  // K0bar
		Dpm->SetUnstable(3122);  // Lambda0
		Dpm->SetUnstable(-3122); // anti-Lambda0
		Dpm->SetUnstable(221);   // eta*/
		Dpm->SetUnstable(3222);  // Sigma+
		Dpm->SetUnstable(-3222); // anti-Sigma-
		Dpm->SetUnstable(3334);  // Omega-
		primGen->AddGenerator(Dpm);
	}
	
	// FTF Generator
	if(UseFtf)
	{
		bool noelastic = true;
		if (Decfile=="FTF1") noelastic=false;
		PndFtfDirect *Ftf = new PndFtfDirect("anti_proton", "G4_H", 1, "ftfp", Mom, 0, noelastic); 
		primGen->AddGenerator(Ftf);
	}

	// EvtGen Generator
	if(UseEvtGenDirect)
	{
		TString Resonance=getInitialResonance(Decfile);
		Resonance.ReplaceAll("pbp","pbarpSystem");

		PndEvtGenDirect *EvtGen = new PndEvtGenDirect(Resonance, Decfile.Data(), Mom);
		EvtGen->SetStoreTree(kTRUE);
		primGen->AddGenerator(EvtGen);
	}

	// ------------- switch off the transport of particles
	primGen->DoTracking(kFALSE);

	//---------------------Create and Set the Field(s)----------
	PndMultiField *fField= new PndMultiField("AUTO");
	fRun->SetField(fField);

	// Setup the Fast Simulation Task
	//-----------------------------
	PndFastSim* fastSim = new PndFastSim(persist);
		
	// increasing verbosity increases the amount of console output (mainly for debugging)
	fastSim->SetVerbosity(0);

	// set PANDA event filters
	//-----------------------------
	if (usePndEventFilter)
	{
		cout <<"Using FairEventFilter"<<endl;
		primGen->SetFilterMaxTries(100000);
		
		//FairEvtFilterOnSingleParticleCounts* chrgFilter = new FairEvtFilterOnSingleParticleCounts("chrgFilter");
		//chrgFilter->AndMinCharge(4, FairEvtFilter::kCharged);
		//primGen->AndFilter(chrgFilter);
		 		 		
		//FairEvtFilterOnCounts* neutFilter = new FairEvtFilterOnCounts("neutFilter");
		//neutFilter->AndMaxCharge(4, FairEvtFilter::kNeutral);
		//primGen->AndFilter(neutFilter);
		
		PndEvtFilterOnInvMassCounts* eeInv= new PndEvtFilterOnInvMassCounts("eeInvMFilter");
		//eeInv->SetVerbose();//highest commenting level of the FairEvtFilterOnCounts
		eeInv->SetPdgCodesToCombine( 13, -13);
		eeInv->SetMinMaxInvMass( 2.5, 3.3 );
		eeInv->SetMinMaxCounts(1,10000);
 		primGen->AndFilter(eeInv);  //add filter to fFilterList		
	}
	
	// set event filters
	//-----------------------------
	if (useEventFilter)
	{
	      // Filters are:
	      // -----------
	      // fastSim->SetMultFilter(type, min, max); 
	      // requires min <= mult <= max
	      
	      // available types are:
	      
	      //  "+"   : positive charged particles
	      //  "-"   : negative charged particles
	      //  "gam" : gammas
	      //  "pi0" : pi0 candidates ( -> 2 gammas); mass window 0.135 +- 0.03 GeV
	      //  "eta" : eta candidates ( -> 2 gammas); mass window 0.547 +- 0.04 GeV 
	      //  "ks"  : K_S candidates ( -> pi+ pi-);  mass window 0.497 +- 0.04 GeV

 	      //fastSim->SetMultFilter("ks",   1,1000);  // at least 1 KS
	      
// 	      fastSim->SetMultFilter("+",   2,1000);  // at least 2 trk+
// 	      fastSim->SetMultFilter("-",   2,1000);  // at least 2 trk-
// 	      fastSim->SetMultFilter("gam", 0,   4);  // at most 4 gammas

	      // fastSim->SetInvMassFilter(comb, m_min, m_max, mult);
	      
	      // requires at least mult combined candidates with m_min < m < m_max
	      
	      // comb is a TString describing the combinatoric
	      // - particle codes are: e+ e- mu+ mu- pi+ pi- k+ k- p+ p- gam pi0 ks eta
	      // - codes must be separated with a single blank
	      // - for charged final states only the mass is set; no pdg code selection is done! 
	      // - optional a 'cc' added at the end of also takes into account charge conjugation
	      
	      // Examples: 
	      // - ("k+ k-", 0.98, 1.1, 2)       : forms K+ K- candidate and requires >=2 in the given window
	      // - ("ks k+ pi- cc", 2.8, 3.2,1 ) : forms ks k+ pi- / ks k- pi+ cands and req. at least one in window
	      
	      fastSim->SetInvMassFilter("mu+ mu-",2.5,3.3,1);  // look for J/psi -> e+ e- candidate
	}

	// enable the merging of neutrals if they have similar direction
	//-----------------------------
	fastSim->MergeNeutralClusters(mergeNeutrals);

	// enable bremsstahlung loss for electrons
	//-----------------------------
	fastSim->EnableElectronBremsstrahlung(electronBrems);

	//enable the producting of parametrized neutral (hadronic) split offs
	// generate electro-magnetic / hadronic split offs in the EMC? switch off when running w/o EMC

	if (enableSplitoff)
		fastSim->EnableSplitoffs(splitpars.Data());

	fastSim->SetUseFlatCov(true);
	// -----------------------------------------------------------------------------------
	//Tracking: Set up in parts of theta coverage. All modelled by PndFsmSimpleTracker.
	// Mind: Numbers on resolution (pRes,thtRes,phiRes) and efficiency are guessed
	// -----------------------------------------------------------------------------------
	if (SwMvdGem) // MVD and GEM are enabled; combined tracking available
	{
		// - (Full Panda Tracking: STT MVD GEM FTS)

		fastSim->AddDetector("ScSttAlone",  "thtMin=145.  thtMax=159.5 ptmin=0.1 pmin=0.0 pRes=0.04 thtRes=0.006 phiRes=0.007 efficiency=0.25");
		fastSim->AddDetector("ScSttMvd",    "thtMin=20.9  thtMax=145.  ptmin=0.1 pmin=0.0 pRes=0.02 thtRes=0.001 phiRes=0.001 efficiency=0.85");
		fastSim->AddDetector("ScSttMvdGem", "thtMin=7.8   thtMax=20.9  ptmin=0.1 pmin=0.0 pRes=0.02 thtRes=0.001 phiRes=0.001 efficiency=0.85");
		fastSim->AddDetector("ScMvdGem",    "thtMin=5.    thtMax=7.8   ptmin=0.1 pmin=0.0 pRes=0.03 thtRes=0.001 phiRes=0.001 efficiency=0.60");
	}
	else // MVD and GEM are disabled; only STT tracking in central region
	{
		// - STT alone:
		fastSim->AddDetector("ScSttAlone",  "thtMin=133.6 thtMax=159.5 ptmin=0.1 pmin=0.0 pRes=0.04 thtRes=0.006 phiRes=0.007 efficiency=0.25");
		fastSim->AddDetector("ScSttAlone2", "thtMin=20.9  thtMax=133.6 ptmin=0.1 pmin=0.0 pRes=0.04 thtRes=0.006 phiRes=0.007 efficiency=0.80");
		fastSim->AddDetector("ScSttAlone3", "thtMin=7.8   thtMax=20.9  ptmin=0.1 pmin=0.0 pRes=0.04 thtRes=0.006 phiRes=0.007 efficiency=0.25");
	}

	if (SwFwdSpec) // Fwd spectrometer enabled -> use Fwd tracking system
	{
		fastSim->AddDetector("ScFts",       "thtMin=0.    thtMax=5.    ptmin=0.0 pmin=0.5 pRes=0.05  thtRes=0.002 phiRes=0.002 efficiency=0.80");
	}

	// -----------------------------------------------------------------------------------
	// Vertexing
	// -----------------------------------------------------------------------------------
	if (SwMvdGem) // MVD and GEM are enabled -> better vertexing in central region
	{
		fastSim->AddDetector("ScVtxMvd",   "thtMin=5. thtMax=145. ptmin=0.1 vtxRes=0.005 efficiency=1."); // efficiency=1: all tracks found in trackers will get a vertex information
		fastSim->AddDetector("ScVtxNoMvd", "thtMin=0. thtMax=5.   ptmin=0.0 vtxRes=0.05  efficiency=1."); // efficiency=1: all tracks found in trackers will get a vertex information
	}
	else // MVD and GEM are disabled -> no good vertexing at all
	{
		fastSim->AddDetector("ScVtxNoMvd", "thtMin=0. thtMax=160. ptmin=0.1 vtxRes=0.1 efficiency=1."); // efficiency=1: all tracks found in trackers will get a vertex information
	}
	// -----------------------------------------------------------------------------------
	// EM Calorimeters w/ default parameters
	// (don't have to be set, just to list the available parameters
	// -----------------------------------------------------------------------------------

	fastSim->AddDetector("EmcFwCap", "thtMin=10.0 thtMax=22.0 Emin=0.01 dist=2.5");
	fastSim->AddDetector("EmcBwCap", "thtMin=142.0 thtMax=160.0 Emin=0.01 dist=0.7");

	if (SwEmcBar)
	{
		// EmcBarrel also allows to set phiMin and phiMax and can be added multiple times as EmcBarrel1, EmcBarrel2, etc.
		// Should be made constistent with EmcPidBarrel below
		fastSim->AddDetector("EmcBarrel","thtMin=22.0 thtMax=142.0 Emin=0.01 barrelRadius=0.5");
	}

	if (SwFwdSpec) // Fwd spectrometer enabled -> use Fwd EMC
	{
		fastSim->AddDetector("EmcFS",    "thtMin=0.05 thtMax=10.0 aPar=0.013 bPar=0.0283 Emin=0.01 dist=8.2");
	}

	// -----------------------------------------------------------------------------------
	// PID
	// -----------------------------------------------------------------------------------

	// PID detectors being always in: STT, MUO Barrel, EMC FwdCap, EMC BwdCap
	//Note: A dEdX parametrization from 2008
	fastSim->AddDetector("SttPid","thtMin=7.8 thtMax=159.5 ptmin=0.1 dEdxRes=0.15 efficiency=1.");
	fastSim->AddDetector("ScMdtPidBarrel", "thtMin=10.0 thtMax=130.0 pmin=0.5 efficiency=0.95 misId=0.01");
	fastSim->AddDetector("ScEmcPidFwCap",  "thtMin=10.0  thtMax=22.0  ptmin=0.0 pmin=0.0 efficiency=1.0");
	fastSim->AddDetector("ScEmcPidBwCap",  "thtMin=142.0 thtMax=160.0  ptmin=0.0 pmin=0.0 efficiency=1.0");

	if (SwMvdGem) // MVD and GEM are enabled -> MVD PID available
	{
		//Note: A Bethe-Bloch-Landau-Gauss Prametrization from 2008
		fastSim->AddDetector("MvdPid","thtMin=5.  thtMax=133.6 ptmin=0.1  dEdxResMulti=1. efficiency=1.");
	}

	if (SwEmcBar) // EMC Barrel enable -> EMC barrel PID available
	{
		fastSim->AddDetector("ScEmcPidBarrel", "thtMin=22.0  thtMax=142.0 ptmin=0.2 pmin=0.0 efficiency=1.0");
	}

	if (SwDrc) // Barrel DIRC enabled
	{
		fastSim->AddDetector("DrcBarrel","thtMin=22.0 thtMax=140.0 dthtc=0.01 nPhotMin=5 effNPhotons=0.075");
	}

	if (SwDsc) // Disc DIRC enabled
	{
		fastSim->AddDetector("DrcDisc","thtMin=5.0 thtMax=22.0 dthtc=0.01 nPhotMin=5 effNPhotons=0.075");
	}

	if (SwFwdSpec) // Fwd spectrometer enabled -> use RICH, FwdMUO and EMC FS
	{
		fastSim->AddDetector("ScEmcPidFS",     "thtMin=0.5   thtMax=10.0  ptmin=0.0 pmin=0.5 efficiency=1.0");
		fastSim->AddDetector("Rich","angleXMax=5.0 angleYMax=10.0 dthtc=0.01 nPhotMin=5 effNPhotons=0.075");
		fastSim->AddDetector("ScMdtPidForward","thtMin=0.0  thtMax=10.0  pmin=0.5 efficiency=0.95 misId=0.01");
	}


	fRun->AddTask(fastSim);
	
	
	// ***********************
	// *** SoftTriggerTask ***
	// ***********************
	
	if (runST)
	{	
		// this file contains the trigger line definitions
		TString triggercfg   = TString(gSystem->Getenv("VMCWORKDIR"))+"/softrig/triggerlines_fsim.cfg";
		
		PndSoftTriggerTask *stTask = new PndSoftTriggerTask(Mom, 0, 0, triggercfg);
		stTask->SetFastSimDefaults();
		stTask->SetQAEvent();
		fRun->AddTask(stTask);
	}
	
	// ***********************
	// *** SoftTriggerTask ***
	// ***********************

    
    
	// *****************************
	// *** PndSimpleCombinerTask ***
	// *****************************
		
	if (doreco)
	{
		PndSimpleCombinerTask *scTask = new PndSimpleCombinerTask(anadecay, anaparms+":algo=PidChargedProbability",Mom, run, runmode);
		scTask->SetPidAlgo("PidChargedProbability");
		fRun->AddTask(scTask);
	}
	
	// *****************************
	// *** PndSimpleCombinerTask ***
	// *****************************
	
	
	
	// *****************************
	// *** PndParticleQATask ***
	// *****************************
	
	if (partQA)
	{
		PndParticleQATask *partQaTask = new PndParticleQATask(kTRUE,chrg,neut,mc); // particle QA task for FastSim
		fRun->AddTask(partQaTask);
	}
	
	// *****************************
	// *** PndParticleQATask ***
	// *****************************
	
	
	//-------------------------  Initialize the RUN  -----------------
	fRun->Init();
	
	//-------------------------  Run the Simulation  -----------------
	fRun->Run(nEvents);
	
	//-------------------------  Write Filter Info to File -----------
	if (usePndEventFilter) primGen->WriteEvtFilterStatsToRootFile(); 
	
	//------------------------Print some info and exit----------------
	timer.Stop();
	Double_t rtime = timer.RealTime();
	Double_t ctime = timer.CpuTime();
	printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);
}

void getRange(TString par, double &min, double &max)
{
	par.ReplaceAll(" ","");
	par = par(par.Index("(")+1, par.Length()-par.Index("(")-2);
	
	TString smin=par, smax=par;
	
	if (par.Contains(",")) 
	{
		smin = par(0,par.Index(","));
		smax = par(par.Index(",")+1,1000);
	}
	
	min = smin.Atof();
	max = smax.Atof();
	
	//if (min>max) {double tmp=min;min=max;max=tmp;}
}

TString getInitialResonance(TString &fEvtGenFile)
{
  
  TString IniRes="";
  
  if (fEvtGenFile.Contains(":")) // is the initial resonance provide as <decfile>.dec:iniRes ? 
  {
    IniRes = fEvtGenFile(fEvtGenFile.Index(":")+1,1000);
    fEvtGenFile = fEvtGenFile(0,fEvtGenFile.Index(":"));
  }
  
  if (IniRes=="") // we need to search the decay file
  {
    std::ifstream fs(fEvtGenFile.Data());	
    char line[250];
  
    while (fs)
    {
      fs.getline(line,249);
      TString s(line);
      s.ReplaceAll("\r","");
      if (IniRes=="" && s.Contains("Decay "))
      {
        if (s.Contains("#")) s=s(0,s.Index("#"));
        s.ReplaceAll("Decay ","");
        s.ReplaceAll(" ","");
        IniRes = s;
      }	 
    } 
    fs.close();
  }
  
  return IniRes;
}
