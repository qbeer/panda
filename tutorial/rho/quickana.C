// ********************************************
// Quick analysis using PndSimpleCombinerTask/PndSimpleCombiner
// ********************************************

// The parameters are
// -------------------
// USAGE
// quickana.C( <input>, <mom>, <decay>, [nevt], [parms], [fastsim], [runST], [runnum], [mode] )
//    <input>   : input file name with PndPidCandidates
//    <mom>     : pbar momentum; negative values are interpreted as -E_cm
//    <decay>   : the decay pattern to be reconstructed, e.g. 'phi -> K+ K-; D_s+ -> phi pi-'
//    [nevt]    : number of events; default: 0 = all
//    [parms]   : parameters for the analysis, e.g. 'mwin=0.4:mwin(phi)=0.1:emin=0.1:pmin=0.1:qamc'
//    [fastsim] : set true, if running fast sim (sets the PID algos properly); default: false'
//    [runST]   : if 'true' runs Software Trigger (default: false)
//    [runnum]  : integer run number (default: 0)
//    [mode]    : arbitrary mode number (default: 0)
// -------------------

void quickana(TString Fname="", double Mom=0, TString anadecay="", int nevts=0, TString anaparms="", bool fastsim=false, bool runST=false, int run=0, int runmode=0)
{
	if (Fname=="" || Mom==0) 
	{
		cout << "USAGE:\n";
		cout << "quickana.C( <input>, <mom>, <decay>, [nevt], [parms], [fastsim], [runST], [runnum], [mode] )\n\n";
		cout << "   <input>   : input file name with PndPidCandidates\n";
		cout << "   <mom>     : pbar momentum; negative values are interpreted as -E_cm\n";
		cout << "   <decay>   : the decay pattern to be reconstructed, e.g. 'phi -> K+ K-; D_s+ -> phi pi-'\n";
		cout << "   [nevt]    : number of events; default: 0 = all\n";
		cout << "   [parms]   : parameters for the analysis, e.g. 'mwin=0.4:mwin(phi)=0.1:emin=0.1:pmin=0.1:qamc'; 'qapart' runs particle QA task\n";
		cout << "   [fastsim] : set true, if running fast sim (sets the PID algos properly); default: false'\n";
		cout << "   [runST]   : if 'true' runs Software Trigger (default: false)\n";
		cout << "   [runnum]  : integer run number (default: 0)\n";
		cout << "   [mode]    : arbitrary mode number (default: 0)\n\n";
		return;
	}
	
	// do some reconstruction ?
	bool doreco  = (anadecay != "" || anaparms.Contains("nevt"));

	// do particle QA?
	bool partQA  = (anaparms.Contains("qapart"));
	bool mc      = !(anaparms.Contains("!mc")) && !(anaparms.Contains("qamc"));
	bool neut    = !(anaparms.Contains("!neut"));
	bool chrg    = !(anaparms.Contains("!chrg"));
	
	// if Mom<0, interprete as -E_cm
	double mp = 0.938272;
	
	// if mom<0, it's -E_cm -> compute mom
	if (Mom<0)
	{
		double X = (Mom*Mom-2*mp*mp)/(2*mp);
		Mom = sqrt(X*X-mp*mp);
	}
	
	// PID algorithm for the PndSimpleCombinerTask (for Eventshape variables)
	TString pidalgo = "PidAlgoEmcBayes;PidAlgoDrc;PidAlgoDisc;PidAlgoStt;PidAlgoMdtHardCuts;PidAlgoRich;PidAlgoSciT";
	if (fastsim) pidalgo = "PidChargedProbability";
	
	// allow shortcuts
	anadecay.ReplaceAll("pbp","pbarpSystem");
	anadecay.ReplaceAll("pbp0","pbarpSystem0");
	
	// Prevent generator from throwing a lot of warnings
	//TLorentzVector fIni(0,0,Mom,0.938272+sqrt(Mom*Mom+0.938272*0.938272));
	TDatabasePDG::Instance()->AddParticle("pbarpSystem","pbarpSystem",3,kFALSE,0.1,0, "",88888);
	TDatabasePDG::Instance()->AddParticle("pbarpSystem0","pbarpSystem0",3,kFALSE,0.1,0, "",88880);
	
	// *** set this to your output path
	TString OutFile = Fname;//(Fname.Last('/')+1,Fname.Length()); // cut away input path
	OutFile.ReplaceAll(".root","_ana.root");
	
	// *** the output file for FairRunAna
	TString InFile  = Fname;
	if (!InFile.EndsWith(".root")) InFile+="_fast.root";
					
	// *** initialization
	FairLogger::GetLogger()->SetLogToFile(kFALSE);

	FairRunAna* fRun = new FairRunAna();
	fRun->SetGenerateRunInfo(kFALSE);
	fRun->SetSource(new FairFileSource(InFile));
	fRun->SetOutputFile(OutFile);

	// *** take constant field; needed for PocaVtx
	RhoCalculationTools::ForceConstantBz(20.0);
	
	// ***********************
	// *** SoftTriggerTask ***
	// ***********************
	
	if (runST)
	{	
		// this file contains the trigger line definitions
		TString      triggercfg = TString(gSystem->Getenv("VMCWORKDIR"))+"/softrig/triggerlines.cfg";       // fullsim trigger definitions 		
		if (fastsim) triggercfg = TString(gSystem->Getenv("VMCWORKDIR"))+"/softrig/triggerlines_fsim.cfg";  // fastsim trigger definitions	
		
		PndSoftTriggerTask *stTask = new PndSoftTriggerTask(Mom, 0, run, triggercfg);
		
		if (fastsim) stTask->SetFastSimDefaults();
		else         stTask->SetFullSimDefaults();
				
		fRun->AddTask(stTask);
	}
	
	// --------------------------------
	// *** Analysis Task ***
	// --------------------------------

	// *****************************
	// *** PndSimpleCombinerTask ***
	// *****************************
	
	if (doreco)
	{
		if (fastsim) anaparms+=":algo="+pidalgo;
		PndSimpleCombinerTask *scTask = new PndSimpleCombinerTask(anadecay, anaparms, Mom, run, runmode);
		scTask->SetPidAlgo(pidalgo);
		fRun->AddTask(scTask);
	}
	
	// *****************************
	// *** PndParticleQATask ***
	// *****************************
	
	if (partQA)
	{
		PndParticleQATask *partQaTask = new PndParticleQATask(fastsim,chrg,neut,mc); // particle QA task
		fRun->AddTask(partQaTask);
	}
		

	// *** and run analysis
	fRun->Init(); 
	fRun->Run(0,nevts);	
}
