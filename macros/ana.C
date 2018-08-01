class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;
class RhoTuple;

void ana(int nevts = 10000)
{
    TString parAsciiFile = "all.par";
    TString prefix = "ppbar_psi2s2pi_jpsi2pi_2mu";
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

    TFile *out = TFile::Open(prefix + "_simple_output_ana.root", "RECREATE");

    RhoTuple *npbarp = new RhoTuple("npbarp", "pbarp System");

    PndAnalysis *theAnalysis = new PndAnalysis();
    if (nevts == 0)
        nevts = theAnalysis->GetEntries();

    RhoCandList pbarp, muplus, muminus, piplus, piminus, jpsi, psi2s, jpsiFit, psi2sFit;

    double m0_jpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass(); // Get nominal PDG mass of the J/psi
    RhoMassParticleSelector *jpsiMassSel = new RhoMassParticleSelector("jpsi", m0_jpsi, 1.0);

    TString pidSelection = "PidAlgoEmcBayes;PidAlgoDrc;PidAlgoDisc;PidAlgoStt;PidAlgoMdtHardCuts";

    // changed it accordingly to sim_complete beam momentum, inital lorentz vector of pbarpSystem
    // (p_x, p_y, p_z, E/c) style
    TLorentzVector ini1(0, 0, 9.231552, 9.279084);
    TLorentzVector ini2(0, 0, 0, 0.938);
    TLorentzVector ini(ini1 + ini2);

    std::cout << "**********************************************" << std::endl;
    std::cout << ini[2] << " " << ini[3] << std::endl;
    std::cout << "**********************************************" << std::endl;

    PndRhoTupleQA qa(theAnalysis, ini.P());

    int i = 0;
    while (theAnalysis->GetEvent() && i++ < nevts)
    {

        theAnalysis->FillList(muplus, "MuonAllPlus", pidSelection);
        theAnalysis->FillList(muminus, "MuonAllMinus", pidSelection);
        theAnalysis->FillList(piplus, "PionAllPlus", pidSelection);
        theAnalysis->FillList(piminus, "PionAllMinus", pidSelection);

        jpsi.Combine(muplus, muminus);
        jpsi.SetType(443); // checked - j/psi

        psi2s.Combine(jpsi, piplus, piminus);
        psi2s.SetType(100443); // checked - psi(2s)

        pbarp.Combine(psi2s, piplus, piminus);
        pbarp.SetType(88888); // checked - pbarpSystem

        for (int j = 0; j < pbarp.GetLength(); ++j)
        {
            qa.qaComp("pbarp_", pbarp[j], npbarp);

            RhoDecayTreeFitter dtf(pbarp[j], ini);
            dtf.setMassConstraint(pbarp[j]->Daughter(0)->Daughter(0));

            int check = dtf.Fit();

            double chiSqr = dtf.chiSquare();

            if (chiSqr > 0 && chiSqr < 100)
            {
                npbarp->Column("fitflag", (Float_t)check, -999.9f);
                npbarp->Column("fitchiq", (Float_t)chiSqr, -999.9f);
                qa.qaComp("pbarpfit_", pbarp[j]->GetFit(), npbarp);
            }

            npbarp->DumpData();
        }

        int threads = std::thread::hardware_concurrency();
        std::vector<std::future<void>> resultWaiter(threads);

        auto parallel_fill = [&](RhoCandListIterator beg) {
            RhoCandidate *cand = beg.Current();
            qa.qaComp("pbarp_", , cand, npbarp);
            RhoDecayTreeFitter dtf(cand, ini);
            dtf.setMassConstraint(cand->Daughter(0)->Daughter(0));
            int check = dtf.Fit();
            double chiSqr = dtf.chiSquare();
            if (chiSqr > 0 && chiSqr < 100)
            {
                npbarp->Column("fitflag", (Float_t)check, -999.9f);
                npbarp->Column("fitchiq", (Float_t)chiSqr, -999.9f);
                qa.qaComp("pbarpfit_", cand->GetFit(), npbarp);
            }
            while (cand = beg.Next())
            {
                qa.qaComp("pbarp_", , cand, npbarp);
                RhoDecayTreeFitter dtf(cand, ini);
                dtf.setMassConstraint(cand->Daughter(0)->Daughter(0));
                int check = dtf.Fit();
                double chiSqr = dtf.chiSquare();
                if (chiSqr > 0 && chiSqr < 100)
                {
                    npbarp->Column("fitflag", (Float_t)check, -999.9f);
                    npbarp->Column("fitchiq", (Float_t)chiSqr, -999.9f);
                    qa.qaComp("pbarpfit_", cand->GetFit(), npbarp);
                }
            }
            npbarp->DumpData();
            return;
        };

        std::vector<RhoCandListIterator> iterators(n);

        for (int k = 0; k < n; k++)
        {
            iterators[k] = new RhoCandListIterator(pbarp);
            resultWaiter[k] = std::async(std::launch::async, parallel_fill, iterators[k]);
        }

        std::for_each(resultWaiter.begin(),resultWaiter.end(),[&]( std::future<void>& f ){ std::cout << f.get(); });

    }

    out->cd();

    npbarp->GetInternalTree()->Write();

    out->Save();
}
class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;
class RhoTuple;

void ana_ntp_qa_complete_basic_aolar(int nevts = 10000)
{
    TString parAsciiFile = "all.par";
    TString prefix = "ppbar_psi2s2pi_jpsi2pi_2mu";
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

    TFile *out = TFile::Open(prefix + "_simple_output_ana.root", "RECREATE");

    RhoTuple *npbarp = new RhoTuple("npbarp", "pbarp System");

    PndAnalysis *theAnalysis = new PndAnalysis();
    if (nevts == 0)
        nevts = theAnalysis->GetEntries();

    RhoCandList pbarp, muplus, muminus, piplus, piminus, jpsi, psi2s, jpsiFit, psi2sFit;

    double m0_jpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass(); // Get nominal PDG mass of the J/psi
    RhoMassParticleSelector *jpsiMassSel = new RhoMassParticleSelector("jpsi", m0_jpsi, 1.0);

    TString pidSelection = "PidAlgoEmcBayes;PidAlgoDrc;PidAlgoDisc;PidAlgoStt;PidAlgoMdtHardCuts";

    // changed it accordingly to sim_complete beam momentum, inital lorentz vector of pbarpSystem
    // (p_x, p_y, p_z, E/c) style
    TLorentzVector ini1(0, 0, 9.231552, 9.279084);
    TLorentzVector ini2(0, 0, 0, 0.938);
    TLorentzVector ini(ini1 + ini2);

    std::cout << "**********************************************" << std::endl;
    std::cout << ini[2] << " " << ini[3] << std::endl;
    std::cout << "**********************************************" << std::endl;

    PndRhoTupleQA qa(theAnalysis, ini.P());

    int i = 0;
    while (theAnalysis->GetEvent() && i++ < nevts)
    {

        if (i % 10 == 0)
            std::cout << "EVENT: #" << i << std::endl;

        theAnalysis->FillList(muplus, "MuonAllPlus", pidSelection);
        theAnalysis->FillList(muminus, "MuonAllMinus", pidSelection);
        theAnalysis->FillList(piplus, "PionAllPlus", pidSelection);
        theAnalysis->FillList(piminus, "PionAllMinus", pidSelection);

        jpsi.Combine(muplus, muminus);
        jpsi.SetType(443); // checked - j/psi

        psi2s.Combine(jpsi, piplus, piminus);
        psi2s.SetType(100443); // checked - psi(2s)

        pbarp.Combine(psi2s, piplus, piminus);
        pbarp.SetType(88888); // checked - pbarpSystem

        for (int j = 0; j < pbarp.GetLength(); ++j)
        {
            qa.qaComp("pbarp_", pbarp[j], npbarp);

            RhoDecayTreeFitter dtf(pbarp[j], ini);
            dtf.setMassConstraint(pbarp[j]->Daughter(0)->Daughter(0));

            int check = dtf.Fit();

            double chiSqr = dtf.chiSquare();

            if (chiSqr > 0 && chiSqr < 100)
            {
                npbarp->Column("fitflag", (Float_t)check, -999.9f);
                npbarp->Column("fitchiq", (Float_t)chiSqr, -999.9f);
                qa.qaComp("pbarpfit_", pbarp[j]->GetFit(), npbarp);
            }

            npbarp->DumpData();
        }
        class RhoCandList;
        class RhoCandidate;
        class PndAnaPidSelector;
        class PndAnaPidCombiner;
        class PndAnalysis;
        class RhoTuple;

        void ana_ntp_qa_complete_basic_aolar(int nevts = 10000)
        {
            TString parAsciiFile = "all.par";
            TString prefix = "ppbar_psi2s2pi_jpsi2pi_2mu";
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

            TFile *out = TFile::Open(prefix + "_simple_output_ana.root", "RECREATE");

            RhoTuple *npbarp = new RhoTuple("npbarp", "pbarp System");

            PndAnalysis *theAnalysis = new PndAnalysis();
            if (nevts == 0)
                nevts = theAnalysis->GetEntries();

            RhoCandList pbarp, muplus, muminus, piplus, piminus, jpsi, psi2s, jpsiFit, psi2sFit;

            double m0_jpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass(); // Get nominal PDG mass of the J/psi
            RhoMassParticleSelector *jpsiMassSel = new RhoMassParticleSelector("jpsi", m0_jpsi, 1.0);

            TString pidSelection = "PidAlgoEmcBayes;PidAlgoDrc;PidAlgoDisc;PidAlgoStt;PidAlgoMdtHardCuts";

            // changed it accordingly to sim_complete beam momentum, inital lorentz vector of pbarpSystem
            // (p_x, p_y, p_z, E/c) style
            TLorentzVector ini1(0, 0, 9.231552, 9.279084);
            TLorentzVector ini2(0, 0, 0, 0.938);
            TLorentzVector ini(ini1 + ini2);

            std::cout << "**********************************************" << std::endl;
            std::cout << ini[2] << " " << ini[3] << std::endl;
            std::cout << "**********************************************" << std::endl;

            PndRhoTupleQA qa(theAnalysis, ini.P());

            int i = 0;
            while (theAnalysis->GetEvent() && i++ < nevts)
            {

                if (i % 10 == 0)
                    std::cout << "EVENT: #" << i << std::endl;

                theAnalysis->FillList(muplus, "MuonAllPlus", pidSelection);
                theAnalysis->FillList(muminus, "MuonAllMinus", pidSelection);
                theAnalysis->FillList(piplus, "PionAllPlus", pidSelection);
                theAnalysis->FillList(piminus, "PionAllMinus", pidSelection);

                jpsi.Combine(muplus, muminus);
                jpsi.SetType(443); // checked - j/psi

                psi2s.Combine(jpsi, piplus, piminus);
                psi2s.SetType(100443); // checked - psi(2s)

                pbarp.Combine(psi2s, piplus, piminus);
                pbarp.SetType(88888); // checked - pbarpSystem

                for (int j = 0; j < pbarp.GetLength(); ++j)
                {
                    qa.qaComp("pbarp_", pbarp[j], npbarp);

                    RhoDecayTreeFitter dtf(pbarp[j], ini);
                    dtf.setMassConstraint(pbarp[j]->Daughter(0)->Daughter(0));

                    int check = dtf.Fit();

                    double chiSqr = dtf.chiSquare();

                    if (chiSqr > 0 && chiSqr < 100)
                    {
                        npbarp->Column("fitflag", (Float_t)check, -999.9f);
                        npbarp->Column("fitchiq", (Float_t)chiSqr, -999.9f);
                        qa.qaComp("pbarpfit_", pbarp[j]->GetFit(), npbarp);
                    }

                    npbarp->DumpData();
                }
            }

            out->cd();

            npbarp->GetInternalTree()->Write();

            out->Save();
            class RhoCandList;
            class RhoCandidate;
            class PndAnaPidSelector;
            class PndAnaPidCombiner;
            class PndAnalysis;
            class RhoTuple;

            void ana_ntp_qa_complete_basic_aolar(int nevts = 10000)
            {
                TString parAsciiFile = "all.par";
                TString prefix = "ppbar_psi2s2pi_jpsi2pi_2mu";
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

                TFile *out = TFile::Open(prefix + "_simple_output_ana.root", "RECREATE");

                RhoTuple *npbarp = new RhoTuple("npbarp", "pbarp System");

                PndAnalysis *theAnalysis = new PndAnalysis();
                if (nevts == 0)
                    nevts = theAnalysis->GetEntries();

                RhoCandList pbarp, muplus, muminus, piplus, piminus, jpsi, psi2s, jpsiFit, psi2sFit;

                double m0_jpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass(); // Get nominal PDG mass of the J/psi
                RhoMassParticleSelector *jpsiMassSel = new RhoMassParticleSelector("jpsi", m0_jpsi, 1.0);

                TString pidSelection = "PidAlgoEmcBayes;PidAlgoDrc;PidAlgoDisc;PidAlgoStt;PidAlgoMdtHardCuts";

                // changed it accordingly to sim_complete beam momentum, inital lorentz vector of pbarpSystem
                // (p_x, p_y, p_z, E/c) style
                TLorentzVector ini1(0, 0, 9.231552, 9.279084);
                TLorentzVector ini2(0, 0, 0, 0.938);
                TLorentzVector ini(ini1 + ini2);

                std::cout << "**********************************************" << std::endl;
                std::cout << ini[2] << " " << ini[3] << std::endl;
                std::cout << "**********************************************" << std::endl;

                PndRhoTupleQA qa(theAnalysis, ini.P());

                int i = 0;
                while (theAnalysis->GetEvent() && i++ < nevts)
                {

                    if (i % 10 == 0)
                        std::cout << "EVENT: #" << i << std::endl;

                    theAnalysis->FillList(muplus, "MuonAllPlus", pidSelection);
                    theAnalysis->FillList(muminus, "MuonAllMinus", pidSelection);
                    theAnalysis->FillList(piplus, "PionAllPlus", pidSelection);
                    theAnalysis->FillList(piminus, "PionAllMinus", pidSelection);

                    jpsi.Combine(muplus, muminus);
                    jpsi.SetType(443); // checked - j/psi

                    psi2s.Combine(jpsi, piplus, piminus);
                    psi2s.SetType(100443); // checked - psi(2s)

                    pbarp.Combine(psi2s, piplus, piminus);
                    pbarp.SetType(88888); // checked - pbarpSystem

                    for (int j = 0; j < pbarp.GetLength(); ++j)
                    {
                        qa.qaComp("pbarp_", pbarp[j], npbarp);

                        RhoDecayTreeFitter dtf(pbarp[j], ini);
                        dtf.setMassConstraint(pbarp[j]->Daughter(0)->Daughter(0));

                        int check = dtf.Fit();

                        double chiSqr = dtf.chiSquare();

                        if (chiSqr > 0 && chiSqr < 100)
                        {
                            npbarp->Column("fitflag", (Float_t)check, -999.9f);
                            npbarp->Column("fitchiq", (Float_t)chiSqr, -999.9f);
                            qa.qaComp("pbarpfit_", pbarp[j]->GetFit(), npbarp);
                        }

                        npbarp->DumpData();
                    }
                }

                out->cd();

                npbarp->GetInternalTree()->Write();

                out->Save();
                class RhoCandList;
                class RhoCandidate;
                class PndAnaPidSelector;
                class PndAnaPidCombiner;
                class PndAnalysis;
                class RhoTuple;

                void ana_ntp_qa_complete_basic_aolar(int nevts = 10000)
                {
                    TString parAsciiFile = "all.par";
                    TString prefix = "ppbar_psi2s2pi_jpsi2pi_2mu";
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

                    TFile *out = TFile::Open(prefix + "_simple_output_ana.root", "RECREATE");

                    RhoTuple *npbarp = new RhoTuple("npbarp", "pbarp System");

                    PndAnalysis *theAnalysis = new PndAnalysis();
                    if (nevts == 0)
                        nevts = theAnalysis->GetEntries();

                    RhoCandList pbarp, muplus, muminus, piplus, piminus, jpsi, psi2s, jpsiFit, psi2sFit;

                    double m0_jpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass(); // Get nominal PDG mass of the J/psi
                    RhoMassParticleSelector *jpsiMassSel = new RhoMassParticleSelector("jpsi", m0_jpsi, 1.0);

                    TString pidSelection = "PidAlgoEmcBayes;PidAlgoDrc;PidAlgoDisc;PidAlgoStt;PidAlgoMdtHardCuts";

                    // changed it accordingly to sim_complete beam momentum, inital lorentz vector of pbarpSystem
                    // (p_x, p_y, p_z, E/c) style
                    TLorentzVector ini1(0, 0, 9.231552, 9.279084);
                    TLorentzVector ini2(0, 0, 0, 0.938);
                    TLorentzVector ini(ini1 + ini2);

                    std::cout << "**********************************************" << std::endl;
                    std::cout << ini[2] << " " << ini[3] << std::endl;
                    std::cout << "**********************************************" << std::endl;

                    PndRhoTupleQA qa(theAnalysis, ini.P());

                    int i = 0;
                    while (theAnalysis->GetEvent() && i++ < nevts)
                    {

                        if (i % 10 == 0)
                            std::cout << "EVENT: #" << i << std::endl;

                        theAnalysis->FillList(muplus, "MuonAllPlus", pidSelection);
                        theAnalysis->FillList(muminus, "MuonAllMinus", pidSelection);
                        theAnalysis->FillList(piplus, "PionAllPlus", pidSelection);
                        theAnalysis->FillList(piminus, "PionAllMinus", pidSelection);

                        jpsi.Combine(muplus, muminus);
                        jpsi.SetType(443); // checked - j/psi

                        psi2s.Combine(jpsi, piplus, piminus);
                        psi2s.SetType(100443); // checked - psi(2s)

                        pbarp.Combine(psi2s, piplus, piminus);
                        pbarp.SetType(88888); // checked - pbarpSystem

                        for (int j = 0; j < pbarp.GetLength(); ++j)
                        {
                            qa.qaComp("pbarp_", pbarp[j], npbarp);

                            RhoDecayTreeFitter dtf(pbarp[j], ini);
                            dtf.setMassConstraint(pbarp[j]->Daughter(0)->Daughter(0));

                            int check = dtf.Fit();

                            double chiSqr = dtf.chiSquare();

                            if (chiSqr > 0 && chiSqr < 100)
                            {
                                npbarp->Column("fitflag", (Float_t)check, -999.9f);
                                npbarp->Column("fitchiq", (Float_t)chiSqr, -999.9f);
                                qa.qaComp("pbarpfit_", pbarp[j]->GetFit(), npbarp);
                            }

                            npbarp->DumpData();
                        }
                    }

                    out->cd();

                    npbarp->GetInternalTree()->Write();

                    out->Save();
                    class RhoCandList;
                    class RhoCandidate;
                    class PndAnaPidSelector;
                    class PndAnaPidCombiner;
                    class PndAnalysis;
                    class RhoTuple;

                    void ana_ntp_qa_complete_basic_aolar(int nevts = 10000)
                    {
                        TString parAsciiFile = "all.par";
                        TString prefix = "ppbar_psi2s2pi_jpsi2pi_2mu";
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

                        TFile *out = TFile::Open(prefix + "_simple_output_ana.root", "RECREATE");

                        RhoTuple *npbarp = new RhoTuple("npbarp", "pbarp System");

                        PndAnalysis *theAnalysis = new PndAnalysis();
                        if (nevts == 0)
                            nevts = theAnalysis->GetEntries();

                        RhoCandList pbarp, muplus, muminus, piplus, piminus, jpsi, psi2s, jpsiFit, psi2sFit;

                        double m0_jpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass(); // Get nominal PDG mass of the J/psi
                        RhoMassParticleSelector *jpsiMassSel = new RhoMassParticleSelector("jpsi", m0_jpsi, 1.0);

                        TString pidSelection = "PidAlgoEmcBayes;PidAlgoDrc;PidAlgoDisc;PidAlgoStt;PidAlgoMdtHardCuts";

                        // changed it accordingly to sim_complete beam momentum, inital lorentz vector of pbarpSystem
                        // (p_x, p_y, p_z, E/c) style
                        TLorentzVector ini1(0, 0, 9.231552, 9.279084);
                        TLorentzVector ini2(0, 0, 0, 0.938);
                        TLorentzVector ini(ini1 + ini2);

                        std::cout << "**********************************************" << std::endl;
                        std::cout << ini[2] << " " << ini[3] << std::endl;
                        std::cout << "**********************************************" << std::endl;

                        PndRhoTupleQA qa(theAnalysis, ini.P());

                        int i = 0;
                        while (theAnalysis->GetEvent() && i++ < nevts)
                        {

                            if (i % 10 == 0)
                                std::cout << "EVENT: #" << i << std::endl;

                            theAnalysis->FillList(muplus, "MuonAllPlus", pidSelection);
                            theAnalysis->FillList(muminus, "MuonAllMinus", pidSelection);
                            theAnalysis->FillList(piplus, "PionAllPlus", pidSelection);
                            theAnalysis->FillList(piminus, "PionAllMinus", pidSelection);

                            jpsi.Combine(muplus, muminus);
                            jpsi.SetType(443); // checked - j/psi

                            psi2s.Combine(jpsi, piplus, piminus);
                            psi2s.SetType(100443); // checked - psi(2s)

                            pbarp.Combine(psi2s, piplus, piminus);
                            pbarp.SetType(88888); // checked - pbarpSystem

                            for (int j = 0; j < pbarp.GetLength(); ++j)
                            {
                                qa.qaComp("pbarp_", pbarp[j], npbarp);

                                RhoDecayTreeFitter dtf(pbarp[j], ini);
                                dtf.setMassConstraint(pbarp[j]->Daughter(0)->Daughter(0));

                                int check = dtf.Fit();

                                double chiSqr = dtf.chiSquare();

                                if (chiSqr > 0 && chiSqr < 100)
                                {
                                    npbarp->Column("fitflag", (Float_t)check, -999.9f);
                                    npbarp->Column("fitchiq", (Float_t)chiSqr, -999.9f);
                                    qa.qaComp("pbarpfit_", pbarp[j]->GetFit(), npbarp);
                                }

                                npbarp->DumpData();
                            }
                        }

                        out->cd();

                        npbarp->GetInternalTree()->Write();

                        out->Save();
                    