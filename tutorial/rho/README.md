# For detailed documentation of this tutorial directory, please consult [web resource](http://panda-wiki.gsi.de/cgi-bin/view/Computing/PandaRootRhoTutorial)

Run examples with:

## Simulation and reconstruction

```bash
root -l -b -q tut_sim.C         # full simulation of 100 events;    ==> signal_sim.root, signal_par.root
root -l -b -q tut_aod.C         # reconstruction of the 100 events; ==> signal_pid.root

./tut_runall.sh 100             # run both above

root -l -b -q tut_fastsim.C     # fast simulation of 1000 events;   ==> signal_fast.root
```

## Analysis

```bash
root -l -b -q tut_ana.C         # analysis of full sim events;      ==> signal_ana.root
root -l -b -q tut_ana_fast.C    # analysis of fast sim events       ==> signal_ana_fast.root

root -l -b -q tut_ana_pid.C     # analysis, PID by hand             ==> signal_ana_pid.root
root -l -b -q tut_ana_comb.C    # analysis, combinatorics by hand   ==> signal_ana_comb.root
root -l -b -q tut_ana_mcmatch.C # analysis, MC truth match          ==> signal_ana_mcmatch.root
root -l -b -q tut_ana_mclist.C  # print out MC truth list infos     ==> console text output 
root -l -b -q tut_ana_fit.C     # analysis, fitting                 ==> signal_ana_fit.root
root -l -b -q tut_ana_ntp.C     # analysis, n-Tuple output          ==> signal_ana_ntp.root
```

## Analysis in task

-> Add line 'add_subdirectory(tutorials/rho)' in $VMCWORKDIR/CMakeLists.txt somewhere after line 260ff
-> cd your/build/directory; make; cd -

```bash
root -l -b -q tut_ana_task.C    # analysis running in a task        ==> signal_ana_task.root
```

# Quick analysis tools

## Running on full sim output

```bash
root -l -b -q 'quickana.C("signal_pid.root", 6.232, "J/psi->mu+ mu-;pbarpSystem->J/psi pi+ pi-", 0, "fit4c:fitvtx:mwin(J/psi)=0.8")'    #  ==> signal_pid_ana.root

# Running on fast sim output
root -l -b -q 'quickana.C("signal_fast.root", 6.232, "J/psi->mu+ mu-;pbarpSystem->J/psi pi+ pi-", 0, "fit4c:fitvtx:mwin(J/psi)=0.8",1)' #  ==> signal_fast_ana.root

# Simulating and analysing fast simulation
root -l -b -q 'quickfsimana.C("signal", "pp_jpsi2pi_jpsi_mumu.dec", 6.232, "J/psi -> mu+ mu-; pbarpSystem -> J/psi pi+ pi-", 1000, "fit4c:fitvtx:mwin=0.8")' # ==> signal_0_ana.root  (1000 ev)
root -l -b -q 'quickfsimana.C("bkg", "DPM", 6.232, "J/psi -> mu+ mu-; pbarpSystem -> J/psi pi+ pi-", 10000, "fit4c:fitvtx:mwin=0.8")'                        # ==> bkg_o_ana.root    (10000 ev)

# Submit quickfsimana jobs to kronos
sbatch -a1-10 jobquickfa_kronos.sh MySig 'quickfsimana.C("PREFIX","pp_jpsi2pi_jpsi_mumu.dec",7.0,"J/psi->e+ e-;pbp0->J/psi pi+ pi-",10,"fit4c:fitvtx:mwin=0.8",0,RUN)'
sbatch -a1-10 jobquickfa_kronos.sh MyBkg 'quickfsimana.C("PREFIX","DPM",7.0,"J/psi->e+ e-;pbp0->J/psi pi+ pi-",10,"fit4c:fitvtx:mwin=0.8",0,RUN)'
```

