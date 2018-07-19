# Installation

* I couldn't install the software on my own computer properly yet, I'll update this readme as soon as I'll be able to do that
* since I have a GSI account I could install it on GSI's computer farm

* first of all, I needed to ssh to the my account

```bash
$ ssh aolar@lx-pool.gsi.de
```

* this command selects a computer out of the lx-pool and drops me on that, I was prompted to give my password here
* from there, I needed to ssh to the khronos servers since they have Fairsoft and Fairroot set up properly

```bash
$ ssh aolar@kronos.hpc.gsi.de
```

* no password was requiered from here
* here I started to set up Fairsoft and Fairroot accordingly

```bash
$ export SIMPATH=/cvmfs/fairroot.gsi.de/fairsoft/oct17_root6/
$ export FAIRROOTPATH=/cvmfs/fairroot.gsi.de/fairroot/v-17.10a_fairsoft-oct17_root6/
```

* from here I was almost done, I used git to download pandaroot, I copied the line from the wiki page and executed the command in my home dir

```bash
$ git clone https://pandaatfair.githost.io/PandaRootGroup/PandaRoot.git ./pandaroot
```

* after this step I executed the following commands

```bash
$ cd pandaroot
$ mkdir buildPanda
$ cd buildPanda
$ cmake -DUSE_PATH_INFO=TRUE ../
$ make -j4
```

* after these steps I had the right tools in the pandaroot dir to execute a full simulation with reconstruction and digitalization

```bash
$ cd ~/pandaroot/macro/run
$ root -l sim_complete.C
$ root -l digi_complete.C
$ root -l reco_complete.C
$ root -l pid_complete.C
$ root -l ana_complete.C
```

* after each root command a `.q` was needed to quit *ROOT*