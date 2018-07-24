# Time based simulation

* Events described in MC files, these include a *background*, and any number of *signals*

* Time structure is calculated in **digitization**

* PANDA detects a continuous beam with Poisson statistics, there is a very short time between each event so the detectors are still busy when a new hit should be encountered

## Time variables for an event

* *event time* is the actual mean time of an event, it can be accessed by $FairRootManager::Instance()->GetEventTime()$ function

  * important to note that there seems to be an *event manager* that can access all important data of an event
