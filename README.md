# HN_simulation
An analysis framework to study the sensitivity of the LHCb experiment to heavy neutrino (HN) candidates.
A set of signal samples with user-defined HN masses will be generated, as well as the dominant Standard Model backgrounds. Various kinematic plots are produced for signal-background discrimination. Using the theoretical calculations of the HN decay time and branching ratio as a function of HN mass and coupling, the expected number of HN candidates per year as a function of HN mass and coupling is produced, providing the expected reach of LHCb in Run3+4.

## Getting started
By default, the code decays the HN to a muon and a tau lepton, and decays the tau to a muon and neutrinos, but this can be modified. Three scenarios of muonic decay of tau, decay of the tau to a charged pion and a neutrino, and gen-level tau information (equivalent to not decaying) are implemented. To run the muon scenario, edit first few lines of `signal_muon/signal_muon_production.sh` to provide your desired list of HN mass points (by default, masses=(2.2 3.0 4.0 5.0 4.5)). Then run `source muon.sh [nEvents] [generate]` to run the full analysis. For the next trials, the last paramater (`generate`) can be omitted to skip the regeneration of samples and just redo the analysis. The core analysis code can be found in `signal_muon/mastercode_signal_muon.py`.

