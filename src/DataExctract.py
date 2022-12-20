
# Define packages
import pickle
import math, sys

# Load data-------------------------
FolderName = sys.argv[0].split('/')[:-1]
FolderName = "/".join(FolderName)
FolderName2="C:/Users/re208655/OneDrive - Knights - University of Central Florida/PhD studies/Reserach/Reconfiguration/Distributed/Simulation/Clusters/Version7/TestNetwork/123Bus(OpenDSS)/WorkSpace1.pickle";
FileFolderName = FolderName + "/Figures/Scenario1/WorkSpace1.pickle";
pickle_in = open(FolderName2,"rb")
DDSR = pickle.load(pickle_in)

#Restoration time steps---------------------
RestorationTimes = list(DDSR['X'].keys())

#Clusters------------------------------------
Clusters = list(list(DDSR['X'][RestorationTimes[0]].keys()))

# Extract total loads-------------------
TotalLoads = []
for t in RestorationTimes:
    TotalLoads.append(0)

for t in RestorationTimes:
    for cluster in Clusters:
        TotalLoads[RestorationTimes.index(t)] += DDSR['X'][t][cluster]['P_LoadTotal']





