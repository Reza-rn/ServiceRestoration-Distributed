#*******************************************************************************************************************
__author__ = "Reza Roofegari nejad"                                                                          #******
__email__ = "r.roofegari@knights.ucf.edu"                                                                    #******
__website__ = "http://www.ece.ucf.edu/~rezarn/"                                                              #******
__copyright__ = "Copyright 2019, Scalable/Secure Cooperative Algorithms and Framework for Extremely-high " \
                "Penetration Solar Integration (SolarExPert), ENERGISE Project"                              #******
__ProjectWebsite__ = "https://www.cs.ucf.edu/~qu/MA-OpenDSS.php"                                             #******
#*******************************************************************************************************************
# This code is to manage clusters and find the neighbor bus of each cluster within a defined network
# ***************************************************************************************************************
# Import network and all data
# ***************************************************************************************************************
from OpenDSSCkt import *
from scipy.io import loadmat

#print(AllBuses)

# ***************************************************************************************************************
# Clustering Functions
# ***************************************************************************************************************
# ******************************************************************************************
#  *** This function is used to read clusters from MATLAB
# ******************************************************************************************
def FindClusters(loadmat):
    InitialClustering = loadmat('C:/Users/re208655/OneDrive - Knights - University of Central Florida/PhD studies/Reserach/Reconfiguration/Distributed/Simulation/Clusters/8500Node/Version1_Polish_MIP_PVSystem/TestNetwork/8500Node(OpenDSS)/Clusters.mat')

    NumClusters = len(InitialClustering['Clusters'])

    Clusters = {}
    for i in range(NumClusters):
        CL_Name = 'Cluster' + str(i + 1)
        ClusterLength = len(InitialClustering['Clusters'][i])
        # Find empty elements and reduce from Cluster length
        m = 0
        for k in range(ClusterLength):
            #if InitialClustering['Clusters'][i][k] is []:
            if isinstance(InitialClustering['Clusters'][i][k][0], str) == False:
                m = m + 1
        ClusterLength = ClusterLength - m

        for j in range(ClusterLength):
            if CL_Name not in Clusters:
                Clusters[CL_Name] = [InitialClustering['Clusters'][i][j][0]]
            else:
                Clusters[CL_Name].append(InitialClustering['Clusters'][i][j][0])


    return [Clusters, NumClusters]



# ******************************************************************************************
#  *** This function is used to find cluster of each bus and build BusClsuters dictionary
# ******************************************************************************************
def FindBusClsuters(Clusters):
    Bus_Clusters = {}

    for CL, Buses in Clusters.items():
        for bus in Buses:
            Bus_Clusters[bus] = CL

    return Bus_Clusters

# ******************************************************************************************
#  *** This function is used to find neighbor buses of each cluster
# ******************************************************************************************
def FindClustersNeighborBuses(Clusters, ParentBus, ChildrenBuses):
    ClustersNeighbor = {}
    ClustersNeighborRelations = {}
    BusesinClusters = {}
    ClusterbusVariables = {}

    for cluster, Buses in Clusters.items():
        ClustersNeighbor[cluster] = {}
        ClustersNeighborRelations[cluster] = {}
        ClusterbusVariables[cluster] = []
        # Find Parent Bus
        for bus in Buses:
            BusesinClusters[bus] = cluster
            ClusterbusVariables[cluster].append(bus)
            if bus in ParentBus:
                for ParentTemp in ParentBus[bus]:
                    if ParentTemp not in Buses:
                        if bus not in ClustersNeighbor[cluster]:
                            ClustersNeighbor[cluster][bus] = [ParentTemp]
                        else:
                            ClustersNeighbor[cluster][bus].append(ParentTemp)

                        if ParentTemp not in ClusterbusVariables[cluster]:
                            ClusterbusVariables[cluster].append(ParentTemp) # Add bus to the ClusterbusVariables
                        # Also add bus to the ClustersNeighborRelations
                        if 'ParentNeighborBuses' not in ClustersNeighborRelations[cluster]:
                            ClustersNeighborRelations[cluster]['ParentNeighborBuses'] = [ParentTemp]
                        else:
                            ClustersNeighborRelations[cluster]['ParentNeighborBuses'].append(ParentTemp)

        # Find Children Buses
            if bus in ChildrenBuses:
                for child in ChildrenBuses[bus]:
                    if child not in Buses:
                        if bus not in ClustersNeighbor[cluster]:
                            ClustersNeighbor[cluster][bus] = [child]
                        else:
                            ClustersNeighbor[cluster][bus].append(child)
                        # Also add bus to the ClustersNeighborRelations
                        if 'ChildrenNeighborBuses' not in ClustersNeighborRelations[cluster]:
                            ClustersNeighborRelations[cluster]['ChildrenNeighborBuses'] = [child]
                        else:
                            ClustersNeighborRelations[cluster]['ChildrenNeighborBuses'].append(child)

    return [ClustersNeighbor, ClustersNeighborRelations, BusesinClusters, ClusterbusVariables]

# ******************************************************************************************
#  *** This function is used to find loads in each cluster
# ******************************************************************************************
def FindClusterLoads(Clusters, LoadToBus):
    ClusterLoads = {}
    ClusterLoadsNonDispatchable = {}
    for CL, Buses in Clusters.items():
        ClusterLoads[CL] = []
        ClusterLoadsNonDispatchable[CL] = []
        for bus in Buses:
            if bus in LoadToBus:
                # All Loads
                for load in LoadToBus[bus]:
                    ClusterLoads[CL].append(load)
                # Non-Dispatchable Loads
                    if LoadType[load] == 0:
                        ClusterLoadsNonDispatchable[CL].append(load)

    return [ClusterLoads, ClusterLoadsNonDispatchable]

# ******************************************************************************************
#  *** This function is used to find PD Elements within each cluster
# ******************************************************************************************
def FindClusterPDElements(Clusters, BusNeighborPDElements, ClustersNeighborBuses, PDElementsConnections,
                          ClustersNeighborRelations, RemovedPDElements, BusesInClusters):
    ClusterPD = {}
    ClusterBoundaryPD={}
    clusterLineVariables = {}
    ClusterPDReduced = {}
    for CL, Buses in Clusters.items():
        ClusterPD[CL] = []
        ClusterBoundaryPD[CL] = {}
        clusterLineVariables[CL] = []
        ClusterPDReduced[CL] = []
        ClustersNeighborRelations[CL]['ParentPDElements'] = []
        ClustersNeighborRelations[CL]['ChildrenPDElements'] = []
        for bus in Buses:
            if bus in BusNeighborPDElements:
                for PDElement in BusNeighborPDElements[bus]:
                    if PDElement not in ClusterPD[CL]:
                        ClusterPD[CL].append(PDElement)
                        clusterLineVariables[CL].append(PDElement)
                    if PDElement not in ClusterPDReduced[CL]:
                        if PDElement not in RemovedPDElements:
                            ClusterPDReduced[CL].append(PDElement)
                # Find Boundary PD Element
                    if bus in ClustersNeighborBuses[CL]:
                    # Find the next terminal of PD Element------
                        if PDElement in PDElementsConnections:
                            for i in PDElementsConnections[PDElement]:
                                if i != bus:
                                    bus2 = i
                    # --------------------------
                        if bus2 not in Clusters[CL]:
                            ClusterBoundaryPD[CL][PDElement] = BusesInClusters[bus2]

                            # find parent or children PD Elements and add it to ClustersNeighborRelations
                            # Parent PD Elements-------
                            if 'ParentNeighborBuses' in ClustersNeighborRelations[CL]:
                                if bus2 in ClustersNeighborRelations[CL]['ParentNeighborBuses']:
                                    ClustersNeighborRelations[CL]['ParentPDElements'].append(PDElement)
                            # Children PD Elements----
                            if 'ChildrenNeighborBuses' in ClustersNeighborRelations[CL]:
                                if bus2 in ClustersNeighborRelations[CL]['ChildrenNeighborBuses']:
                                    ClustersNeighborRelations[CL]['ChildrenPDElements'].append(PDElement)
                                    clusterLineVariables[CL].remove(PDElement)

    return [ClusterPD, ClusterBoundaryPD, ClustersNeighborRelations, clusterLineVariables, ClusterPDReduced]

# ******************************************************************************************
#  *** This function is used to find PD Elements within each cluster
# ******************************************************************************************
def FindClusterSourceBus(Clusters, SourceBus):
    ClustersSourcebus = {}
    for CL, Buses in Clusters.items():
        ClustersSourcebus[CL] = []
        for bus in Buses:
            if bus in SourceBus:
                ClustersSourcebus[CL].append(bus)

    return ClustersSourcebus

# ******************************************************************************************
#  *** This function is used to find PVs within each cluster
# ******************************************************************************************
def FindClusterPVs(Clusters, PVBus):
    ClustersPV = {}
    for CL, Buses in Clusters.items():
        ClustersPV[CL] = []
        for bus in Buses:
            if bus in PVBus:
                if PVBus[bus] != []:
                    for PV in PVBus[bus]:
                        ClustersPV[CL].append(PV)

    return ClustersPV

# ******************************************************************************************
#  *** This function is used to find switchable lines within each cluster
# ******************************************************************************************
def FindClusterSwitchableLines(ClustersPDElements, Switches_Lines, ClustersBoundaryPDElements):
    ClusterSwitchLines = {}
    ClustersBoundarySwitchLines = {}
    for CL, PDElements in ClustersPDElements.items():
        ClusterSwitchLines[CL] = []
        ClustersBoundarySwitchLines[CL] = []
        for line in PDElements:
            if line in Switches_Lines:
                ClusterSwitchLines[CL].append(line)
                if line in ClustersBoundaryPDElements[CL]:
                    ClustersBoundarySwitchLines[CL].append(line)

    return [ClusterSwitchLines, ClustersBoundarySwitchLines]

# ******************************************************************************************
#  *** This function is used to Modify Cluster Neighbor Relations Data structure
# ******************************************************************************************
def ModifyClusterNeighborRelations(ClustersNeighborRelations, BusesInClusters, ClustersBoundaryPDElements):

    for cluster, DataStr in ClustersNeighborRelations.items():

        # Find parent clusters
        if 'ParentNeighborBuses' in ClustersNeighborRelations[cluster]:
            for bus in ClustersNeighborRelations[cluster]['ParentNeighborBuses']:
                if 'ParentClusters' not in ClustersNeighborRelations[cluster]:
                    ClustersNeighborRelations[cluster]['ParentClusters'] = [BusesInClusters[bus]]
                else:
                    if cluster not in ClustersNeighborRelations[cluster]['ParentClusters']:
                        ClustersNeighborRelations[cluster]['ParentClusters'].append(BusesInClusters[bus])

        # Find children Clusters
        if 'ChildrenNeighborBuses' in ClustersNeighborRelations[cluster]:
            for bus in ClustersNeighborRelations[cluster]['ChildrenNeighborBuses']:
                if 'ChildrenClusters' not in ClustersNeighborRelations[cluster]:
                    ClustersNeighborRelations[cluster]['ChildrenClusters'] = [BusesInClusters[bus]]
                else:
                    if cluster not in ClustersNeighborRelations[cluster]['ChildrenClusters']:
                        ClustersNeighborRelations[cluster]['ChildrenClusters'].append(BusesInClusters[bus])

    return ClustersNeighborRelations

# ******************************************************************************************
#  *** This function is used to find capacito banks in each cluster
# ******************************************************************************************
def FindClusterCapacitorBanks(Clusters, CapacitorBus):
    ClusterCB = {}
    for CL, Buses in Clusters.items():
        ClusterCB[CL] = []
        for bus in Buses:
            if bus in CapacitorBus:
                if CapacitorBus[bus] != []:
                    for CB in CapacitorBus[bus]:
                        ClusterCB[CL].append(CB)

    return ClusterCB
# ******************************************************************************************
#  *** This function is used to draw clusters
# ******************************************************************************************
def PlotClusters(Clusters, PDElementsConnections, XY_Coordinates, plt, BusNameEnable):

    Colors = ['g', 'r', 'c', 'm', 'y', 'k','gray','orange','purple','brown','blue']
    NumberofColor = len(Colors)
    plt.figure(1)
    # Draw Lines---------------
    for PDElements, Buses in PDElementsConnections.items():
        if (XY_Coordinates[Buses[0]][0] != 0) and (XY_Coordinates[Buses[0]][1] != 0):
            if (XY_Coordinates[Buses[1]][0] != 0) and (XY_Coordinates[Buses[1]][1] != 0):
                TempX = [XY_Coordinates[Buses[0]][0], XY_Coordinates[Buses[1]][0]]
                TempY = [XY_Coordinates[Buses[0]][1], XY_Coordinates[Buses[1]][1]]
                NumPhases = sum(PDElementPhaseMatrix[PDElements])
                if NumPhases == 3:
                    plt.plot(TempX, TempY, 'b-', linewidth=0.8)
                elif NumPhases == 2:
                    plt.plot(TempX, TempY, 'b--', linewidth=0.8)
                else:
                    plt.plot(TempX, TempY, 'b:', linewidth=0.8)

    # Draw nodes with different colors to separate clusters
    ColorIndex = 0
    ClusterNumber = 1
    for Cluster, Buses in Clusters.items():
        ClusterColor = Colors[ColorIndex]
        for Bus in Buses:
            if (XY_Coordinates[Bus][0] != 0) and (XY_Coordinates[Bus][1] != 0):
                plt.scatter(XY_Coordinates[Bus][0], XY_Coordinates[Bus][1], c=ClusterColor, s=10)
                if BusNameEnable:
                    plt.text(XY_Coordinates[Bus][0] - 10, XY_Coordinates[Bus][1] - 100, '{}'.format(Bus),fontsize=5)

        # Cluster Number
        plt.text(XY_Coordinates[Buses[0]][0] - 1200, XY_Coordinates[Buses[0]][1] - 500, 'CL{}'.format(ClusterNumber), fontsize=8, color=ClusterColor)
        ClusterNumber += 1

        ColorIndex += 1
        if ColorIndex >= NumberofColor:
            ColorIndex = 0

    # ---------------------
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()



# ***************************************************************************************************************
# End of Functions
# ***************************************************************************************************************

# ***************************************************************************************************************
# Read Clustering Results from Matlab
# ***************************************************************************************************************
Clusters = FindClusters(loadmat)
NumClusters = Clusters[1]
Clusters = Clusters[0]

# ***************************************************************************************************************
# Find neighboring buses in each cluster
# ***************************************************************************************************************
ClustersNeighborBuses = FindClustersNeighborBuses(Clusters, ParentBusAll, ChildrenBuses)
ClusterBusVariables = ClustersNeighborBuses[3]
BusesInClusters = ClustersNeighborBuses[2]
ClustersNeighborRelations = ClustersNeighborBuses[1]
ClustersNeighborBuses = ClustersNeighborBuses[0]
# ***************************************************************************************************************
# Find all loads in each cluster
# ***************************************************************************************************************
ClustersLoads = FindClusterLoads(Clusters, LoadToBus)
ClustersLoadsNonDispatchable = ClustersLoads[1]
ClustersLoads = ClustersLoads[0]


# ***************************************************************************************************************
# Find all PD Elements within each cluster
# ***************************************************************************************************************
ClustersPDElements = FindClusterPDElements(Clusters, BusNeighborPDElements, ClustersNeighborBuses, PDElementsConnections,
                                           ClustersNeighborRelations, RemovedPDElements, BusesInClusters)
ClustersPDElementsReduced = ClustersPDElements[4]
ClusterCurrentVariables = ClustersPDElements[3]
ClustersNeighborRelations = ClustersPDElements[2]  # part of modification to find parent and children PD Elements is done here
ClustersBoundaryPDElements = ClustersPDElements[1]
ClustersPDElements = ClustersPDElements[0]
# ***************************************************************************************************************
# Find Source buses within each cluster
# ***************************************************************************************************************
ClustersSourceBus = FindClusterSourceBus(Clusters, SourceBus)

# ***************************************************************************************************************
# Find PVs within each cluster
# ***************************************************************************************************************
ClustersPVs = FindClusterPVs(Clusters, PVBus)

# ***************************************************************************************************************
# Find Switchable lines within each cluster
# ***************************************************************************************************************
ClustersSwitchableLines = FindClusterSwitchableLines(ClustersPDElements, Switches_Lines, ClustersBoundaryPDElements)
ClustersBoundarySwitchableLines = ClustersSwitchableLines[1]
ClustersSwitchableLines = ClustersSwitchableLines[0]
# ***************************************************************************************************************
# Modify ClusterNeighborRelations
# ***************************************************************************************************************
ClustersNeighborRelations = ModifyClusterNeighborRelations(ClustersNeighborRelations, BusesInClusters, ClustersBoundaryPDElements)

# ***************************************************************************************************************
# Find all PD Elements within each cluster
# ***************************************************************************************************************

# ***************************************************************************************************************
# Find capacitor banks of each cluster
# ***************************************************************************************************************
ClustersCBs = FindClusterCapacitorBanks(Clusters, CapacitorBus)

# ***************************************************************************************************************
# Plot Clusters
# ***************************************************************************************************************
PlotClusters(Clusters, PDElementsConnections, XY_Coordinates, plt, 0)
# ----------------------------------------------------------------------------------------------------------------
# Results
#print(Clusters)

