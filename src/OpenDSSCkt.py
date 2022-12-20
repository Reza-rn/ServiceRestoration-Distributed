#*******************************************************************************************************************
__author__ = "Reza Roofegari nejad"                                                                          #******
__email__ = "r.roofegari@knights.ucf.edu"                                                                    #******
__website__ = "http://www.ece.ucf.edu/~rezarn/"                                                              #******
__copyright__ = "Copyright 2019, Scalable/Secure Cooperative Algorithms and Framework for Extremely-high " \
                "Penetration Solar Integration (SolarExPert), ENERGISE Project"                              #******
__ProjectWebsite__ = "https://www.cs.ucf.edu/~qu/MA-OpenDSS.php"                                             #******
#*******************************************************************************************************************
from OpenDSSInit import *
from Functions import *
import math, sys
import numpy as np
import matplotlib.pyplot as plt
#plt.rcParams['figure.dpi'] = 400#200
import time
# ****************************************************
# * Start of time
# ****************************************************
StartTime = time.time ()

# ****************************************************
# * DSS Object
# ****************************************************
# create DSS object of the circuit
Ckt = DSS()


# Load IEEE 123-node circuit
FolderName = sys.argv[0].split('/')[:-1]
FolderName = "/".join(FolderName)
Ckt.dssText.Command = r"Compile '"+FolderName+ "/TestNetwork/8500Node(OpenDSS)/Master.dss"
#Ckt.dssText.Command = r"Compile 'C:\Users\re-208655\OneDrive - Knights - University of Central Florida\PhD studies\Reserach\Reconfiguration\Distributed\Simulation\Clusters\Version1\TestNetwork\123Bus(OpenDSS)\IEEE123Master.dss"

# EnergyMeter
Ckt.dssText.command = "New Energymeter.m1 Line.ln5815900-1 1"

#==================================================================================
# Set Maximum iterations
Ckt.dssText.command = "Set Maxiterations=30"

# Solve the Circuit
Ckt.Solve()
if Ckt.dssSolution.Converged:
    #print("The Circuit Solved Successfully")
    print("The OpenDSS Connection Established Successfully")
else:
    print('Unsuccessful Connection')

# ****************************************************
# * End defining circuit
# ****************************************************
# ************************************************
# Get all the list of all buses, lines and loads
AllLines = Ckt.dssCircuit.Lines.AllNames
AllBuses = Ckt.dssCircuit.AllBusNames
AllLoads = Ckt.dssCircuit.Loads.AllNames
# print(AllBuses)
CktIncMatrix = Ckt.IncMatrix()
AllPDElements = Ckt.FindAllPDElements()


# ****************************************************
# * Find parent Buses of each Bus
# ****************************************************
# Calculate PD Elements connections
PDElementsConnections = PDElementsConnectionsFind(CktIncMatrix, AllPDElements,
                                                  Ckt.dssCircuit.Solution.IncMatrixcols)

# Find Parent set
ParentBus = FindParentBus(PDElementsConnections)
ParentBusAll = FindParentBusAll(PDElementsConnections)

#Reduced PD Elements
RemovedPDElements = FindRemovedPDElemnts(PDElementsConnections)
AllPDElementsReduced = RemovedPDElements[0]
PDElementsConnectionsReduced = RemovedPDElements[1]
ReplacedRemovedPDElements = RemovedPDElements[2]
RemovedPDElements = RemovedPDElements[3]

# ****************************************************
# * Find children Buses of each Bus
# ****************************************************
ChildrenBuses = FindChildrenBuses(PDElementsConnections)

# ****************************************************
# * Find neighbor PD Elements of each bus
# ****************************************************
BusNeighborPDElements = FindBusNeighborPDElements(AllBuses, PDElementsConnections)

# ****************************************************
# * Find neighbor Buses each Bus
# ****************************************************
BusNeighborBuses = FindBusNeighborBuses(AllBuses, BusNeighborPDElements, PDElementsConnections)

# ****************************************************
# * Find children lines of each line
# ****************************************************
ChildrenLines = FindChildrenLines(PDElementsConnections)

# ****************************************************
# * Find PD Element Phase Matrix
# ****************************************************
PDElementPhaseMatrix = FindPDElementPhaseMatrix(Ckt, AllPDElements)

# ****************************************************
# * Find Buses Phase Matrix
# ****************************************************
BusesPhaseMatrix = FindBusesPhaseMatrix(Ckt, AllBuses)
# print(BusesPhaseMatrix['300'])


# *************************************************************
# * Find S_Base and V_Base and Z_Base and circuit parameters
# *************************************************************
#SourceBus = Ckt.FindSourceBus()
SourceBus = ['sourcebus']
SubPDElement = FindFirstPDElement(PDElementsConnections, SourceBus)  # find connected PD element to the source bus
Source={}; Source[SourceBus[0]] = SubPDElement

# Add another source
#SourceBus.append('sourcebus2')
Source = FindSourcebusConnection(Source, SourceBus, PDElementsConnections, ParentBus)

S_Base = Ckt.FindS_Base()
V_Base = Ckt.FindV_Base()
Z_Base = (V_Base**2)*1000/S_Base  # Base Impedance-region 1
I_Base = S_Base/(math.sqrt(3)*V_Base)  # Base current
L_Base = I_Base**2        # square of current base;
PhasesNum = 3  # distribution networks phase numbers
PhaseSequence = ['a', 'b', 'c']

# *************************************************************
# * Loads Data
# *************************************************************
# Load connections
LoadsConnection = FindLoadsConnections(Ckt, S_Base)

# Load to Bus Connection matrix
LoadToBus = FindLoadToBus(AllBuses, LoadsConnection)

# Loads Maximum capacity
Loads_max = FindLoadsCapacity(Ckt, S_Base)

# Modify secondary loads (Not needed for 123 node)
#[AllLoadsModi, LoadToBusModi, LoadsConnectionModi, Loads_maxModi] = ModifySecondaryLoad(CktNum, Ckt, AllLoads,
#                                            ParentBus,PDElementsConnections, LoadsConnection, Loads_max, LoadToBus)

# Loads priority or weight
LoadsWeight = FindLoadsWeight(AllLoads)#, s65a=1.5, s65b=1.5, s65c=1.5, s48=1.5)


# ************************************************************************
# Load Type matrix
# ***********************************************************************
# 1 mean non-dispatchable
# pass load name to function if you want to change it into dispatchable
DispatchableLoads = FindDisptachableLoads(AllLoads, 200, 300)#0, 1177)
LoadType = DefineLoadType(AllLoads, DispatchableLoads)

# ************************************************************************
# Lines
# ***********************************************************************
# Lines Length (The unit length of ckt 123-nodes is kft)
LinesLength = FindLinesLength(Ckt)

# Lines Capacity
LinesSquareCapacity = FindLinesSquareCapacity(Ckt, I_Base)

# Line Resistance (Ohm/units)
LinesResistance = FindLinesResistance(Ckt, AllLines)

# Line Reactance (Ohm/units)
LinesReactance = FindLinesReactance(Ckt, AllLines)

# ************************************************************************
# PD Elements = Lines + Transformers
# ***********************************************************************
# All PD Elements Length
PDElementsLength = FindPDElementLength(AllPDElements, LinesLength)

# All PD Elements Capacity
PDElementsSquareCapacity = FindPDElementsSquareCapacity(AllPDElements, LinesSquareCapacity)

# Integrate all PD Elements Resistance
PDElementsR = FindPDElementsR(AllPDElements, PDElementPhaseMatrix, LinesResistance, PDElementsLength, Z_Base)

# Integrate all PD Elements Reactance
PDElementsX = FindPDElementsX(AllPDElements, PDElementPhaseMatrix, LinesReactance, PDElementsLength, Z_Base)

# PD Elements R_hat resistance:
PDElementsR_hat = CalculatePDElements_R_hat(PDElementsR, PDElementsX, math.sqrt, np)

# PD Elements X_hat resistance:
PDElementsX_hat = CalculatePDElements_X_hat(PDElementsX, PDElementsR, math.sqrt, np)

# PD Elements X_hat resistance:
PDElementsZ_hat = CalculatePDElements_Z_hat(PDElementsR, PDElementsX, math.sqrt)


# *****************************************************************************
# PV
# *****************************************************************************
# define PV matrix Bus connection which first is PV capacity
PVconnection = FindPVConnection(Ckt, S_Base)#, Busm1026706=[333, 1], Busm1047592=[333, 1])
PVBus = FindPVtoBus(AllBuses, PVconnection)

# *****************************************************************************
# Distribution Lines' capacity
# *****************************************************************************
# Switchable lines, all lines name like Line.sw... considers as switchable, also other lines can be added to end of
# this function to be considered as switchable
SwitchableLines = FindSwitchableLines(AllPDElements) #, 'Line.110')

#Normlly open switches can have higher wieght
NormallyCloseSwitches = ['Line.swx8223_48332_sw', 'Line.swv7173_48332_sw', 'Line.swv9287_48332_sw',
                             'Line.swa8735_48332_sw', 'Line.swl5437_48332_sw', 'Line.swln4625713_sw',
                             'Line.swln4641075_sw',
                             'Line.swl5491_48332_sw', 'Line.swln4625696_sw', 'Line.swln4586093_sw',
                             'Line.swl5659_48332_sw',
                             'Line.swln3693186_sw', 'Line.swv7313_48332_sw', 'Line.swa8611_48332_sw',
                             'Line.swln4625876_sw',
                             'Line.swv9111_48332_sw', 'Line.swx8271_48332_sw', 'Line.swv7041_48332_sw',
                             'Line.swl5397_48332_sw', 'Line.swln247171_sw', 'Line.swl5523_48332_sw',
                             'Line.sw2002200004641085_sw', 'Line.swxj171_48332_sw', 'Line.swa8869_48332_sw',
                             'Line.sw2002200004991174_sw', 'Line.swl9191_48332_sw', 'Line.swln0247162_sw',
                             'Line.swln247160_sw', 'Line.swln293471_sw', 'Line.swa8645_48332_sw',
                             'Line.swv9109_48332_sw',
                             'Line.swg9343_48332_sw', 'Line.swl5565_48332_sw', 'Line.swa333_48332_sw',
                             'Line.swln4625680_sw', 'Line.sw2002200004868472_sw', 'Line.swl9407_48332_sw']

SwitchableLineWeight = DefineSwitchWeight(SwitchableLines, NormallyCloseSwitches, 10, 1) #third entry is weight of NC Sw.

# Switches for Lines
Switches_Lines = FindSwitches(SwitchableLines)

# *****************************************************************************
# CapacitorBanks
# *****************************************************************************
CapacitorConnection = FindCapacitorConnection(Ckt, S_Base)
CapacitorBus = FindCapacitorToBus(AllBuses, CapacitorConnection)

# *****************************************************************************
# X_Y coordinates
# *****************************************************************************
XY_Coordinates = FindXY_Coordinates(Ckt)
#print(XY_Coordinates)

# Modify coordinates for switches in IEEE 8500 node system
XY_Coordinates = ModifySwitchesCoordinates(XY_Coordinates, SwitchableLines, PDElementsConnections, ParentBus,
                                               ChildrenBuses)
# *****************************************************************************
#print(AllBuses)

EndNetworkLoadingTime = time.time()



