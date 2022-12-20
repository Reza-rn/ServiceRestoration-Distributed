# This is the Simplified Distribution network reconfiguration program
# In this version, fault is isolated through switches (This version is much more debugged than version 3)
# ***************************************************************************************************************
# Restoration constants
# ***************************************************************************************************************
from OpenDSSCkt import *
from gurobipy import *

# ***************************************************************************************************************
# Optimization constants
# ***************************************************************************************************************
M = 100

# ***************************************************************************************************************
# Restoration constants
# ***************************************************************************************************************
# Restoration Time-----------------------
RestorationTime = SetRestorationTimes(1)

# Source Capacity-----------------------
SourceCapacity = ArrangeSourceCapacity(SourceBus, S_Base, [100000000, 3000, 3000])
USource = CalculateUsource(SourceBus, SourceCapacity, RestorationTime, 1.04)

# Voltages-------------------------
V_std = 1
V_Tolerance = 0.05
V_max = V_std + V_Tolerance
V_min = V_std - V_Tolerance

U_max = V_max * V_max
U_min = V_min * V_min
# ---------------------------------
# Faulted Lines
if CktNum == 4:
    FaultedLines = []#'Line.ln5774453-1', 'Line.ln6047578-1', 'Line.c10_ln5774453-1', 'Line.c3_ln6504020-2']
elif CktNum == 3:
    FaultedLines = []#'Line.ln5774453-1', 'Line.ln6047578-1', 'Line.157132']
elif CktNum == 2:
    FaultedLines = []#'Line.ln5774453-1', 'Line.ln6047578-1']
elif CktNum == 1:
    FaultedLines = []#'Line.l55']
else:
    FaultedLines = []
# 'Line.ln6350541-1']#'Line.ln6413954-1']#'Line.ln6139641-1']

# ***************************************************************************************************************
# Manipulation of Circuit Data
# ***************************************************************************************************************

# ***************************************************************************************************************
# Create a Model
# ***************************************************************************************************************
DSR = Model('DSR')
DSR.Params.timelimit = 3600.0
# ***************************************************************************************************************
# Create Variables
# ***************************************************************************************************************
# Continuous Variables!!!!!!!!!!!!!*****************************************************************************
# Substation Active Power---------------------------------------------------------------------
P_sub = DSR.addVars(PhaseSequence, SourceBus, lb=0, name='P_sub')
# Substation Reactive Power---------------------------------------------------------------------
Q_sub = DSR.addVars(PhaseSequence, SourceBus, lb=0, name='Q_sub')
# Buses Square Voltage---------------------------------------------------------------------
U_Bus = DSR.addVars(PhaseSequence, AllBuses, lb=0, name='U')
# Loads Active power-----------------------------------------------------------------------
P_Load = DSR.addVars(AllLoadsModi, lb=0, name='P_Load')
# Loads Reactive power---------------------------------------------------------------------
Q_Load = DSR.addVars(AllLoadsModi, lb=0, name='Q_Load')
# Distribution lines square current--------------------------------------------------------
l_Line = DSR.addVars(PhaseSequence, AllPDElements, lb=0, name='l_Line')
# Distribution lines Active power----------------------------------------------------------
P_Line = DSR.addVars(PhaseSequence, AllPDElements, lb=-5, name='P_Line')
# Distribution lines Active power----------------------------------------------------------
Q_Line = DSR.addVars(PhaseSequence, AllPDElements, lb=-5, name='Q_Line')
# PV Active power--------------------------------------------------------------------------
P_PV = DSR.addVars(PhaseSequence, PVconnection, lb=0, name='P_PV')
# PV Reactive power--------------------------------------------------------------------------
Q_PV = DSR.addVars(PhaseSequence, PVconnection, lb=-5, name='Q_PV')
# Total Active Power--------------------------------------------------------------------------
P_LoadTotoal = DSR.addVar(lb=0, name='P_LoadTotal')
# Total Generation--------------------------------------------------------------------------
P_GenTotoal = DSR.addVar(lb=0, name='P_GenTotal')

# Binary Variables!!!!!!!!!!!!!**************************************************************************
# Bus binary variables--------------------------------------------------------------------------
x_Bus = DSR.addVars(AllBuses, vtype=GRB.BINARY, name='x_Bus')
# Load binary variables--------------------------------------------------------------------------
x_Load = DSR.addVars(AllLoadsModi, vtype=GRB.BINARY, name='x_Load')
# Line binary variables--------------------------------------------------------------------------
#Because of spaning tree constraints not need to be binary
x_Line = DSR.addVars(Switches_Lines, lb=0, name='x_Line')
# Line flow direction--------------------------------------------------------------------------
FlowDirection = ['d1','d2']
Beta = DSR.addVars(AllPDElementsReduced, FlowDirection, vtype=GRB.BINARY, name='Beta')

# ***************************************************************************************************************
# Objective Function
# ***************************************************************************************************************
# objective functions:
# Restored Loads----------------
obj1 = quicksum(LoadsWeight[j]*P_Load[j] for j in AllLoadsModi)

# Line Energization
obj2 = quicksum(SwitchableLineWeight[j]*x_Line[j] for j in Switches_Lines)

# Loss minimization
# obj3 = -CalculateTotalLoss(RestorationTime, AllPDElements, PhaseSequence, l_Line, PDElementsR)

# multi-objective using weighted sum method
if PVconnection:
    if PVconnection['PVBusm1026706'][1] <= 400/S_Base:
        c1 = 0.9997
    else:
        c1 = 0.997
else:
    c1 = 0.99
c2 = 1 - c1

obj = c1*obj1 + c2*obj2

DSR.setObjective(obj1, GRB.MAXIMIZE)

# multi-objective using Gurobi built-in options
# Primary objective:
'''DSR.setObjectiveN(obj1, 2, 1)
# Alternative, lower priority objectives:
DSR.setObjectiveN(obj2, 1, 0)
DSR.ModelSense = GRB.MAXIMIZE
'''
# ***************************************************************************************************************
# Constraints
# ***************************************************************************************************************
# Load capacity
for j in AllLoadsModi:
    DSR.addConstr(P_Load[j] == x_Load[j]*Loads_maxModi[j][0])
    DSR.addConstr(Q_Load[j] == x_Load[j]*Loads_maxModi[j][1])

# ----------------------------------------------------------------------------------------------------------------
# Total Load
for t in RestorationTime:
    DSR.addConstr(P_LoadTotoal == quicksum(P_Load[j] for j in AllLoadsModi))

# ----------------------------------------------------------------------------------------------------------------
# Total Generation capacity
DSR.addConstr(P_GenTotoal == quicksum(P_sub[i, j] for i in PhaseSequence for j in SourceBus)+
                                quicksum(P_PV[i, j] for i in PhaseSequence for j in PVconnection))

DSR.addConstr(P_LoadTotoal <= P_GenTotoal)

# ----------------------------------------------------------------------------------------------------------------
# Total substation power capacity
for j in SourceBus:
    DSR.addConstr(P_sub['a', j] + P_sub['b', j] + P_sub['c', j] <= SourceCapacity[j])
'''
# ----------------------------------------------------------------------------------------------------------------
# Power Factor or reactive power constraints
for t in RestorationTime:
    for j in AllLoadsModi:
        DSR.addConstr(Q_Load[j, t] <= 0.6197*P_Load[j, t])

# pick-up all reactive power loads at the last step
for j in AllLoadsModi:
    DSR.addConstr(Q_Load[j, RestorationTime[-1]] == Loads_maxModi[j][1])
'''
# ----------------------------------------------------------------------------------------------------------------
# Voltage
# Un-enegized buses cannot be detected to make their voltage zero!!!!!!!
for j in AllBuses:
    if j != SourceBus:
        for i in PhaseSequence:
            if BusesPhaseMatrix[j][PhaseSequence.index(i)] == 1:
                DSR.addConstr(U_Bus[i, j] <= U_max * x_Bus[j])
                DSR.addConstr(x_Bus[j] * U_min <= U_Bus[i, j])

            elif BusesPhaseMatrix[j][PhaseSequence.index(i)] == 0:
                DSR.addConstr(U_Bus[i, j] == 0)

# Substation Voltage
for j in SourceBus:
    for i in PhaseSequence:
        DSR.addConstr(U_Bus[i, j] == USource[j]['t1'])

# ----------------------------------------------------------------------------------------------------------------
# Branch Power capacity (can be covered by power balance equations)
for j in AllPDElements:
    #Replace Removed PD Elements with representatives for lines binaries
    TempLine = (j + '.')[:-1]
    if j in RemovedPDElements:
        TempLine = RemovedPDElements[j]

    for i in PhaseSequence:
        if PDElementPhaseMatrix[j][PhaseSequence.index(i)] == 0:
            DSR.addConstr(P_Line[i, j] == 0)
            DSR.addConstr(Q_Line[i, j] == 0)
        else:
            if j in Switches_Lines:
                DSR.addConstr(P_Line[i, j] <= 3*x_Line[TempLine])
                DSR.addConstr(-3*x_Line[TempLine] <= P_Line[i, j])
                DSR.addConstr(Q_Line[i, j] <= 3*x_Line[TempLine])
                DSR.addConstr(-3*x_Line[TempLine] <= Q_Line[i, j])
            else:
                DSR.addConstr(P_Line[i, j] <= 3)
                DSR.addConstr(-3 <= P_Line[i, j])
                DSR.addConstr(Q_Line[i, j] <= 3)
                DSR.addConstr(-3 <= Q_Line[i, j])

# ----------------------------------------------------------------------------------------------------------------
# Generation Constraints
# PV Constraints*****************************************************************
# Maximum PV inverter capacity-------------------
for j, k in PVconnection.items():
    for i in PhaseSequence:
        DSR.addConstr(P_PV[i, j] <= x_Bus[k[0]]*k[1])
        #Reactive power (add %20 more capacity to inverter to compensate reactive power)
        DSR.addConstr(Q_PV[i, j] <= 0.2*k[1])
        DSR.addConstr(-0.2*k[1] <= Q_PV[i, j])

# ----------------------------------------------------------------------------------------------------------------
# Power flow constraints:
# Power Balance--------------------------------------------------
for j in AllPDElements:
    for i in PhaseSequence:
        # Active Power-----------------------------------
        DSR.addConstr(P_Line[i, j] == (ActivePowerBalance(j, i, P_Load, P_Line, P_PV, P_sub,
                        PDElementsConnections, PDElementPhaseMatrix, LoadToBusModi, LoadsConnectionModi, ChildrenLines,
                                            PVBus, SourceBus, Source))*PDElementPhaseMatrix[j][PhaseSequence.index(i)])
        # Reactive Power-----------------------------------
        DSR.addConstr(Q_Line[i, j] == (ReactivePowerBalance(j, i, Q_Load, Q_Line, Q_PV, Q_sub,
                   PDElementsConnections, PDElementPhaseMatrix, LoadToBusModi, LoadsConnectionModi, ChildrenLines,
                                            PVBus, SourceBus, Source))*PDElementPhaseMatrix[j][PhaseSequence.index(i)])

# Voltage equation-------------------------------------------------
for j in AllPDElements:
    # Replace Removed PD Elements with representatives for lines binaries
    TempLine = (j + '.')[:-1]
    if j in RemovedPDElements:
        TempLine = RemovedPDElements[j]
    for i in PhaseSequence:
        if PDElementPhaseMatrix[j][PhaseSequence.index(i)] == 1:
            if j in Switches_Lines:
                Bus_Temp = PDElementsConnections[j][1]
                ParentBus_Temp = PDElementsConnections[j][0]
                if BusesPhaseMatrix[Bus_Temp][PhaseSequence.index(i)] == 1:
                    DSR.addConstr(U_Bus[i, Bus_Temp] - (VoltageRHS(i, j, ParentBus_Temp, U_Bus, P_Line, Q_Line,
                                                                        PDElementsR_hat, PDElementsX_hat))
                                                                        <= (1 - x_Line[TempLine]) * M)

                    DSR.addConstr(U_Bus[i, Bus_Temp] - (VoltageRHS(i, j, ParentBus_Temp, U_Bus, P_Line, Q_Line,
                                                                        PDElementsR_hat, PDElementsX_hat))
                                                                         >= -(1 - x_Line[TempLine]) * M)

            else:
                Bus_Temp = PDElementsConnections[j][1]
                ParentBus_Temp = PDElementsConnections[j][0]
                if BusesPhaseMatrix[Bus_Temp][PhaseSequence.index(i)] == 1:
                    DSR.addConstr(U_Bus[i, Bus_Temp] == (VoltageRHS(i, j, ParentBus_Temp, U_Bus, P_Line, Q_Line,
                                                                   PDElementsR_hat, PDElementsX_hat)))

# ***************************************************************************************************************
# Connectivity and Sequencing Constraints:*********************************************************

# Lines connectivity ====================================
# Out of Service Lines (Faulted Lines)----------
for j in FaultedLines:
    Bus_Temp = PDElementsConnections[j][1]
    ParentBus_Temp = PDElementsConnections[j][0]
    for i in PhaseSequence:
        DSR.addConstr(P_Line[i, j] == 0)
        DSR.addConstr(Q_Line[i, j] == 0)
        DSR.addConstr(U_Bus[i, Bus_Temp] == 0)
        DSR.addConstr(U_Bus[i, ParentBus_Temp] == 0)


# Spanning tree constraints=========================================
# Substation bus is not parent bus of any bus
for j in SourceBus:
    if SourceBus.index(j) == 0:
        DSR.addConstr(Beta[Source[j], 'd2'] == 0)
    else:
        DSR.addConstr(Beta[Source[j], 'd1'] == 0)

# A line represents a single relation
for j in AllPDElementsReduced:
    if j not in FaultedLines:
        if j in Switches_Lines:
            DSR.addConstr(Beta[j, 'd1'] + Beta[j, 'd2'] == x_Line[j])
        else:
            DSR.addConstr(Beta[j, 'd1'] + Beta[j, 'd2'] == 1)

# Every bus has one or less parent bus
for j in AllBuses:
    DSR.addConstr(FindParentFlow(j, Beta, PDElementsConnectionsReduced, BusNeighborBuses) <= 1)

# ***************************************************************************************************************
# Perform Optimization
# ***************************************************************************************************************
DSR.optimize()

# ***************************************************************************************************************
# Optimization Results
# *****************************************************************************************************
# Write optimization problem in text file-----------------------------------
DSR.write('out.lp')

# Draw circuit---------------------------------------------------------------
PlotCircuit(XY_Coordinates, AllBuses, AllPDElementsReduced, PDElementsConnectionsReduced, DSR, plt,
             PVBus, SwitchableLines, AllLoadsModi, LoadsConnectionModi, PDElementPhaseMatrix, FaultedLines, SourceBus)

'''for v in DSR.getVars():
    if v.varName.split('[')[0] == 'x_Load':
        print('%s %g' % (v.varName, v.x))

for v in DSR.getVars():
    if v.varName.split('[')[0] == 'x_Line':
        if SwitchableLines[v.varName.split('[')[1].split(']')[0]] == 1:
            print('%s %g' % (v.varName, v.x))

for v in DSR.getVars():
    if v.varName.split('[')[0] == 'P_Line':
        if v.varName.split('[')[1].split(',')[1].split(']')[0] == Source[SourceBus[0]]:
            print('%s %g' % (v.varName, v.x))

'''
for v in DSR.getVars():
    if v.varName.split('[')[0] == 'P_LoadTotal':
        print('%s %g' % (v.varName, v.x))
print(S_Base)
'''for v in DSR.getVars():
    if v.varName.split('[')[0] == 'U':
        print('%s %g' % (v.varName, v.x))'''


'''for v in DSR.getVars():
    if v.varName.split('[')[0] == 'U':
        if v.varName.split('[')[1].split(',')[1].split(']')[0] == 'e206209':
            print('%s %g' % (v.varName, v.x))

for v in DSR.getVars():
    if v.varName.split('[')[0] == 'U':
        if v.varName.split('[')[1].split(',')[1].split(']')[0] == 'l3235258':
            print('%s %g' % (v.varName, v.x))'''

for v in DSR.getVars():
    if v.varName.split('[')[0] == 'U':
        if v.varName.split('[')[1].split(',')[1].split(']')[0] == 'c2_c8_node':
            print('%s %g' % (v.varName, v.x))
        elif v.varName.split('[')[1].split(',')[1].split(']')[0] == 'masternode':
            print('%s %g' % (v.varName, v.x))
        elif v.varName.split('[')[1].split(',')[1].split(']')[0] == 'c0_c6_node':
            print('%s %g' % (v.varName, v.x))


for v in DSR.getVars():
    if v.varName.split('[')[0] == 'x_Line':
        if v.x ==0:
            print('%s %g' % (v.varName, v.x))

EndTime = time.time()
# Voltage Constraint
'''
for v in DSR.getVars():
    if v.varName.split('[')[0] == 'P_Line':
        print('%s %g' % (v.varName, v.x))

for v in DSR.getVars():
    if v.varName.split('[')[0] == 'P_Load':
        print('%s %g' % (v.varName, v.x))

for v in DSR.getVars():
    if v.varName.split('[')[0] == 'P_LoadTotal':
        print('%s %g' % (v.varName, v.x))


for v in DSR.getVars():
    if v.varName.split('[')[0] == 'U':
        print('%s %g' % (v.varName, np.sqrt(v.x)))
'''

print(EndTime-StartTime)
print(EndNetworkLoadingTime-StartTime)



