#*******************************************************************************************************************
__author__ = "Reza Roofegari nejad"                                                                          #******
__email__ = "r.roofegari@knights.ucf.edu"                                                                    #******
__website__ = "http://www.ece.ucf.edu/~rezarn/"                                                              #******
__copyright__ = "Copyright 2019, Scalable/Secure Cooperative Algorithms and Framework for Extremely-high " \
                "Penetration Solar Integration (SolarExPert), ENERGISE Project"                              #******
__ProjectWebsite__ = "https://www.cs.ucf.edu/~qu/MA-OpenDSS.php"                                             #******
#*******************************************************************************************************************
# 1. This is Distibuted Service restoration framework. In this model heuristic ADMM is used.
# 2. In this version, fault is isolated through switches
# 3. In this version fisrt the problem is relaxed and then when the risidual is less than a minimum the projection terms
# are included in the problem. in fact, the first relaxed version worked as a warm start for the rest of the problem.
# ***************************************************************************************************************
# Restoration constants
# ***************************************************************************************************************
# from Clustering import *
import Clustering as Circuit
from Functions import *
import gurobipy as Solver
from math import sqrt
from DSR_Function import *
import pickle, time
from DSR_Polish_Function import *
import matplotlib.pyplot as plt
#plt.rcParams['figure.dpi'] = 400 #200

# ***************************************************************************************************************
# Optimization constants
# ***************************************************************************************************************
M = 100

# ***************************************************************************************************************
# Restoration constants
# ***************************************************************************************************************
# Restoration Time-----------------------
RestorationTime = SetRestorationTimes(3)

# Source Capacity-----------------------
#SourceCapacity = ArrangeSourceCapacity(Circuit.SourceBus, Circuit.S_Base, [10000, 3000, 3000])
#SourceCapacity = {'sourcebus':[1000]}
SourceCapacity = {'sourcebus':[1000,5000,20000]}

SourceCapacity = FindPerUnitofSourceCapacity(SourceCapacity, Circuit.S_Base)
USource = CalculateUsource(SourceCapacity, RestorationTime, sourcebus=1)

# Voltages-------------------------

V_std = 1
V_Tolerance = 0.05
V_max = V_std + V_Tolerance
V_min = V_std - V_Tolerance

U_max = V_max * V_max
U_min = V_min * V_min
# ---------------------------------
# Faulted Lines
#FaultedLines = ['Line.ln5594239-1']
FaultedLines = ['Line.ln5503576-1']#'Line.ln5774453-1', 'Line.ln6047578-1']
IsolatingSwitches = ['Line.swa8645_48332_sw','Line.swa8611_48332_sw']#'Line.sw4', 'Line.sw5', 'Line.sw8']
Oscillating_Vars = []#'f739845','293471','f739845','293471','f739841','l0247160','f739842','l0247162','e182744','l2801895','293471']
# ***************************************************************************************************************
# ADMM Parameters and Constants and Variables
# ***************************************************************************************************************
rho = 1
rho_Polish = 10

#-----------------
# Dynamic penalty parameter constants
Dynamic_rho = 1 # Set equal to 1 if want to have dynamic rho
Dynamic_rho_Polish = 1
Thau_inc = 2
Thau_dec = 2
mu_penalty = 10

#-----------------------------------------------
iteration = 0
inner_iteration = 1 #projection iteration count
iteration_max = 1500

# Proximal operator parameters------------------------------
Enable_Prox = 0  # This parameter determines if drive phase should be existed or directly after relax go into projection
t_Prox_Max = 10000
t_prox_iteraion_stop = 200 # The number of iteration that prox should be stopped before projection
t_prox_iteraion_min = 200
t_Prox = [0, 0]
c1 = 0.7 # constant for showing growth of t_prox (used for normal operation)
c2 = 20

# Projected ADMM trigger, if this variable becomes 1 the projected ADMM starts
Projected_ADMM_trigger = 0

# polish phase --------------------------------------------
Polish_Trigger = 1
Polish_Maximum_iteration = 1000
Polish_status = 0 # This shows if we are in the polish state
Just_Polish = 1 # if this become 1 then we will have just polish phase without reconfiguration
#-----------------------------------------------
# Residual
PrimalResidual_All = [1]
DualResidual_All = [1]

PrimalResidual = 5
DualResidual = 5

Binary_Residual_Status = 0
PrimalBinaryResidual = 1
DualBinaryResidual = 1
# Projection Residuals
PrimalResidual_Projection = 5
DualResidual_Projection = 5

# Residual Bound
#Residual_bound = 0.005*sqrt(Circuit.NumClusters)
Residual_bound = 0.0005*sqrt(Circuit.NumClusters)
Residual_bound_Projection = 0.001*sqrt(Circuit.NumClusters)
#Residual_bound_Projection = 0.0001*sqrt(Circuit.NumClusters)
Residual_bound_Prox = 0.0001
Residual_bound_Polish = 0.0005*sqrt(Circuit.NumClusters)


#-----------------------------------------------
# build Lagrange multipliers matrix
mu = Build_mu(Circuit, RestorationTime)
#-----------------------------------------------
# X and Z variables
X = Build_X(Circuit, RestorationTime)  # Continous Variables
Y = Build_Y(Circuit, RestorationTime)  # Relaxed Binary Variables

# Consensus and projected variables
X_bar = Build_X_bar(Circuit, mu, RestorationTime, Circuit.PhaseSequence)
if not (Just_Polish):
    Z = Build_Z(Circuit, mu, RestorationTime, 0)  # The fourth entry is the initial value of switches and etc
else:  # This section is used when using Just Polish phase
    Z = Build_Z(Circuit, mu, RestorationTime, 1)
    # Manually change status of some switches:
    Opened_Switches = ['Line.swwf586_48332_sw','Line.swwg127_48332_sw','Line.swwd701_48332_sw','Line.swwf856_48332_sw',
                       'Line.swv7995_48332_sw','Line.swa8645_48332_sw','Line.swa8611_48332_sw']#,'Line.swtieline']
    Z = Manipulate_Z(Z, Circuit, RestorationTime, Opened_Switches)

#-------------------------------------------
# Some Saved items during evolution of ADMM
SavedItems = Build_VariableSave(RestorationTime, X, Y, Z)

# ***************************************************************************************************************
# ADMM loop starts here
# ***************************************************************************************************************
StartDDSR_Time = time.time()
if not(Just_Polish):
    while ((((PrimalResidual > Residual_bound_Projection) and (DualResidual > Residual_bound_Projection)) or
           (not Projected_ADMM_trigger) ) and (iteration <= iteration_max)):

        iteration = iteration + 1
        IterationIndex = 'iteration{}'.format(iteration)
        # ----------------------------------------------------------------------------------------------------
        # First step of consensus ADMM (Each agent solves its own problem and update data)------------------------------
        for cluster in Circuit.Clusters:
            DSRCluster = DSR_Cluster(Circuit, Solver, RestorationTime, SourceCapacity, USource, U_max, U_min, FaultedLines,
                                     M, cluster, mu, X, Y, X_bar, Z, rho, iteration, IterationIndex, Projected_ADMM_trigger,
                                    IsolatingSwitches)

        #------------------------------------------------------------------------------------------------------
        # Second step of consensus ADMM (Projection for Z matrix and consensus update for X_bar)

        # Consensus update for X_bar matrix-------------
        X_bar = Update_X_bar(X_bar, X, mu, iteration, Circuit.ClustersNeighborBuses, Circuit.BusesInClusters,
                             Circuit.ClustersBoundaryPDElements, Circuit.ClustersNeighborRelations, Circuit.PhaseSequence)

        # Projection update for Z matrix-------------
        Z = Update_Z(Z, Y, mu, iteration, Binary_Projection, Circuit.ClustersNeighborRelations, Circuit.BusesInClusters,
                     Circuit.ClustersBoundaryPDElements, Circuit.ClustersBoundarySwitchableLines,
                     Circuit.ClustersNeighborBuses, Projected_ADMM_trigger, Prox_Oper, t_Prox[iteration])

        # ------------------------------------------------------------------------------------------------------
        # Third step of consensus ADMM (Updating Lagrange multipliers)
        mu = Update_mu(mu, X, X_bar, Y, Z, iteration, Polish_status)

        # ------------------------------------------------------------------------------------------------------
        # Calculate Residuals

        # Primal Residual--------------------
        PrimalResidual = Calculate_PrimalResidual(X, Y, X_bar, Z, iteration, mu, sqrt, Projected_ADMM_trigger,
                                Circuit.ClustersBoundarySwitchableLines, Circuit.ClustersBoundaryPDElements,
                                                  Polish_status, 0, Oscillating_Vars)
        PrimalResidual_All.append(PrimalResidual)

        # Dual Residual----------------------
        DualResidual = Calculate_DualResidual(X_bar, Z, iteration, mu, sqrt, rho, Projected_ADMM_trigger,
                              Circuit.ClustersBoundarySwitchableLines, Circuit.ClustersBoundaryPDElements,
                                              Polish_status, 0, Oscillating_Vars)
        DualResidual_All.append(DualResidual)

        # Primal Binary Residual--------------------
        #if Binary_Residual_Status:
        #    PrimalBinaryResidual = Calculate_PrimalResidual(X, Y, X_bar, Z, iteration, mu, sqrt, Projected_ADMM_trigger,
        #                            Circuit.ClustersBoundarySwitchableLines, Circuit.ClustersBoundaryPDElements,
        #                                                Polish_status, 0)
        # Dual Binary Residual----------------------
        #    DualBinaryResidual = Calculate_DualResidual(X_bar, Z, iteration, mu, sqrt, rho, Projected_ADMM_trigger,
        #                                      Circuit.ClustersBoundarySwitchableLines, Circuit.ClustersBoundaryPDElements,
        #                                      Polish_status, 0)



        # if neccessary change Process into projected ADMM
        if not Projected_ADMM_trigger:
            if (((iteration <= iteration_max - t_prox_iteraion_stop) and (t_Prox[iteration] <= t_Prox_Max))
                 and ((PrimalResidual >= Residual_bound_Prox) or (inner_iteration <= t_prox_iteraion_min))):
                if t_Prox[iteration] == 0:
                    if (PrimalResidual <= Residual_bound) and (DualResidual <= Residual_bound):
                        if Enable_Prox:
                            t_Prox.append(0.01)
                            PrimalResidual = 1
                            DualResidual = 1
                            Binary_Residual_Status = 1
                        else:
                            Projected_ADMM_trigger = 1
                            PrimalResidual = 1
                            DualResidual = 1
                            t_Prox.append(t_Prox[iteration])
                            Binary_Residual_Status = 0

                    else:
                        t_Prox.append(t_Prox[iteration])

                else:
                    if (PrimalResidual != 0)and(DualResidual != 0):
                        t_Prox.append(t_Prox[iteration] + c1*((1/PrimalResidual)+(1/DualResidual)) )
                    elif PrimalResidual == 0:
                        t_Prox.append(t_Prox[iteration] + c1*(1/DualResidual)*10)
                    elif DualResidual == 0:
                        t_Prox.append(t_Prox[iteration] + c1*(1/PrimalResidual)*10)
                    else:
                        t_Prox.append(t_Prox[iteration] + 20)
                    inner_iteration += 1


            else:
                if Projected_ADMM_trigger == 0:
                    Projected_ADMM_trigger = 1
                    PrimalResidual = 0.1
                    DualResidual = 0.1
                    t_Prox.append(t_Prox[iteration])
                    Binary_Residual_Status = 0
                else:
                    t_Prox.append(t_Prox[iteration])
        else:
            t_Prox.append(t_Prox[iteration])

        # Save some costumized variables
        SavedItems = ADMM_VariableSave(SavedItems, X, Y, Z, iteration, 1, Circuit.LoadType, Circuit.Loads_max,
                                       Circuit.S_Base, Polish_status, RestorationTime, Circuit.LoadsConnection,
                                       Circuit.BusesInClusters)

        # Dynamic Penalty Parameter of ADMM-------
        if Dynamic_rho == 1:
            if PrimalResidual > (mu_penalty*DualResidual):
                rho = rho * Thau_inc
            elif DualResidual > (mu_penalty*PrimalResidual):
                rho = rho / Thau_dec

            if rho < 0.4:
                rho = 0.1
            elif rho > 4:
                rho = 4

# ***************************************************************************************************
# Start polish phase*********************************************************************************
# ***************************************************************************************************
Last_iter = iteration + 0  # Drive phase last iteration
if Polish_Trigger:
    Polish_status = 1 # This shows if we are in the polish state
    Polish_iteration = 0
    PrimalResidual_Polish = 1
    DualResidual_Polish = 1
    while (((PrimalResidual_Polish > Residual_bound_Polish) and (DualResidual_Polish > Residual_bound_Polish)) and
        (Polish_iteration <= Polish_Maximum_iteration)):

        iteration = iteration + 1
        IterationIndex = 'iteration{}'.format(iteration)
        Polish_iteration += 1

        # ----------------------------------------------------------------------------------------------------
        # First step of Polish consensus ADMM (Each agent solves its own problem and update data)------------------------------
        for cluster in Circuit.Clusters:
            DSRCluster = DSR_Polish_Cluster(Circuit, Solver, RestorationTime, SourceCapacity, USource, U_max, U_min,
                                FaultedLines, M, cluster, mu, X, X_bar, Z, rho_Polish, iteration, IterationIndex,
                                            Projected_ADMM_trigger, IsolatingSwitches, Last_iter)

        # ------------------------------------------------------------------------------------------------------
        # Second step of Polish consensus ADMM(update consensus variable of X_bar)
        X_bar = Update_X_bar(X_bar, X, mu, iteration, Circuit.ClustersNeighborBuses, Circuit.BusesInClusters,
                             Circuit.ClustersBoundaryPDElements, Circuit.ClustersNeighborRelations, Circuit.PhaseSequence)

        # ------------------------------------------------------------------------------------------------------
        # Third step of Polish consensus ADMM (Updating Lagrange multipliers)
        mu = Update_mu(mu, X, X_bar, Y, Z, iteration, Polish_status)

        # ------------------------------------------------------------------------------------------------------
        # Calculate Residuals

        # Primal Residual--------------------
        PrimalResidual_Polish = Calculate_PrimalResidual(X, Y, X_bar, Z, iteration, mu, sqrt, Projected_ADMM_trigger,
                                Circuit.ClustersBoundarySwitchableLines, Circuit.ClustersBoundaryPDElements,
                                                         Polish_status, 0, Oscillating_Vars)
        PrimalResidual_All.append(PrimalResidual_Polish)

        # Dual Residual----------------------
        DualResidual_Polish = Calculate_DualResidual(X_bar, Z, iteration, mu, sqrt, rho_Polish, Projected_ADMM_trigger,
                            Circuit.ClustersBoundarySwitchableLines, Circuit.ClustersBoundaryPDElements, Polish_status,
                                                     0, Oscillating_Vars)
        DualResidual_All.append(DualResidual_Polish)

        # Save some costumized variables
        SavedItems = ADMM_VariableSave(SavedItems, X, Y, Z, iteration, Last_iter, Circuit.LoadType, Circuit.Loads_max,
                                       Circuit.S_Base, Polish_status, RestorationTime,Circuit.LoadsConnection,
                                       Circuit.BusesInClusters)

        # Dynamic Penalty Parameter of ADMM-------
        if Dynamic_rho_Polish == 1:
            if PrimalResidual_Polish > (mu_penalty * DualResidual_Polish):
                rho_Polish = rho_Polish * Thau_inc
            elif DualResidual_Polish > (mu_penalty * PrimalResidual_Polish):
                rho_Polish = rho_Polish / Thau_dec

            if rho_Polish < 2:
                rho_Polish = 2
            elif rho_Polish > 10:
                rho_Polish = 10


EndDDSR_Time = time.time()
# ***************************************************************************************************************
# Optimization Results
# ***************************************************************************************************************

# residual plot
Plot_Residuals(PrimalResidual_All, DualResidual_All, iteration, plt)

# t-Prox plot
Plot_t_Prox(t_Prox[1:], Last_iter, plt, c1)

# Total load plot
P_TotalLoad = Plot_TotalLoad(X, Z, iteration, Last_iter, plt, Circuit.np, RestorationTime, Circuit.S_Base, Circuit.Loads_max, Circuit.LoadType)

#Circuit plot
Plot_Circuit(X, Y, X_bar, Z, plt, Circuit, RestorationTime[-1], iteration, Last_iter, FaultedLines)

# this function is used to find that points of problem
#Find_Points_Primal(X, Y, X_bar, Z, iteration, mu, sqrt, Projected_ADMM_trigger, Last_iter,Circuit. ClustersBoundarySwitchableLines, Circuit.ClustersBoundaryPDElements)

#Find_Points_Dual(X_bar, Z, iteration, mu, sqrt, rho, Projected_ADMM_trigger, Last_iter, Circuit. ClustersBoundarySwitchableLines, Circuit.ClustersBoundaryPDElements)


#************************************************************
#print total time for DDSR procedure
print('DDSR_TotalTime = {}'.format(EndDDSR_Time-StartDDSR_Time))

# Save results into a matlab file*************************************************************************************
#WorkSpace = {'X':X, 'Y':Y, 'Z':Z, 'X_bar':X_bar, 'mu':mu, 'iteration':iteration, 'Last_iter':Last_iter,
#             'Total_Load':P_TotalLoad, 'PrimalResidual_All':PrimalResidual_All, 'DualResidual_All':DualResidual_All}

#pickle_out = open("WorkSpace1.pickle","wb")
#pickle.dump(WorkSpace, pickle_out)
#pickle_out.close()



#plt.rcParams['figure.dpi'] = 400 #200

