# This function is used to polish DSR problem results from previous phase for each cluster in a Distribution Network

def DSR_Polish_Cluster(Circuit, Solver ,RestorationTime, SourceCapacity, USource, U_max, U_min, FaultedLines, M, ClusterIndex,
                mu, X, X_bar, Z, rho_Polish, iteration, IterationIndex, Projected_ADMM_trigger, IsolatingSwitches, Last_iter):
    # ***************************************************************************************************************
    # Create a Model
    # ***************************************************************************************************************
    DSR = Solver.Model('DSR')
    DSR.Params.timelimit = 3600.0;
    #DSR.Params.OutputFlag = 0;
    # ***************************************************************************************************************
    # Create Variables
    # ***************************************************************************************************************
    # Continuous Variables!!!!!!!!!!!!!*****************************************************************************
    # Loads Active power---------------------------------------------------------------------------------------------
    P_Load = DSR.addVars(Circuit.ClustersLoads[ClusterIndex], RestorationTime, lb=0, name='P_Load')
    # Loads Reactive power-------------------------------------------------------------------------------------------
    Q_Load = DSR.addVars(Circuit.ClustersLoads[ClusterIndex], RestorationTime, lb=0, name='Q_Load')
    # Buses Square Voltage-------------------------------------------------------------------------------------------
    U_Bus = DSR.addVars(Circuit.PhaseSequence, Circuit.ClusterBusVariables[ClusterIndex], RestorationTime, lb=0, name='Voltage')
    # Distribution lines square current------------------------------------------------------------------------------
    #l_Line = DSR.addVars(Circuit.PhaseSequence, Circuit.ClusterCurrentVariables[ClusterIndex], RestorationTime, lb=0, name='l_Line')
    # Distribution lines Active power--------------------------------------------------------------------------------
    P_Line = DSR.addVars(Circuit.PhaseSequence, Circuit.ClustersPDElements[ClusterIndex], RestorationTime, lb=-5, name='P_Line')
    # Distribution lines Active power--------------------------------------------------------------------------------
    Q_Line = DSR.addVars(Circuit.PhaseSequence, Circuit.ClustersPDElements[ClusterIndex], RestorationTime, lb=-5, name='Q_Line')
    # Substation Active Power----------------------------------------------------------------------------------------
    P_sub = DSR.addVars(Circuit.PhaseSequence, Circuit.ClustersSourceBus[ClusterIndex], RestorationTime, lb=0, name='P_sub')
    # Substation Reactive Power--------------------------------------------------------------------------------------
    Q_sub = DSR.addVars(Circuit.PhaseSequence, Circuit.ClustersSourceBus[ClusterIndex], RestorationTime, lb=0, name='Q_sub')
    # Capacitor banks------------------------------------------------------------------------------------------------
    Q_cap = DSR.addVars(Circuit.PhaseSequence, Circuit.ClustersCBs[ClusterIndex], RestorationTime, lb=0, name='Q_cap')
    # PV Active power--------------------------------------------------------------------------
    P_PV = DSR.addVars(Circuit.PhaseSequence, Circuit.ClustersPVs[ClusterIndex], RestorationTime, lb=0, name='P_PV')
    # PV Reactive power--------------------------------------------------------------------------
    Q_PV = DSR.addVars(Circuit.PhaseSequence, Circuit.ClustersPVs[ClusterIndex], RestorationTime, lb=-5, name='Q_PV')
    # Total Active Power--------------------------------------------------------------------------
    P_LoadTotoal = DSR.addVars(RestorationTime, lb=0, name='P_LoadTotal')
    # Total Generation--------------------------------------------------------------------------
    P_GenTotoal = DSR.addVars(RestorationTime, lb=0, name='P_GenTotal')

    # Binary Variables!!!!!!!!!!!!!**************************************************************************
    # Loads binary variables--------------------------------------------------------------------------
    x_Load = DSR.addVars(Circuit.ClustersLoadsNonDispatchable[ClusterIndex], RestorationTime, vtype=Solver.GRB.BINARY, name='x_Load')
    # ***************************************************************************************************************
    # Objective Function
    # ***************************************************************************************************************
    # objective functions:
    # Restored Loads----------------
    obj1 = Solver.quicksum(Circuit.LoadsWeight[j] * P_Load[j, t] for j in Circuit.ClustersLoads[ClusterIndex]
                           for t in RestorationTime)

    # Line Energization-------------
    #obj2 = Solver.quicksum(Circuit.SwitchableLineWeight[j] * Z[Last_iter][t][ClusterIndex]['x_Switch'][j]
    #                       for j in Circuit.ClustersSwitchableLines[ClusterIndex] for t in RestorationTime)

    # Loss minimization-----------------
    # By current****
    #obj3 = -Circuit.CalculateTotalLossByCurrent(RestorationTime, Circuit.ClusterCurrentVariables[ClusterIndex],
     #                                           Circuit.PhaseSequence, l_Line, Circuit.PDElementsR)

    # By Power*****
    #obj3 = Circuit.CalculateTotalLossByPower(RestorationTime, Circuit.ClustersPDElements[ClusterIndex],
    #                    Circuit.PhaseSequence, P_Line, Q_Line, Circuit.PDElementsR, USource[Circuit.SourceBus[0]], Solver)

    # multi-objective using weighted sum method
    c1 = 1
    c2 = 1-c1
    c3 = 1 - c1 - c2
    #obj_f = -c1*obj1 - c2*obj2 + c3*obj3
    #obj_f = -c1*obj1 - c2*obj2
    obj_f = -c1*obj1

    # Main objective function
    obj = Solver.QuadExpr()
    obj.add(obj_f)

    # ADMM terms of the objective function------------------------------------
    obj_Voltage = Solver.QuadExpr()
    obj_P_Line = Solver.QuadExpr()
    obj_Q_Line = Solver.QuadExpr()

    for t in RestorationTime:
        # Voltage---------------------------------------------------------------
        # Voltage Augmented terms are same in both relaxed and projected versions
        for bus in mu[0][RestorationTime[0]][ClusterIndex]['Voltage']:
            for phase in Circuit.PhaseSequence:
                obj_Voltage.add((U_Bus[phase, bus, t] - X_bar[iteration-1][t][ClusterIndex]['Voltage'][bus][phase] +
                                      mu[iteration-1][t][ClusterIndex]['Voltage'][bus][phase]) *
                                (U_Bus[phase, bus, t] - X_bar[iteration - 1][t][ClusterIndex]['Voltage'][bus][phase] +
                                      mu[iteration - 1][t][ClusterIndex]['Voltage'][bus][phase]))

        # Add to the objective function
        obj.add(obj_Voltage, rho_Polish/2)

        # Line power-------------------------------------------------------------
        # Line power Augmented terms are same in both relaxed and projected versions
        for PDElement in Circuit.ClustersBoundaryPDElements[ClusterIndex]:
            for phase in Circuit.PhaseSequence:
                # P_Line
                obj_P_Line.add( (P_Line[phase, PDElement, t] - X_bar[iteration-1][t][ClusterIndex]
                    ['P_Line'][PDElement][phase] + mu[iteration-1][t][ClusterIndex]['P_Line'][PDElement][phase]) *
                    (P_Line[phase, PDElement, t] - X_bar[iteration - 1][t][ClusterIndex]
                    ['P_Line'][PDElement][phase] + mu[iteration - 1][t][ClusterIndex]['P_Line'][PDElement][phase]) )

                # Q_Line
                obj_Q_Line.add( (Q_Line[phase, PDElement, t] - X_bar[iteration-1][t][ClusterIndex]
                    ['Q_Line'][PDElement][phase] + mu[iteration-1][t][ClusterIndex]['Q_Line'][PDElement][phase]) *
                    (Q_Line[phase, PDElement, t] - X_bar[iteration - 1][t][ClusterIndex]
                    ['Q_Line'][PDElement][phase] + mu[iteration - 1][t][ClusterIndex]['Q_Line'][PDElement][phase]) )

        # Add to the objective function
        obj.add(obj_P_Line, rho_Polish/2)
        obj.add(obj_Q_Line, rho_Polish/2)


    # Optimize the objective function------------------------------------------
    DSR.setObjective(obj, Solver.GRB.MINIMIZE)

    # ***************************************************************************************************************
    # Constraints
    # ***************************************************************************************************************
    # Load capacity
    for t in RestorationTime:
        for load in Circuit.ClustersLoads[ClusterIndex]:
            if load in Circuit.ClustersLoadsNonDispatchable[ClusterIndex]:
                DSR.addConstr(P_Load[load, t] == x_Load[load, t]*Circuit.Loads_max[load][0])
                DSR.addConstr(Q_Load[load, t] == x_Load[load, t]*Circuit.Loads_max[load][1])

            else:
                DSR.addConstr(P_Load[load, t] <= Circuit.Loads_max[load][0])
                DSR.addConstr(Q_Load[load, t] <= Circuit.Loads_max[load][1])

    # Load energization-----------------------
    '''for t in RestorationTime:
        for load in Circuit.ClustersLoads[ClusterIndex]:
            if load in Circuit.ClustersLoadsNonDispatchable[ClusterIndex]:
                DSR.addConstr(x_Load[load, t] <= x_Bus[Circuit.LoadsConnection[load][0], t])'''

    # ----------------------------------------------------------------------------------------------------------------
    # Sequencing constraint (Once a Load energized, it should not be shedded!)
    for t in RestorationTime:
        if t != RestorationTime[0]:
            for load in Circuit.ClustersLoads[ClusterIndex]:
                DSR.addQConstr(P_Load[load, t] >= P_Load[load, RestorationTime[RestorationTime.index(t) - 1]])

    # ----------------------------------------------------------------------------------------------------------------
    # Load and Available power capacity of each cluster
    # Total Load of the cluster------------
    for t in RestorationTime:
        DSR.addConstr(P_LoadTotoal[t] == Solver.quicksum(P_Load[load, t] for load in Circuit.ClustersLoads[ClusterIndex]))

    # Available Power of the cluster------------
    for t in RestorationTime:
        DSR.addConstr(P_GenTotoal[t] == Solver.quicksum(P_sub[phase, source, t] for phase in Circuit.PhaseSequence
                                                        for source in Circuit.ClustersSourceBus[ClusterIndex])
                                     +  Solver.quicksum(P_PV[phase, PV, t] for phase in Circuit.PhaseSequence
                                                        for PV in Circuit.ClustersPVs[ClusterIndex])

                                     +  Solver.quicksum(P_Line[phase, PDElement, t] for phase in Circuit.PhaseSequence
                                for PDElement in Circuit.ClustersNeighborRelations[ClusterIndex]['ParentPDElements'])

                                     -  Solver.quicksum(P_Line[phase, PDElement, t] for phase in Circuit.PhaseSequence
                                for PDElement in Circuit.ClustersNeighborRelations[ClusterIndex]['ChildrenPDElements']))

    # Available Generation and Load constraint--
    for t in RestorationTime:
        DSR.addConstr(P_LoadTotoal[t] <= P_GenTotoal[t])

    # ----------------------------------------------------------------------------------------------------------------
    # Total substation power capacity
    for t in RestorationTime:
        for Sourcebus in Circuit.ClustersSourceBus[ClusterIndex]:
            # Source node capaciy limit
            #DSR.addConstr(P_sub['a', Sourcebus, t] + P_sub['b', Sourcebus, t] + P_sub['c', Sourcebus, t] <=
            #              SourceCapacity[Sourcebus][RestorationTime.index(t)])

            DSR.addConstr(P_sub['a', Sourcebus, t] + P_sub['b', Sourcebus, t] + P_sub['c', Sourcebus, t] <=
                           SourceCapacity[Sourcebus][RestorationTime.index(t)])
            #DSR.addConstr(P_sub['a', Sourcebus, t] <= (1 / 3) * SourceCapacity[Sourcebus][RestorationTime.index(t)])
            #DSR.addConstr(P_sub['b', Sourcebus, t] <= (1 / 3) * SourceCapacity[Sourcebus][RestorationTime.index(t)])
            #DSR.addConstr(P_sub['c', Sourcebus, t] <= (1 / 3) * SourceCapacity[Sourcebus][RestorationTime.index(t)])

            # Source related PDElement capacity limit-------------
            # Maximum capacity
            DSR.addConstr(P_Line['a', Circuit.Source[Sourcebus], t] + P_Line['b', Circuit.Source[Sourcebus], t] +
                     P_Line['c', Circuit.Source[Sourcebus], t] <= SourceCapacity[Sourcebus][RestorationTime.index(t)])
            # DSR.addConstr(P_Line['a', Circuit.Source[Sourcebus], t] <=
            #                (1/3)*SourceCapacity[Sourcebus][RestorationTime.index(t)])
            # DSR.addConstr(P_Line['b', Circuit.Source[Sourcebus], t] <=
            #                (1 / 3) * SourceCapacity[Sourcebus][RestorationTime.index(t)])
            # DSR.addConstr(P_Line['c', Circuit.Source[Sourcebus], t] <=
            #                (1 / 3) * SourceCapacity[Sourcebus][RestorationTime.index(t)])

            # Minimum reverse power(here reverse power flow is not allowed)
            DSR.addConstr(P_Line['a', Circuit.Source[Sourcebus], t] >= 0)
            DSR.addConstr(P_Line['b', Circuit.Source[Sourcebus], t] >= 0)
            DSR.addConstr(P_Line['c', Circuit.Source[Sourcebus], t] >= 0)

    # ----------------------------------------------------------------------------------------------------------------
    # Voltage
    # Un-enegized buses cannot be detected to make their voltage zero!!!!!!!
    for t in RestorationTime:
        for bus in Circuit.ClusterBusVariables[ClusterIndex]:
            if bus not in Circuit.ClustersSourceBus[ClusterIndex]:
                for phase in Circuit.PhaseSequence:
                    if Circuit.BusesPhaseMatrix[bus][Circuit.PhaseSequence.index(phase)] == 1:
                        DSR.addConstr(U_Bus[phase, bus, t] <= U_max*Z[Last_iter][t][ClusterIndex]['x_Bus'][bus])
                        DSR.addConstr(U_min*Z[Last_iter][t][ClusterIndex]['x_Bus'][bus] <= U_Bus[phase, bus, t])

                    elif Circuit.BusesPhaseMatrix[bus][Circuit.PhaseSequence.index(phase)] == 0:
                        DSR.addConstr(U_Bus[phase, bus, t] == 0)

    # Substation Voltage
    for t in RestorationTime:
        for SourceBus in Circuit.ClustersSourceBus[ClusterIndex]:
            for phase in Circuit.PhaseSequence:
                DSR.addConstr(U_Bus[phase, SourceBus, t] == USource[SourceBus][t])

    # ----------------------------------------------------------------------------------------------------------------
    # Branch Power capacity (can be covered by power balance equations)
    for t in RestorationTime:
        for PDElement in Circuit.ClustersPDElements[ClusterIndex]:
            for phase in Circuit.PhaseSequence:
                if Circuit.PDElementPhaseMatrix[PDElement][Circuit.PhaseSequence.index(phase)] == 0:
                    DSR.addConstr(P_Line[phase, PDElement, t] == 0)
                    DSR.addConstr(Q_Line[phase, PDElement, t] == 0)
                else:
                    if PDElement in Circuit.ClustersSwitchableLines[ClusterIndex]:
                        DSR.addConstr(P_Line[phase, PDElement, t] <= 3*Z[Last_iter][t][ClusterIndex]['x_Switch'][PDElement])
                        DSR.addConstr(-3*Z[Last_iter][t][ClusterIndex]['x_Switch'][PDElement] <= P_Line[phase, PDElement, t])
                        DSR.addConstr(Q_Line[phase, PDElement, t] <= 3 * Z[Last_iter][t][ClusterIndex]['x_Switch'][PDElement])
                        DSR.addConstr(-3 * Z[Last_iter][t][ClusterIndex]['x_Switch'][PDElement] <= Q_Line[phase, PDElement, t])
                    else:
                        DSR.addConstr(P_Line[phase, PDElement, t] <= 3)
                        DSR.addConstr(-3 <= P_Line[phase, PDElement, t])
                        DSR.addConstr(Q_Line[phase, PDElement, t] <= 3)
                        DSR.addConstr(-3 <= Q_Line[phase, PDElement, t])
    # ----------------------------------------------------------------------------------------------------------------
    # Capacitor banks capacity
    for t in RestorationTime:
        for CB in Circuit.ClustersCBs[ClusterIndex]:
            for phase in Circuit.PhaseSequence:
                if Circuit.CapacitorConnection[CB][1][Circuit.PhaseSequence.index(phase)] != 0:
                    DSR.addConstr(Q_cap[phase, CB, t] <= Circuit.CapacitorConnection[CB][1][Circuit.PhaseSequence.index(phase)])

    # ----------------------------------------------------------------------------------------------------------------
    # Generation Constraints
    # PV Constraints*****************************************************************
    # chose one of linearized or quadratic constraints
    # Approximate Linearized ------
    # Maximum PV inverter capacity-------------------
    for t in RestorationTime:
        for PV in Circuit.ClustersPVs[ClusterIndex]:
            for phase in Circuit.PhaseSequence:
                DSR.addConstr(P_PV[phase, PV, t] <= Z[Last_iter][t][ClusterIndex]['x_Bus'][Circuit.PVconnection[PV][0]] *
                              Circuit.PVconnection[PV][1][Circuit.PhaseSequence.index(phase)])
                # Reactive power (add %20 more capacity to inverter to compensate reactive power)
                DSR.addConstr(Q_PV[phase, PV, t] <= 0.2*Circuit.PVconnection[PV][1][Circuit.PhaseSequence.index(phase)])
                DSR.addConstr(-0.2*Circuit.PVconnection[PV][1][Circuit.PhaseSequence.index(phase)] <= Q_PV[phase, PV, t])

    # quadratic constraints-----
    '''for t in RestorationTime:
        for PV in Circuit.ClustersPVs[ClusterIndex]:
            for phase in Circuit.PhaseSequence:
                DSR.addQConstr(P_PV[phase, PV, t]*P_PV[phase, PV, t] + Q_PV[phase, PV, t]*Q_PV[phase, PV, t] <=
                Z[Last_iter][t][ClusterIndex]['x_Bus'][Circuit.PVconnection[PV][0]]*
                               (Circuit.PVconnection[PV][1][Circuit.PhaseSequence.index(phase)])*
                               (Circuit.PVconnection[PV][1][Circuit.PhaseSequence.index(phase)]))'''

    # ----------------------------------------------------------------------------------------------------------------
    # Power flow constraints:
    # Power Balance--------------------------------------------------
    for t in RestorationTime:
        for PDElement in Circuit.ClustersPDElements[ClusterIndex]:
            if PDElement not in Circuit.ClustersNeighborRelations[ClusterIndex]['ChildrenPDElements']:
                if PDElement not in Circuit.ClustersSwitchableLines[ClusterIndex]:
                    # if line is not switchable
                    for phase in Circuit.PhaseSequence:
                        # Active Power-----------------------------------
                        DSR.addConstr(P_Line[phase, PDElement, t] == (Circuit.ActivePowerBalance(t, PDElement, phase,
                            P_Load, P_Line, P_PV, P_sub, Circuit.PDElementsConnections, Circuit.PDElementPhaseMatrix, Circuit.LoadToBus,
                                Circuit.LoadsConnection, Circuit.ChildrenLines, Circuit.PVBus, Circuit.ClustersSourceBus[ClusterIndex],
                            Circuit.Source)*Circuit.PDElementPhaseMatrix[PDElement][Circuit.PhaseSequence.index(phase)]))

                        # Reactive Power-----------------------------------
                        DSR.addConstr(Q_Line[phase, PDElement, t] == (Circuit.ReactivePowerBalance(t, PDElement, phase,
                                Q_Load, Q_Line, Q_PV, Q_sub, Q_cap, Circuit.PDElementsConnections, Circuit.PDElementPhaseMatrix, Circuit.LoadToBus,
                                    Circuit.LoadsConnection, Circuit.ChildrenLines, Circuit.PVBus, Circuit.ClustersSourceBus[ClusterIndex],
                                    Circuit.Source, Circuit.CapacitorBus, Circuit.CapacitorConnection)
                                            *Circuit.PDElementPhaseMatrix[PDElement][Circuit.PhaseSequence.index(phase)]))
                else:
                    # if line is switchable
                    if Z[Last_iter][t][ClusterIndex]['x_Switch'][PDElement] == 1: # if switchable line is connected
                        for phase in Circuit.PhaseSequence:
                            # Active Power-----------------------------------
                            DSR.addConstr( P_Line[phase, PDElement, t] == (Circuit.ActivePowerBalance(t, PDElement, phase,
                            P_Load, P_Line, P_PV, P_sub, Circuit.PDElementsConnections, Circuit.PDElementPhaseMatrix, Circuit.LoadToBus,
                            Circuit.LoadsConnection, Circuit.ChildrenLines, Circuit.PVBus, Circuit.ClustersSourceBus[ClusterIndex],
                            Circuit.Source)*Circuit.PDElementPhaseMatrix[PDElement][Circuit.PhaseSequence.index(phase)]))

                            # Reactive Power-----------------------------------
                            DSR.addConstr(Q_Line[phase, PDElement, t] == (Circuit.ReactivePowerBalance(t, PDElement, phase,
                                Q_Load, Q_Line, Q_PV, Q_sub, Q_cap, Circuit.PDElementsConnections, Circuit.PDElementPhaseMatrix, Circuit.LoadToBus,
                                    Circuit.LoadsConnection, Circuit.ChildrenLines, Circuit.PVBus, Circuit.ClustersSourceBus[ClusterIndex],
                                    Circuit.Source, Circuit.CapacitorBus, Circuit.CapacitorConnection)
                                            *Circuit.PDElementPhaseMatrix[PDElement][Circuit.PhaseSequence.index(phase)]))

    # Voltage equation-------------------------------------------------
    for t in RestorationTime:
        for PDElement in Circuit.ClusterCurrentVariables[ClusterIndex]:
            for phase in Circuit.PhaseSequence:
                if Circuit.PDElementPhaseMatrix[PDElement][Circuit.PhaseSequence.index(phase)] == 1:
                    Bus_Temp = Circuit.PDElementsConnections[PDElement][1]
                    ParentBus_Temp = Circuit.PDElementsConnections[PDElement][0]
                    # if line is switchable
                    if PDElement in Circuit.ClustersSwitchableLines[ClusterIndex]:
                        if Z[Last_iter][t][ClusterIndex]['x_Switch'][PDElement] == 1:  # if switchable line is connected
                            if Circuit.BusesPhaseMatrix[Bus_Temp][Circuit.PhaseSequence.index(phase)] == 1:
                                DSR.addConstr(U_Bus[phase, Bus_Temp, t] == (Circuit.VoltageRHS(t, phase, PDElement,
                                ParentBus_Temp, U_Bus, P_Line, Q_Line, Circuit.PDElementsR_hat, Circuit.PDElementsX_hat)))

                                # DSR.addConstr(U_Bus[phase, Bus_Temp, t] - (Circuit.VoltageRHS(t, phase, PDElement,
                                # ParentBus_Temp, U_Bus, P_Line, Q_Line, Circuit.PDElementsR_hat, Circuit.PDElementsX_hat))
                                #      <= (1-Z[Last_iter][t][ClusterIndex]['x_Switch'][PDElement]) * M)
                                #
                                # DSR.addConstr(U_Bus[phase, Bus_Temp, t] - (Circuit.VoltageRHS(t, phase, PDElement,
                                # ParentBus_Temp, U_Bus, P_Line, Q_Line, Circuit.PDElementsR_hat, Circuit.PDElementsX_hat))
                                #      >= -(1 - Z[Last_iter][t][ClusterIndex]['x_Switch'][PDElement]) * M)

                    # if line is not switchable
                    else:
                        if Circuit.BusesPhaseMatrix[Bus_Temp][Circuit.PhaseSequence.index(phase)] == 1:
                            DSR.addConstr(U_Bus[phase, Bus_Temp, t] == (Circuit.VoltageRHS(t, phase, PDElement,
                            ParentBus_Temp, U_Bus, P_Line, Q_Line, Circuit.PDElementsR_hat, Circuit.PDElementsX_hat)))

    # ***************************************************************************************************************
    # Connectivity and Sequencing Constraints:*********************************************************

    # Lines connectivity ====================================
    # Out of Service Lines (Faulted Lines)----------
    for t in RestorationTime:
        for fault_line in FaultedLines:
            if fault_line in Circuit.ClustersPDElements[ClusterIndex]:
                Bus_Temp = Circuit.PDElementsConnections[fault_line][1]
                ParentBus_Temp = Circuit.PDElementsConnections[fault_line][0]
                #if faulted lin connect to children cluster we do not need the voltage of children cluster node to be zero
                if fault_line in Circuit.ClustersNeighborRelations[ClusterIndex]['ChildrenPDElements']:
                    for phase in Circuit.PhaseSequence:
                        DSR.addConstr(P_Line[phase, fault_line, t] == 0)
                        DSR.addConstr(Q_Line[phase, fault_line, t] == 0)
                        #DSR.addConstr(U_Bus[phase, ParentBus_Temp, t] == 0)
                #if faulted line does not connect to the children cluster
                else:
                    for phase in Circuit.PhaseSequence:
                        DSR.addConstr(P_Line[phase, fault_line, t] == 0)
                        DSR.addConstr(Q_Line[phase, fault_line, t] == 0)
                        #DSR.addConstr(U_Bus[phase, Bus_Temp, t] == 0)
                        #DSR.addConstr(U_Bus[phase, ParentBus_Temp, t] == 0)


    # ***************************************************************************************************************
    # Optimization Results
    # *****************************************************************************************************
    # Write optimization problem in text file-----------------------------------
    #DSR.write('DSR_Function1.lp')

    # ***************************************************************************************************************
    # Perform Optimization
    # ***************************************************************************************************************
    DSR.optimize()

    # ***************************************************************************************************************
    # Optimization Results
    # *****************************************************************************************************

    X = Circuit.Update_X_Polish(DSR, X, ClusterIndex)
    # Update loads status in Z for Polish Phase
    Z = Circuit.Update_Z_Polish(DSR, Z, ClusterIndex, Last_iter)

    return [X, Z]

