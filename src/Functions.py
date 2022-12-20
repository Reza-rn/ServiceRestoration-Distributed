# ******************************************************************************************
#  *** This function returns connection of PD elements in terms of the buses
# ******************************************************************************************
def PDElementsConnectionsFind(CktIncMatrix, AllPDElements, AllBuses):
    i = 0;
    PDElementsConnections = {}
    while i <= len(CktIncMatrix) - 6:
        if CktIncMatrix[i + 2] == 1:
            To = AllBuses[CktIncMatrix[i + 1]]
            In = AllBuses[CktIncMatrix[i + 4]]
        else:
            To = AllBuses[CktIncMatrix[i + 4]]
            In = AllBuses[CktIncMatrix[i + 1]]
        PDElementsConnections[AllPDElements[CktIncMatrix[i]]] = (To, In)
        i = i + 6
    return PDElementsConnections

# ******************************************************************************************
#  *** This function is used to find removed PD Elements because of same connection buses
# ******************************************************************************************
def FindRemovedPDElemnts(PDElementsConnections):

    Removed={}
    RemovedConnection = []
    PDElementsConnectionsReduced = PDElementsConnections.copy()

    for i, j in PDElementsConnections.items():
        for k, m in PDElementsConnections.items():
            if i != k:
                if j == m:
                    if m not in RemovedConnection:
                        RemovedConnection.append(m)
                        Removed[k] = i
                        del PDElementsConnectionsReduced[k]

    AllPDElementsReduced = list(PDElementsConnectionsReduced.keys())

    ReplacedRemoved={}
    for i, j in Removed.items():
        if j in ReplacedRemoved:
            ReplacedRemoved[j].append(i)
        else:
            ReplacedRemoved[j] = i

    return [AllPDElementsReduced, PDElementsConnectionsReduced, ReplacedRemoved, Removed]

# ******************************************************************************************
# *** This function is to find parent buses of each bus using PDElementsConnections
# ******************************************************************************************
def FindParentBus(PDElementsConnections):
    ParentBusRelation={}
    for PDElelement in PDElementsConnections:
        bus = PDElementsConnections[PDElelement][1]
        ParentBusRelation[bus] = PDElementsConnections[PDElelement][0]
    return ParentBusRelation

# ******************************************************************************************
# *** This function is to find all parent buses of each bus using PDElementsConnections
# ******************************************************************************************
def FindParentBusAll(PDElementsConnections):
    ParentBusRelation={}
    for PDElelement in PDElementsConnections:
        bus = PDElementsConnections[PDElelement][1]
        if bus in ParentBusRelation:
            ParentBusRelation[bus].append(PDElementsConnections[PDElelement][0])
        else:
            ParentBusRelation[bus] = [PDElementsConnections[PDElelement][0]]
    return ParentBusRelation


# ******************************************************************************************
# *** This function is to find children buses of each bus using Parent matrix
# ******************************************************************************************
def FindChildrenBuses(PDElementsConnections):
    Children = {}
    for PDElelement in PDElementsConnections:
        bus = PDElementsConnections[PDElelement][0]
        if bus not in Children:
            Children[bus] = [PDElementsConnections[PDElelement][1]]
        else:
            Children[bus].append(PDElementsConnections[PDElelement][1])

    return Children

# ******************************************************************************************
# *** This function is used to find Neighbor PD Elements of each bus
# ******************************************************************************************
def FindBusNeighborPDElements(AllBuses, PDElementsConnections):
    Neighbor={}
    for i in AllBuses:
        for j, k in PDElementsConnections.items():
            if (k[0] == i) or (k[1] == i):
                if i in Neighbor:
                    Neighbor[i].append(j)
                else:
                    Neighbor[i] = [j]

    return Neighbor

# ******************************************************************************************
# *** This function is used to find Neighbor Buses of each Bus
# ******************************************************************************************
def FindBusNeighborBuses(AllBuses, BusNeighborPDElements, PDElementsConnections):
    Neighbor={}
    for i in AllBuses:
        for j in BusNeighborPDElements[i]:
            if PDElementsConnections[j][0] == i:
                if i in Neighbor:
                    if PDElementsConnections[j][1] not in Neighbor[i]:
                        Neighbor[i].append(PDElementsConnections[j][1])
                else:
                    Neighbor[i] = [PDElementsConnections[j][1]]

            elif PDElementsConnections[j][1] == i:
                if i in Neighbor:
                    if PDElementsConnections[j][0] not in Neighbor[i]:
                        Neighbor[i].append(PDElementsConnections[j][0])
                else:
                    Neighbor[i] = [PDElementsConnections[j][0]]

    return Neighbor
# ******************************************************************************************
# *** This function is organized Neighbor PD Elemtns for considering a constraint
# ******************************************************************************************
def OrganizeNeighborPDElements(time, Bus, x_Line, x_MT, x_PV, BusNeighborPDElements, MTBus, PVBus):
    Neighbor = 0
    for i in BusNeighborPDElements[Bus]:
        Neighbor = Neighbor + x_Line[i, time]

    # If MT connects to this Bus
    if MTBus[Bus] != []:
        Neighbor = Neighbor + x_MT[MTBus[Bus][0], time]

    # If PV connects to this Bus
    if PVBus[Bus] != []:
        Neighbor = Neighbor + x_PV[PVBus[Bus][0], time]

    return Neighbor

# ******************************************************************************************
# *** This function is organized Neighbor PD Elemtns for considering a constraint
# ******************************************************************************************
def OrganizeNeighborPDElementsWithoutDERs(time, Bus, x_Line, BusNeighborPDElements):
    Neighbor = 0
    for i in BusNeighborPDElements[Bus]:
        Neighbor = Neighbor + x_Line[i, time]

    return Neighbor

# ***************************************************************************************************
# *** This function is to find children lines of each lines using PD Elements connection matrix
# ***************************************************************************************************
def FindChildrenLines(PDElementsConnections):
    Children={}
    for i, j in PDElementsConnections.items():
        TempChildren=[]
        for k, m in PDElementsConnections.items():
            if i is not k:
                if (j[1] is m[0]) or (j[1] is m[1]):
                    TempChildren.append(k)
        Children[i] = TempChildren

    return Children

# ***************************************************************************************************
# *** This function returns Phase Connection matrix of all PD Elements
# ***************************************************************************************************

def FindPDElementPhaseMatrix(Ckt,AllPDElements):
    PDElementPhaseConnection = {}
    for i in AllPDElements:
        Element = i.split('.')[0]
        Name = i.split('.')[1]
        PDElementPhaseConnection[i] = [0, 0, 0]
        Ckt.dssCircuit.SetActiveElement('{}.{}'.format(Element, Name))
        NumPhases = Ckt.dssCircuit.ActiveCktElement.NumPhases
        if NumPhases != 0:
            if Element == 'Line':
                for j in range(NumPhases):
                    PDElementPhaseConnection[i][Ckt.dssCircuit.ActiveCktElement.NodeOrder[j] - 1] = 1

            if Element == 'Transformer':
                for j in range(NumPhases):
                    PDElementPhaseConnection[i][Ckt.dssCircuit.ActiveCktElement.NodeOrder[j] - 1] = 1

            if Element == 'Reactor':
                for j in range(NumPhases):
                    PDElementPhaseConnection[i][Ckt.dssCircuit.ActiveCktElement.NodeOrder[j] - 1] = 1

    return PDElementPhaseConnection

# ***************************************************************************************************
# *** This function returns Phase Connection matrix of all buses
# ***************************************************************************************************
def FindBusesPhaseMatrix(Ckt,AllBuses):
    BusesPhaseMatrix={}
    for i in AllBuses:
        BusesPhaseMatrix[i]=[0,0,0]
        Ckt.dssCircuit.SetActiveBus(i)
        NumNodes = Ckt.dssCircuit.ActiveBus.NumNodes
        if NumNodes > 3:
            NumNodes=3
        if NumNodes != 0:
            for j in range(NumNodes):
                BusesPhaseMatrix[i][Ckt.dssCircuit.ActiveBus.Nodes[j] - 1] = 1


    return BusesPhaseMatrix

# ***************************************************************************************************
# *** This function finds connection matrix of each load connection
# ***************************************************************************************************
def FindLoadsConnections(Ckt, S_Base):
    LoadToBus={}
    stop = Ckt.dssCircuit.Loads.First
    while stop:
        name = Ckt.dssCircuit.Loads.Name
        Ckt.dssCircuit.SetActiveElement('Load.{}'.format(name))
        P = (Ckt.dssCircuit.Loads.kW)/S_Base
        Q = (Ckt.dssCircuit.Loads.kvar)/S_Base
        if Ckt.dssCircuit.Loads.IsDelta:
            Phases = [0, 0, 0]
            NumPhases = Ckt.dssCircuit.ActiveCktElement.NumPhases
            BusConnection = Ckt.dssCircuit.ActiveCktElement.BusNames[0].split('.')
            BusName = BusConnection[0]

            if NumPhases == 3:
                Phases = [1, 1, 1]
            elif NumPhases == 2:
                Phases[int(BusConnection[1]) - 1] = 1
                Phases[int(BusConnection[2]) - 1] = 1
            elif NumPhases == 1:
                Phases[int(BusConnection[1]) - 1] = 1


            LoadToBus[name] = [BusName, NumPhases, Phases, P, Q]

        else:
            Phases = [0, 0, 0]
            NumPhases = Ckt.dssCircuit.ActiveCktElement.NumPhases
            BusConnection = Ckt.dssCircuit.ActiveCktElement.BusNames[0].split('.')
            BusName = BusConnection[0]
            if NumPhases == 3:
                Phases = [1, 1, 1]
            else:
                for i in range(NumPhases):
                    Phases[int(BusConnection[i+1])-1] = 1
            LoadToBus[name] = [BusName, NumPhases, Phases, P, Q]

        stop = Ckt.dssCircuit.Loads.Next

    return LoadToBus
# ***************************************************************************************************
# *** This function finds connection matrix of each load to the bus
# ***************************************************************************************************
def FindLoadToBus(AllBuses, LoadConnections):
    LoadToBus = {}
    for i in AllBuses:
        for k, m in LoadConnections.items():
            if i == m[0]:
                if i in LoadToBus.keys():
                    LoadToBus[i].append(k)
                else:
                    LoadToBus[i] = [k]

    return LoadToBus

# ***************************************************************************************************
# *** Define Loads weight
# ***************************************************************************************************
def FindLoadsWeight(AllLoads, **kwargs):
    LoadsWieght = {}
    for i in AllLoads:
        LoadsWieght[i] = 1

    if len(kwargs):
        for i, v in kwargs.items():
            Load = '{}'.format(i)
            if Load in LoadsWieght:
                LoadsWieght[Load] = v

    return LoadsWieght

# ***************************************************************************************************
# *** This function converts source capacity into p.u.
# ***************************************************************************************************
def FindPerUnitofSourceCapacity(SourceCapacityOld, S_Base):
    SourceCapacity = {}

    for source, capacityArray in SourceCapacityOld.items():
        for capacity in capacityArray:
            if source not in SourceCapacity:
                SourceCapacity[source] = [capacity/S_Base]
            else:
                SourceCapacity[source].append(capacity/S_Base)

    return SourceCapacity

# ***************************************************************************************************
# *** This function finds maximum capacity of the loads
# ***************************************************************************************************
def FindLoadsCapacity(Ckt,S_Base):
    LoadsCapacity={}
    Ckt.dssCircuit.Loads.First
    LoadName=Ckt.dssCircuit.Loads.Name
    Power = [(Ckt.dssCircuit.Loads.kW)/S_Base,(Ckt.dssCircuit.Loads.kvar)/S_Base]
    LoadsCapacity[LoadName] = Power
    while 1:
        stop = Ckt.dssCircuit.Loads.Next
        if stop == 0:
            break
        LoadName = Ckt.dssCircuit.Loads.Name
        Power = [(Ckt.dssCircuit.Loads.kW)/S_Base, (Ckt.dssCircuit.Loads.kvar)/S_Base]
        LoadsCapacity[LoadName] = Power
    return LoadsCapacity

# ***************************************************************************************************
# *** Define a matrix for power flow equations
# ***************************************************************************************************
def Define_aMatrix(sqrt,np):
    a1 = complex(-0.5 , -sqrt(3)/2)
    a2 = complex(-0.5 , sqrt(3)/2)
    a = np.array([[1, a1, a2],[a2, 1, a1], [a1, a2, 1]])
    return a
# ***************************************************************************************************
# *** Binary Projection function
# ***************************************************************************************************
def Binary_Projection(x):

    if abs(x-1) <= abs(x):
        return 1
    else:
        return 0

# ***************************************************************************************************
# *** Find dispatchable loads based on random selection of AllLoads
# ***************************************************************************************************
def FindDisptachableLoads(AllLoads, Start, End):
    Loads = []
    for load in range(Start, End):
        Loads.append(AllLoads[load])

    return Loads

# ***************************************************************************************************
# *** Define load type matrix
# ***************************************************************************************************
def DefineLoadType(AllLoads, DispatchableLoads):
    LoadType = {}
    for i in AllLoads:
        LoadType[i] = 0

    if len(DispatchableLoads):
        for j in DispatchableLoads:
            LoadType[j] = 1

    return LoadType

# ***************************************************************************************************
# *** Define MT matrix (Not using now!!)
# ***************************************************************************************************
def DefineMTConnection(AllBuses, S_Base, **kwargs):
    MTBus = {}
    for i in AllBuses:
        MTBus[i] = [0, 0]

    if len(kwargs):
        for i, j in kwargs.items():
            i = i[3::]
            P = j[0]/S_Base
            Q = j[1]/S_Base
            BS = j[2]
            MTBus[i] = [P , Q, BS]

    return MTBus
# ***************************************************************************************************
# *** Define MT to Bus matrix
# ***************************************************************************************************
def FindMTconnection(S_Base, **kwargs):
    MTBus = {}

    if len(kwargs):
        for i, j in kwargs.items():
            i = i[3::]
            MTname = 'MTBus{}'.format(i)
            P = j[0] / S_Base
            Q = j[1] / S_Base
            Ramp = j[2] / S_Base
            BS = j[3] if j[3] in j else 1
            MTBus[MTname] = [i, [P, Q], Ramp, BS]

    return MTBus
# ***************************************************************************************************
# *** Find which Bus has MT connection
# ***************************************************************************************************
def FindMTtoBus(AllBuses, MTconnection):
    MTtoBus = {}
    for i in AllBuses:
        MTtoBus[i] = []
    for i, j in MTconnection.items():
        MTtoBus[j[0]].append(i)

    return MTtoBus

# ***************************************************************************************************
# *** Define PV matrix (Not using now!!)
# ***************************************************************************************************
def DefinePVConnection(AllBuses, S_Base, **kwargs):
    PVBus = {}
    for i in AllBuses:
        PVBus[i] = [0, 0]

    if len(kwargs):
        for i, j in kwargs.items():
            i=i[3::]
            P = j[0]/S_Base
            BS = j[1]
            PVBus[i] = [P, BS]

    return PVBus

# ***************************************************************************************************
# *** Define PV to Bus matrix
# ***************************************************************************************************
def FindPVConnection(Ckt123, S_Base, **kwargs):
    PVBus = {}

    if len(kwargs):
        for i, j in kwargs.items():
            i = i[3::]
            PVname = 'PVBus{}'.format(i)
            P = j[0] / S_Base
            BS = j[1] if j[1] in j else 1
            PVBus[PVname] = [i, P, BS]

    # Find PV Systems in OpenDSS circuit
    stop = Ckt123.dssCircuit.PVSystems.First
    while stop:
        Name = Ckt123.dssCircuit.PVSystems.Name
        Ckt123.dssCircuit.SetActiveElement(Name)
        Bus = Ckt123.dssCircuit.ActiveCktElement.BusNames[0].split('.')[0]
        NumPhases = Ckt123.dssCircuit.ActiveCktElement.NumPhases
        if NumPhases == 1:
            ConnectedPhase = Ckt123.dssCircuit.ActiveCktElement.NodeOrder[0]  # Connected Phase in 1-Phase PV System

        if Ckt123.dssCircuit.PVSystems.kVArated:
            Total_PV_Capacity = Ckt123.dssCircuit.PVSystems.kVArated / S_Base
        elif Ckt123.dssCircuit.PVSystems.kW:
            Total_PV_Capacity = Ckt123.dssCircuit.PVSystems.kW / S_Base

        # Build PV connection Matrix
        if NumPhases == 3:
            EachPhaseCapacity = Total_PV_Capacity / NumPhases
            PVBus[Name] = [Bus, [EachPhaseCapacity, EachPhaseCapacity, EachPhaseCapacity], 1]
        elif NumPhases == 1:
            PVBus[Name] = [Bus, [0, 0, 0], 1]
            PVBus[Name][1][ConnectedPhase - 1] = Total_PV_Capacity
        stop = Ckt123.dssCircuit.PVSystems.Next
    return PVBus


# ***************************************************************************************************
# *** Find which Bus has PV connection
# ***************************************************************************************************
def FindPVtoBus(AllBuses, PVconnection):
    PVtoBus = {}
    for i in AllBuses:
        PVtoBus[i] = []
    for i, j in PVconnection.items():
        PVtoBus[j[0]].append(i)

    return PVtoBus

# ***************************************************************************************************
# *** Define switchable PD Elements and lines
# ***************************************************************************************************
def FindSwitchableLines(AllPDElements, *args):
    SwitchableLines = {}
    for i in AllPDElements:
        Name = i.split('.')[1]
        Name = Name[0:2].lower()
        if Name == 'sw':
            SwitchableLines[i] = 1
        else:
            SwitchableLines[i] = 0
    if len(args):
        for i in args:
            SwitchableLines[i] = 1

    return SwitchableLines

# ***************************************************************************************************
# *** This function is used to list all switches
# ***************************************************************************************************
def FindSwitches(SwitchableLines):
    Switches = []
    for i, j in SwitchableLines.items():
        if j == 1:
            Switches.append(i)

    return Switches
# ***************************************************************************************************
# *** Find all lines length
# ***************************************************************************************************
def FindLinesLength(Ckt):
    LinesLength = {}
    stop = Ckt.dssCircuit.Lines.First
    while stop:
        Length = Ckt.dssCircuit.Lines.Length
        Name = Ckt.dssCircuit.Lines.Name
        LinesLength[Name] = Length
        stop = Ckt.dssCircuit.Lines.Next

    return LinesLength

# ***************************************************************************************************
# *** Integrate All PD Elements into lines length
# ***************************************************************************************************
def FindPDElementLength(AllPDElements, LinesLength):
    PDElementsLength = {}
    for i in AllPDElements:
        Name = i.split('.')[1]
        if Name in LinesLength:
            PDElementsLength[i] = LinesLength[Name]
        else:
            PDElementsLength[i] = 10 #ft

    return PDElementsLength

# ***************************************************************************************************
# *** Find all lines square capacity
# ***************************************************************************************************
def FindLinesSquareCapacity(Ckt, I_Base):
    LinesSquareCapacity = {}
    stop = Ckt.dssCircuit.Lines.First
    while stop:
        Name = Ckt.dssCircuit.Lines.Name
        Capacity = Ckt.dssCircuit.Lines.NormAmps
        Capacity_pu = Capacity/I_Base
        SquareCapacity = Capacity_pu * Capacity_pu
        LinesSquareCapacity[Name]= SquareCapacity
        stop = Ckt.dssCircuit.Lines.Next

    return LinesSquareCapacity

# ***************************************************************************************************
# *** Integrate all PD Elements into Lines square capacity
# ***************************************************************************************************
def FindPDElementsSquareCapacity(AllPDElements, LinesSquareCapacity):
    PDElementsSquareCapacity = {}
    for i in AllPDElements:
        Name = i.split('.')[1]
        if Name in LinesSquareCapacity:
            PDElementsSquareCapacity[i] = LinesSquareCapacity[Name]
        else:
            PDElementsSquareCapacity[i] = 1  # pu

    return PDElementsSquareCapacity

# ***************************************************************************************************
# *** Find Capacitor banks connection
# ***************************************************************************************************
def FindCapacitorConnection(Ckt, S_Base):
    CapacitorConnection = {}

    Ckt.dssCircuit.Capacitors.First
    Name = Ckt.dssCircuit.Capacitors.Name
    NumPhases = Ckt.dssCircuit.ActiveElement.NumPhases
    TotalCapacity = Ckt.dssCircuit.Capacitors.kvar
    EachPhaseCapacity = (TotalCapacity/NumPhases) / S_Base

    Ckt.dssCircuit.SetActiveElement(Name)
    BusOrder = Ckt.dssCircuit.ActiveElement.BusNames[0].split('.')
    BusName = BusOrder[0]
    CapacitorConnection[Name] = [BusName, [0, 0, 0]]
    if NumPhases == 3:
        CapacitorConnection[Name][1] = [EachPhaseCapacity, EachPhaseCapacity, EachPhaseCapacity]
    else:
        for j in range(NumPhases):
            Phase = int(BusOrder[j + 1]) - 1
            CapacitorConnection[Name][1][Phase] = EachPhaseCapacity

    NumCapacitors = Ckt.dssCircuit.Capacitors.Count
    for i in range(NumCapacitors-1):
        Ckt.dssCircuit.Capacitors.Next
        Name = Ckt.dssCircuit.Capacitors.Name
        NumPhases = Ckt.dssCircuit.ActiveElement.NumPhases
        TotalCapacity = Ckt.dssCircuit.Capacitors.kvar
        EachPhaseCapacity = (TotalCapacity / NumPhases)/ S_Base

        Ckt.dssCircuit.SetActiveElement(Name)
        BusOrder = Ckt.dssCircuit.ActiveElement.BusNames[0].split('.')
        BusName = BusOrder[0]
        CapacitorConnection[Name] = [BusName, [0, 0, 0]]
        if NumPhases == 3:
            CapacitorConnection[Name][1] = [EachPhaseCapacity, EachPhaseCapacity, EachPhaseCapacity]
        else:
            for j in range(NumPhases):
                Phase = int(BusOrder[j + 1])-1
                CapacitorConnection[Name][1][Phase] = EachPhaseCapacity

    return CapacitorConnection

# ***************************************************************************************************
# *** Find which Bus has capacitor bank
# ***************************************************************************************************
def FindCapacitorToBus(AllBuses, CapacitorConnection):
    CapacitortoBus = {}
    for i in AllBuses:
        CapacitortoBus[i] = []
    for i, j in CapacitorConnection.items():
        CapacitortoBus[j[0]].append(i)

    return CapacitortoBus

# ***************************************************************************************************
# *** Find the PD element which is connected to the source bus
# ***************************************************************************************************
def FindFirstPDElement(PDElementsConnections, SourceBus):
    PDElement = None
    for i, j in PDElementsConnections.items():
        if SourceBus[0].lower() == j[0].lower():
            PDElement = i

    return PDElement

# ***************************************************************************************************
# *** Set restoration times
# ***************************************************************************************************
def SetRestorationTimes(numbers):
    Time = []
    for i in range(numbers):
        Time.append('t{}'.format(i+1))
    return Time

# ***************************************************************************************************
# *** Find Source bus Connections
# ***************************************************************************************************
def FindSourceconnection(PDElementsConnection, SourceBus):
    Source = {}
    if SourceBus != []:
        for i in SourceBus:
            Source[i] = 0

    for k in SourceBus:
        for i, j in PDElementsConnection.items():
            if (j[0] == k) or (j[1] == k):
                Source[k] = i

    return Source

# *******************************Build_VariableSave********************************************************************
# *** Convert elements of a vector to base values
# ***************************************************************************************************
def Convert_pu(x, base):
    output = []
    for i in x:
        output.append(i/base)

    return output

# ***************************************************************************************************
# *** Convert elements of a vector to base values
# ***************************************************************************************************
def ArrangeSourceCapacity(SourceBus, S_Base, SourceCapacity):
    output = {}
    for i in SourceBus:
        output[i] = SourceCapacity[SourceBus.index(i)]/S_Base

    return output

# ***************************************************************************************************
# *** Find Resistance of distribution lines
# ***************************************************************************************************
def FindLinesResistance(Ckt, AllLines):
    LinesResistance={}
    for i in AllLines:
        Ckt.dssCircuit.SetActiveElement('Line.{}'.format(i))
        R = list(Ckt.dssCircuit.Lines.Rmatrix)
        LinesResistance[i] = R

    return LinesResistance

# ***************************************************************************************************
# *** Find Reactance of distribution lines
# ***************************************************************************************************
def FindLinesReactance(Ckt, AllLines):
    LinesReactance={}
    for i in AllLines:
        Ckt.dssCircuit.SetActiveElement('Line.{}'.format(i))
        X = list(Ckt.dssCircuit.Lines.Xmatrix)
        LinesReactance[i] = X

    return LinesReactance


# ***************************************************************************************************
# *** Find Resistance of All PD Elements
# ***************************************************************************************************
def FindPDElementsR(AllPDElements, PDElementPhaseMatrix, LinesResistance, PDElementLength, Z_Base):
    PDElementResistance={}
    for i in AllPDElements:
        Type = i.split('.')[0]
        Name = i.split('.')[1]
        if Type == 'Line':
            Rmatrix = LinesResistance[Name]
            # Multiply by length and Convert to p.u.------------
            if Name[0:1] == 'sw':
                for j in Rmatrix:
                  Rmatrix[Rmatrix.index(j)] = j / Z_Base
            else:
                for j in Rmatrix:
                  Rmatrix[Rmatrix.index(j)] = j * PDElementLength[i] / Z_Base
            # --------------------------------------------------
            FinalRmatrix = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
            Connection = PDElementPhaseMatrix[i]

            if len(Rmatrix) == 9: # Three Phase-------------------------------------------
                FinalRmatrix[0][:] = Rmatrix[0:3]
                FinalRmatrix[1][:] = Rmatrix[3:6]
                FinalRmatrix[2][:] = Rmatrix[6:9]
                PDElementResistance[i] = FinalRmatrix

            elif len(Rmatrix) == 4: # Two Phase-------------------------------------------
                if (Connection[0] == 1) and (Connection[1] == 1):
                    FinalRmatrix[0][0] = Rmatrix[0]
                    FinalRmatrix[0][1] = Rmatrix[1]
                    FinalRmatrix[1][0] = Rmatrix[2]
                    FinalRmatrix[1][1] = Rmatrix[3]

                elif (Connection[0] == 1) and (Connection[2] == 1):
                    FinalRmatrix[0][0] = Rmatrix[0]
                    FinalRmatrix[0][2] = Rmatrix[1]
                    FinalRmatrix[2][0] = Rmatrix[2]
                    FinalRmatrix[2][2] = Rmatrix[3]

                elif (Connection[1] == 1) and (Connection[2] == 1):
                    FinalRmatrix[1][1] = Rmatrix[0]
                    FinalRmatrix[1][2] = Rmatrix[1]
                    FinalRmatrix[2][1] = Rmatrix[2]
                    FinalRmatrix[2][2] = Rmatrix[3]

                PDElementResistance[i] = FinalRmatrix

            elif len(Rmatrix) == 1: # One Phase-------------------------------------------
                if Connection[0] == 1:
                    FinalRmatrix[0][0] = Rmatrix[0]

                elif Connection[1] == 1:
                    FinalRmatrix[1][1] = Rmatrix[0]

                elif Connection[2] == 1:
                    FinalRmatrix[2][2] = Rmatrix[0]

                PDElementResistance[i] = FinalRmatrix

        elif Type == 'Transformer':
            FinalRmatrix = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
            Connection = PDElementPhaseMatrix[i]

            if Connection[0] == 1: # phase a
                FinalRmatrix[0][0] = 0.01 / Z_Base

            if Connection[1] == 1:  # phase b
                FinalRmatrix[1][1] = 0.01 / Z_Base

            if Connection[2] == 1:  # phase c
                FinalRmatrix[2][2] = 0.01 / Z_Base

            PDElementResistance[i] = FinalRmatrix

        elif Type == 'Reactor':
            FinalRmatrix = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
            Connection = PDElementPhaseMatrix[i]

            if Connection[0] == 1: # phase a
                FinalRmatrix[0][0] = 0.001 / Z_Base

            if Connection[1] == 1:  # phase b
                FinalRmatrix[1][1] = 0.001 / Z_Base

            if Connection[2] == 1:  # phase c
                FinalRmatrix[2][2] = 0.001 / Z_Base

            PDElementResistance[i] = FinalRmatrix

    return PDElementResistance

# ***************************************************************************************************
# *** Find Reactance of All PD Elements
# ***************************************************************************************************
def FindPDElementsX(AllPDElements, PDElementPhaseMatrix, LinesReactance, PDElementsLength, Z_Base):
    PDElementReactance = {}
    for i in AllPDElements:
        Type = i.split('.')[0]
        Name = i.split('.')[1]
        if Type == 'Line':
            Xmatrix = LinesReactance[Name]
            # Multiply by length and Convert to p.u.------------
            for j in Xmatrix:
                Xmatrix[Xmatrix.index(j)] = j * PDElementsLength[i] / Z_Base
            # --------------------------------------------------
            FinalXmatrix = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
            Connection = PDElementPhaseMatrix[i]

            if len(Xmatrix) == 9:  # Three Phase-------------------------------------------
                FinalXmatrix[0][:] = Xmatrix[0:3]
                FinalXmatrix[1][:] = Xmatrix[3:6]
                FinalXmatrix[2][:] = Xmatrix[6:9]
                PDElementReactance[i] = FinalXmatrix

            elif len(Xmatrix) == 4:  # Two Phase-------------------------------------------
                if (Connection[0] == 1) and (Connection[1] == 1):
                    FinalXmatrix[0][0] = Xmatrix[0]
                    FinalXmatrix[0][1] = Xmatrix[1]
                    FinalXmatrix[1][0] = Xmatrix[2]
                    FinalXmatrix[1][1] = Xmatrix[3]

                elif (Connection[0] == 1) and (Connection[2] == 1):
                    FinalXmatrix[0][0] = Xmatrix[0]
                    FinalXmatrix[0][2] = Xmatrix[1]
                    FinalXmatrix[2][0] = Xmatrix[2]
                    FinalXmatrix[2][2] = Xmatrix[3]

                elif (Connection[1] == 1) and (Connection[2] == 1):
                    FinalXmatrix[1][1] = Xmatrix[0]
                    FinalXmatrix[1][2] = Xmatrix[1]
                    FinalXmatrix[2][1] = Xmatrix[2]
                    FinalXmatrix[2][2] = Xmatrix[3]

                PDElementReactance[i] = FinalXmatrix

            elif len(Xmatrix) == 1:  # One Phase-------------------------------------------
                if Connection[0] == 1:
                    FinalXmatrix[0][0] = Xmatrix[0]

                elif Connection[1] == 1:
                    FinalXmatrix[1][1] = Xmatrix[0]

                elif Connection[2] == 1:
                    FinalXmatrix[2][2] = Xmatrix[0]

                PDElementReactance[i] = FinalXmatrix

        elif Type == 'Transformer':
            FinalXmatrix = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
            Connection = PDElementPhaseMatrix[i]

            if Connection[0] == 1:  # phase a
                FinalXmatrix[0][0] = 0.01 / Z_Base

            if Connection[1] == 1:  # phase b
                FinalXmatrix[1][1] = 0.01 / Z_Base

            if Connection[2] == 1:  # phase c
                FinalXmatrix[2][2] = 0.01 / Z_Base

            PDElementReactance[i] = FinalXmatrix

        elif Type == 'Reactor':
            FinalXmatrix = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
            Connection = PDElementPhaseMatrix[i]

            if Connection[0] == 1:  # phase a
                FinalXmatrix[0][0] = 0.01 / Z_Base

            if Connection[1] == 1:  # phase b
                FinalXmatrix[1][1] = 0.01 / Z_Base

            if Connection[2] == 1:  # phase c
                FinalXmatrix[2][2] = 0.01 / Z_Base

            PDElementReactance[i] = FinalXmatrix

    return PDElementReactance

# ***************************************************************************************************
# *** Calculate R_hat PD Elements
# ***************************************************************************************************
def CalculatePDElements_R_hat(PDElementsR, PDElementsX, sqrt, np):
    PDElementsR_hat = {}
    # alpah matrix------------------------------------------------
    a1 = complex(-0.5, -sqrt(3)/2)
    a2 = complex(-0.5, sqrt(3)/2)
    a = np.array([[1, a1, a2], [a2, 1, a1], [a1, a2, 1]])
    # ------------------------------------------------------------
    for i, j in PDElementsR.items():
        Matrix1 = a * np.array(j)
        Matrix2 = a * np.array(PDElementsX[i])
        R_hat = Matrix1.real - Matrix2.imag
        PDElementsR_hat[i] = R_hat.tolist()

    return PDElementsR_hat

# ***************************************************************************************************
# *** Calculate X_hat PD Elements
# ***************************************************************************************************
def CalculatePDElements_X_hat(PDElementsX, PDElementsR, sqrt, np):
    PDElementsX_hat = {}
    # alpah matrix------------------------------------------------
    a1 = complex(-0.5, -sqrt(3)/2)
    a2 = complex(-0.5, sqrt(3)/2)
    a = np.array([[1, a1, a2], [a2, 1, a1], [a1, a2, 1]])
    # ------------------------------------------------------------
    for i, j in PDElementsX.items():
        Matrix1 = a * np.array(j)
        Matrix2 = a * np.array(PDElementsR[i])
        X_hat = Matrix1.real + Matrix2.imag
        PDElementsX_hat[i] = X_hat.tolist()

    return PDElementsX_hat

# ***************************************************************************************************
# *** This function is used to calculate Z_hat
# ***************************************************************************************************
def CalculatePDElements_Z_hat(PDElementsR, PDElementsX, sqrt):
    PDElementsZ_hat = {}
    for i in PDElementsR:
        PDElementsZ_hat[i] = [[0, 0, 0],[0, 0, 0],[0, 0, 0]]

    for i, j in PDElementsR.items():
        for k in range(3):
            for m in range(3):
                PDElementsZ_hat[i][k][m] = sqrt(j[k][m]*j[k][m] + (PDElementsX[i][k][m])*(PDElementsX[i][k][m]))

    return PDElementsZ_hat

# ***************************************************************************************************
# *** This function is used in Active power balance euqation
# ***************************************************************************************************
def ActivePowerBalance(t, PDElement, Phase, P_Load, P_Line, P_PV, P_sub, PDElementsConnections,
                       PDElementPhaseMatrix, LoadToBus, LoadsConnection, ChildrenLines, PVBus, SourceBus, Source):
    PowerBalance = 0
    BusName = PDElementsConnections[PDElement][1]
    #Phase index
    if Phase == 'a':
        PhaseIndex = 0
    elif Phase == 'b':
        PhaseIndex = 1
    elif Phase == 'c':
        PhaseIndex = 2
    # Loads-----------------------------------------------------------------
    if BusName in LoadToBus:
        for load in LoadToBus[BusName]:
            if LoadsConnection[load][2][PhaseIndex] == 1:
                PowerBalance = PowerBalance + (P_Load[load, t])/LoadsConnection[load][1]
    # Lines-----------------------------------------------------------------
    if ChildrenLines[PDElement] != []:
        for childPDElement in ChildrenLines[PDElement]:
            if PDElementPhaseMatrix[childPDElement][PhaseIndex] == 1:
                if PDElementsConnections[childPDElement][0] == BusName:
                    PowerBalance = PowerBalance + P_Line[Phase, childPDElement, t]

                else:
                    PowerBalance = PowerBalance - P_Line[Phase, childPDElement, t]

    # Substations--------------------------------------------------------------
    if BusName in SourceBus:
        PowerBalance = PowerBalance - P_sub[Phase, BusName, t]

    if SourceBus != []:
        if PDElement == Source[SourceBus[0]]:
            PowerBalance = PowerBalance - P_sub[Phase, SourceBus[0], t]

    # Loss-----------------------------------------------------------------
    # Loss is ignored for simplification of problem formulation
    # Generation-----------------------------------------------------------
    # PV---------
    if PVBus[BusName] != []:
        for i in PVBus[BusName]:
            PowerBalance = PowerBalance - P_PV[Phase, i, t]

    return PowerBalance

# ***************************************************************************************************
# *** This function is used in Reactive power balance euqation
# ***************************************************************************************************
def ReactivePowerBalance(t, PDElement, Phase, Q_Load, Q_Line, Q_PV,  Q_sub, Q_cap, PDElementsConnections,
                         PDElementPhaseMatrix, LoadToBus, LoadsConnection, ChildrenLines, PVBus, SourceBus, Source,
                         CapacitorBus, CapacitorConnection):
    PowerBalance = 0
    BusName = PDElementsConnections[PDElement][1]
    # Phase index
    if Phase == 'a':
        PhaseIndex = 0
    elif Phase == 'b':
        PhaseIndex = 1
    #elif Phase == 'c':
    else:
        PhaseIndex = 2
    # Loads-----------------------------------------------------------------
    if BusName in LoadToBus:
        for load in LoadToBus[BusName]:
            if LoadsConnection[load][2][PhaseIndex] == 1:
                PowerBalance = PowerBalance + Q_Load[load, t]/LoadsConnection[load][1]

    # Lines------------------------------------------------------------------
    if ChildrenLines[PDElement] != []:
        for childPDElement in ChildrenLines[PDElement]:
            if PDElementPhaseMatrix[childPDElement][PhaseIndex] == 1:
                if PDElementsConnections[childPDElement][0] == BusName:
                    PowerBalance = PowerBalance + Q_Line[Phase, childPDElement, t]

                else:
                    PowerBalance = PowerBalance - Q_Line[Phase, childPDElement, t]

    # Substations--------------------------------------------------------------
    if BusName in SourceBus:
        PowerBalance = PowerBalance - Q_sub[Phase, BusName, t]

    if SourceBus != []:
        if PDElement == Source[SourceBus[0]]:
            PowerBalance = PowerBalance - Q_sub[Phase, SourceBus[0], t]

    # Loss-----------------------------------------------------------------
    # Loss is ignored for simplification of problem formulation
    # Generation-----------------------------------------------------------
    # PV---------
    if PVBus[BusName] != []:
        for i in PVBus[BusName]:
            PowerBalance = PowerBalance - Q_PV[Phase, i, t]

    # Capacitor Bank-----------------------------------------------------------
    if CapacitorBus[BusName] is not []:
        for CB in CapacitorBus[BusName]:
            if CapacitorConnection[CB][1][PhaseIndex] != 0:
                PowerBalance = PowerBalance - Q_cap[Phase, CB, t]

    return PowerBalance

# *****************************************************************************************
# * This function is used to find the right hand side of the voltage drop equation
# *****************************************************************************************
def VoltageRHS(t, Phase, PDElement, ParentBus, U_Bus, P_Line, Q_Line, PDElementsR_hat, PDElementsX_hat):
    # Phase index
    if Phase == 'a':
        PhaseIndex = 0
    elif Phase == 'b':
        PhaseIndex = 1
    elif Phase == 'c':
        PhaseIndex = 2

    VoltageDrop = 0

    # Parent Bus Voltage
    VoltageDrop = VoltageDrop + U_Bus[Phase, ParentBus, t]

    # r_hat*P_line section
    Temp1 = 0
    for i in range(3):
        Temp1 = Temp1 + PDElementsR_hat[PDElement][PhaseIndex][i]*P_Line[Phase, PDElement, t]

    # x_hat*Q_line section
    Temp2 = 0
    for i in range(3):
        Temp2 = Temp2 + PDElementsX_hat[PDElement][PhaseIndex][i] * Q_Line[Phase, PDElement, t]

    # Integrate both r_hat*P_line and x_hat*Q_line in original equation
    VoltageDrop = VoltageDrop - 2*(Temp1 + Temp2)

    # z_hat*l_Line section
    # remove this term because do not have access to line current for simplification

    return VoltageDrop

# *********************************************************************************
# * Calculate total loss of the distribution network by current variable
# *********************************************************************************
def CalculateTotalLossByCurrent(RestorationTime, AllPDElements, PhaseSequence, l_Line, PDElementsR):
    Loss = 0
    for t in RestorationTime:
        for j in AllPDElements:
            for i in PhaseSequence:
                # Phase index
                if i == 'a':
                    PhaseIndex = 0
                elif i == 'b':
                    PhaseIndex = 1
                elif i == 'c':
                    PhaseIndex = 2

                for k in range(3):
                    Loss = Loss + PDElementsR[j][PhaseIndex][k] * l_Line[PhaseSequence[k], j, t]

    return Loss

# *********************************************************************************
# * Calculate total loss of the distribution network by Power variable
# *********************************************************************************
def CalculateTotalLossByPower(RestorationTime, AllPDElements, PhaseSequence, P_Line, Q_Line, PDElementsR, Voltage,
                              Solver):
    Loss = Solver.QuadExpr()
    for t in RestorationTime:
        for PDElement in AllPDElements:
            for phase in PhaseSequence:
                for k in range(3):
                    Loss.add(PDElementsR[PDElement][PhaseSequence.index(phase)][k] *
                             (P_Line[phase, PDElement, t]*P_Line[phase, PDElement, t] +
                              Q_Line[phase, PDElement, t]*Q_Line[phase, PDElement, t])/Voltage[t])

    return Loss

# *******************************************************************************************************
# * This function adjusts source capacity based on provided power by the source
# *******************************************************************************************************
def CalculateUsource(SourceCapacity, RestorationTime, **kwargs):
    Usource = {}
    if len(kwargs):
        for source, voltage in kwargs.items():
            Usource[source] = {}
            for t in RestorationTime:
                if SourceCapacity[source][RestorationTime.index(t)] != 0:
                    Usource[source][t] = voltage
                else:
                    Usource[source][t] = 0

    return Usource

# *******************************************************************************************************
# * This function finds PD Element connected to each sourcebus
# *******************************************************************************************************
def FindSourcebusConnection(Source, SourceBus, PDElementsConnections, ParentBus):

    # A function to find key from value in dictionary
    key_list = list(PDElementsConnections.keys())
    val_list = list(PDElementsConnections.values())

    for i in SourceBus:
        if SourceBus.index(i) != 0:
            Parent = ParentBus[i]
            if (Parent, i) in val_list:
                PDElement = key_list[val_list.index((Parent, i))]
                Source[i] = PDElement

    return Source

# *********************************************************************************
# * Find X-Y coordinates of buses
# *********************************************************************************
def FindXY_Coordinates(ckt):
    Coordinates = {}
    Temp_buses = ckt.AllBuses()

    for i in Temp_buses:
        x = ckt.dssCircuit.Buses(i).x
        y = ckt.dssCircuit.Buses(i).y
        Coordinates[i]=[x, y]

    return Coordinates

# *******************************************************************************************************
# * This function finds neighbor lines for spanning tree constraint of restricting to have one parent
# *******************************************************************************************************
def FindParentFlow(t, Bus, Beta, PDElementsConnections, BusNeighborBuses):

    FlowDirection = 0

    # A function to find key from value in dictionary
    key_list = list(PDElementsConnections.keys())
    val_list = list(PDElementsConnections.values())

    for i in BusNeighborBuses[Bus]:
        if (Bus, i) in val_list:
            Line = key_list[val_list.index((Bus, i))]
            FlowDirection = FlowDirection + Beta[Line, 'd2', t]

        elif (i, Bus) in val_list:
            Line = key_list[val_list.index((i, Bus))]
            FlowDirection = FlowDirection + Beta[Line, 'd1', t]

    return FlowDirection

# *******************************************************************************************************
# * This function defines higher weight for normally close switches to have them close if possible
# *******************************************************************************************************
def DefineSwitchWeight(SwitchableLines, NormallyCloseSwitches, W1, W2):
    Weight = {}
    for i in SwitchableLines:
        if SwitchableLines[i] == 1:
            if i in NormallyCloseSwitches:
                Weight[i] = W1

            else:
                Weight[i] = W2

    return Weight

# *******************************************************************************************************
# * This function defines higher weight for normally close switches to have them close if possible
# *******************************************************************************************************
def ModifySwitchesCoordinates(XY_Coordinates, SwitchableLines, PDElementsConnections, ParentBus, ChildrenBuses):

    for i, j in SwitchableLines.items():
        if j == 1:
            bus1 = PDElementsConnections[i][0]
            bus2 = PDElementsConnections[i][1]
            if (XY_Coordinates[bus1][0] == 0) or (XY_Coordinates[bus1][1] == 0):
                stop = 1
                Parent = (ParentBus[bus1]+'.')[:-1]
                while stop:
                    if (XY_Coordinates[Parent][0] != 0) and (XY_Coordinates[Parent][1] != 0):
                        XY_Coordinates[bus1][0] = XY_Coordinates[Parent][0]
                        XY_Coordinates[bus1][1] = XY_Coordinates[Parent][1]
                        stop = 0
                    Parent = (ParentBus[Parent]+'.')[:-1]

            if (XY_Coordinates[bus2][0] == 0) or (XY_Coordinates[bus2][1] == 0):
                stop = 1
                Child = (ChildrenBuses[bus2][0]+'.')[:-1]
                while stop:
                    if (XY_Coordinates[Child][0] != 0) and (XY_Coordinates[Child][1] != 0):
                        XY_Coordinates[bus2][0] = XY_Coordinates[Child][0]
                        XY_Coordinates[bus2][1] = XY_Coordinates[Child][1]
                        stop = 0
                    Child = (ChildrenBuses[Child][0]+'.')[:-1]

    return XY_Coordinates

# *******************************************************************************************************
# * This function is used to find the service transfomer which the load in secondary is connected
# *******************************************************************************************************
def FindLoad2Secondarytransformer(loadname, ckt, ParentBus, PDElementsConnections):
    # Find bus name and its parent for load
    ckt.dssCircuit.SetActiveElement('Load.{}'.format(loadname))
    bus = ckt.dssCircuit.ActiveCktElement.BusNames[0].split('.')
    BusName = bus[0]
    ParentName = ParentBus[BusName]

    # Find PD Element connected to the bus
    key_list = list(PDElementsConnections.keys())
    val_list = list(PDElementsConnections.values())
    val_index = val_list.index((ParentName, BusName))
    PDElement = key_list[val_index]

    # Find parent PD Elements until reaching service transfomer
    stop = 1
    Type = PDElement.split('.')[0]
    if Type.lower() == 'Transformer'.lower():
        stop = 0

    while stop:
        BusName = (ParentName + '.')[:-1]
        ParentName = (ParentBus[BusName] +'.')[:-1]
        val_index = val_list.index((ParentName, BusName))
        PDElement = key_list[val_index]

        Type = PDElement.split('.')[0]
        if Type.lower() == 'Transformer'.lower():
            stop = 0

    return PDElement
# *******************************************************************************************************
# * This function is used modify vectors for secondary loads
# *******************************************************************************************************
def ModifySecondaryLoad(CktNum, ckt, AllLoads, ParentBus, PDElementsConnections, LoadsConnection, Loads_max, LoadToBus):

    AllLoadsModi = list(AllLoads).copy()
    LoadsConnectionModi = LoadsConnection.copy()
    Loads_maxModi = Loads_max.copy()
    LoadToBusModi = LoadToBus.copy()

    if (CktNum == 2)|(CktNum == 3):
        for i in AllLoads:
            ckt.dssCircuit.SetActiveElement('Load.{}'.format(i))
            voltage = ckt.dssCircuit.ActiveCktElement.VoltagesMagAng[0]
            if voltage <= 600:
                Transfomer = FindLoad2Secondarytransformer(i, ckt, ParentBus, PDElementsConnections)
                Bus = PDElementsConnections[Transfomer][0]
                NewLoad = Transfomer.split('.')[1] + 'load'

                ckt.dssCircuit.SetActiveElement(Transfomer) # Set Transformer as active element
                Phases = [0, 0, 0]
                NumPhases = ckt.dssCircuit.ActiveCktElement.NumPhases
                BusConnection = ckt.dssCircuit.ActiveCktElement.BusNames[0].split('.')
                Phase = BusConnection[1]

                if NumPhases==1:
                    # Remove load from AllLoadsModi
                    AllLoadsModi.remove(i)
                    # Add new load to AllLoadsModi
                    if NewLoad not in AllLoadsModi:
                        AllLoadsModi.append(NewLoad)

                    # Remove load from LoadTobus
                    LoadToBusModi.pop(LoadsConnection[i][0], None)
                    # Add new load to LoadToBusModi
                    if Bus not in LoadToBusModi:
                        LoadToBusModi[Bus] = [NewLoad]
                    else:
                        LoadToBusModi[Bus].append(NewLoad)

                    # Remove load from LoadConnectionsModi
                    LoadsConnectionModi.pop(i, None)
                    # Add or update new load to LoadsConnectionModi
                    Phases[int(Phase)-1] = 1
                    P_Load = LoadsConnection[i][3]
                    Q_Load = LoadsConnection[i][4]
                    if NewLoad not in LoadsConnectionModi:
                        LoadsConnectionModi[NewLoad] = [Bus,NumPhases, Phases, P_Load, Q_Load]
                    else:
                        LoadsConnectionModi[NewLoad][3] = LoadsConnectionModi[NewLoad][3] + P_Load
                        LoadsConnectionModi[NewLoad][4] = LoadsConnectionModi[NewLoad][4] + Q_Load

                    # Remove load from Loads_maxModi
                    Loads_maxModi.pop(i, None)
                    # Add or update new load to Loads_maxModi
                    if NewLoad not in Loads_maxModi:
                        Loads_maxModi[NewLoad] = [P_Load, Q_Load]
                    else:
                        Loads_maxModi[NewLoad][0] = Loads_maxModi[NewLoad][0] + P_Load
                        Loads_maxModi[NewLoad][1] = Loads_maxModi[NewLoad][1] + Q_Load

    return [AllLoadsModi, LoadToBusModi, LoadsConnectionModi, Loads_maxModi]

# *********************************************************************************
# * Plot circuit connections based on the energization
# *********************************************************************************
def PlotCircuit(XY_Coordinates, AllBuses, AllPDElements, PDElementsConnections, DSR, plt, PVBus, SwitchableLines,
                AllLoads, LoadsConnection, PDElementPhaseMatrix, FaultedLines, SourceBus):

    # Draw Load pickup
    ResultLoads = {}
    for i in AllLoads:
        ResultLoads[i] = 0
    for v in DSR.getVars():
        if v.varName.split('[')[0] == 'x_Load':
            ResultLoads[v.varName.split('[')[1][:-1]] = v.x

    for j, k in ResultLoads.items():
        TempBus = LoadsConnection[j][0]
        if (XY_Coordinates[TempBus][0] != 0) and (XY_Coordinates[TempBus][1] != 0):
            if k == 0:
                #plt.scatter(XY_CoordiPlot_TotalLoadnates[TempBus][0], XY_Coordinates[TempBus][1], c='r', marker='.', linewidths=0.0005)
                plt.scatter(XY_Coordinates[TempBus][0], XY_Coordinates[TempBus][1], c='r', s=0.5)#was 0.5
            elif k == 1:
                #plt.scatter(XY_Coordinates[TempBus][0], XY_Coordinates[TempBus][1], c='g', marker='.', linewidths=0.0005)
                plt.scatter(XY_Coordinates[TempBus][0], XY_Coordinates[TempBus][1], c='g', s=0.5)#was 0.5

    # Draw Buses
    for i in AllBuses:
        if (XY_Coordinates[i][0] != 0) and (XY_Coordinates[i][1] != 0):
            if PVBus[i] != []:
                plt.scatter(XY_Coordinates[i][0], XY_Coordinates[i][1], c='k', marker='D', linewidths=1)
                plt.text(XY_Coordinates[i][0] - 100, XY_Coordinates[i][1] - 100, 'PV', fontsize=5)
            if i in SourceBus:
                plt.scatter(XY_Coordinates[i][0], XY_Coordinates[i][1], c='k', marker='*', linewidths=2)
                plt.text(XY_Coordinates[i][0] - 100, XY_Coordinates[i][1] - 100, 'Sub.', fontsize=6)


    # Draw energized Lines----------------------------------------------------------
    fault = 'Fault'
    ResultLines = {}
    for i in AllPDElements:
        ResultLines[i] = 1
    for v in DSR.getVars():
        if v.varName.split('[')[0] == 'x_Line':
            ResultLines[v.varName.split('[')[1][:-1]] = v.x

    # Faulted Lines
    for i in FaultedLines:
        ResultLines[i] = 0

    for j, k in ResultLines.items():
        Temp_Bus1 = PDElementsConnections[j][0]
        Temp_Bus2 = PDElementsConnections[j][1]
        if (XY_Coordinates[Temp_Bus1][0] != 0) and (XY_Coordinates[Temp_Bus1][1] != 0):
            if (XY_Coordinates[Temp_Bus2][0] != 0) and (XY_Coordinates[Temp_Bus2][1] != 0):
                TempX = [XY_Coordinates[Temp_Bus1][0], XY_Coordinates[Temp_Bus2][0]]
                TempY = [XY_Coordinates[Temp_Bus1][1], XY_Coordinates[Temp_Bus2][1]]
                NumPhases = sum(PDElementPhaseMatrix[j])
                if k == 0:
                    if SwitchableLines[j] == 1:
                        plt.plot(TempX, TempY, 'y-', marker=6, linewidth=1, markersize=1.5, markevery=100) # marker size was 2
                    else:
                        if NumPhases == 3:
                            plt.plot(TempX, TempY, 'r-', linewidth=0.8)
                        elif NumPhases == 2:
                            #plt.plot(TempX, TempY, 'r--', linewidth=0.8)
                            plt.plot(TempX, TempY, 'r', linewidth=0.8)
                        else:
                            #plt.plot(TempX, TempY, 'r:', linewidth=0.8)
                            plt.plot(TempX, TempY, 'r-', linewidth=0.8)
                        #plt.scatter(XY_Coordinates[Temp_Bus1][0], XY_Coordinates[Temp_Bus1][1], c='k', marker='x', linewidths=2)
                        plt.scatter(TempX, TempY, c='k', marker='x', linewidths=2)
                        plt.text(TempX, TempY, 'Fault', fontsize=5)
                        

                else:
                    if SwitchableLines[j] == 1:
                        plt.plot(TempX, TempY, 'k-', marker=7, linewidth=1, markersize=2, markevery=100)#marker size was 3


                    else:
                        if NumPhases == 3:
                            plt.plot(TempX, TempY, 'k-', linewidth=0.7)
                        elif NumPhases == 2:
                            plt.plot(TempX, TempY, 'k--', linewidth=0.4)
                        else:
                            plt.plot(TempX, TempY, 'k--', linewidth=0.4)

    plt.xlabel('X')
    plt.ylabel('Y')

    plt.savefig('IEEE_Network.pdf')
    #plt.savefig('IEEE_Network.png')
    plt.show()

# *******************************************************************************************************
# * This function is used to build Lagrange multipliers (mu) of the ADMM method
# *******************************************************************************************************
def Build_mu(Circuit, RestorationTime):
     FlowDirection = ['d1', 'd2']
     mu = {}
     initial = 0 #May change to 1 **************************
     iteration = 0  # initial iteration

     # Initial iteration
     mu[iteration] = {}

     for t in RestorationTime:
         mu[iteration][t] = {}
         for CL, buses in Circuit.Clusters.items():
             mu[iteration][t][CL] = {}

             # Voltage (Lagrange multiplier requires for parent bus of the cluster and buses that connects to the children buses
             mu[iteration][t][CL]['Voltage'] = {}
             for Bus, NeighborBuses in Circuit.ClustersNeighborBuses[CL].items():
                 for NeighborBus in NeighborBuses:
                     if 'ParentNeighborBuses' in Circuit.ClustersNeighborRelations[CL]:
                         if NeighborBus in Circuit.ClustersNeighborRelations[CL]['ParentNeighborBuses']:
                             mu[iteration][t][CL]['Voltage'][NeighborBus] = {}
                             for phase in Circuit.PhaseSequence:
                                 #if Circuit.BusesPhaseMatrix[NeighborBus][Circuit.PhaseSequence.index(phase)] == 1:
                                 mu[iteration][t][CL]['Voltage'][NeighborBus][phase] = initial
                         else:
                            mu[iteration][t][CL]['Voltage'][Bus] = {}
                            for phase in Circuit.PhaseSequence:
                                mu[iteration][t][CL]['Voltage'][Bus][phase] = initial
                     else:
                         mu[iteration][t][CL]['Voltage'][Bus] = {}
                         for phase in Circuit.PhaseSequence:
                             mu[iteration][t][CL]['Voltage'][Bus][phase] = initial

             # Current
             #mu[iteration][t][CL]['l_Line'] = {}
             #for PDElement in Circuit.ClustersBoundaryPDElements[CL]:
             #    mu[iteration][t][CL]['l_Line'][PDElement] = {}
             #    for phase in Circuit.PhaseSequence:
             #        mu[iteration][t][CL]['l_Line'][PDElement][phase] = initial

             # Line power
             mu[iteration][t][CL]['P_Line'] = {}
             mu[iteration][t][CL]['Q_Line'] = {}
             for PDElement in Circuit.ClustersBoundaryPDElements[CL]:
                 mu[iteration][t][CL]['P_Line'][PDElement] = {}
                 mu[iteration][t][CL]['Q_Line'][PDElement] = {}
                 for phase in Circuit.PhaseSequence:
                     mu[iteration][t][CL]['P_Line'][PDElement][phase] = initial
                     mu[iteration][t][CL]['Q_Line'][PDElement][phase] = initial

             # Switchable Lines Status
             mu[iteration][t][CL]['x_Switch'] = {}
             for Switch in Circuit.ClustersSwitchableLines[CL]:
                 mu[iteration][t][CL]['x_Switch'][Switch] = initial

             # Loads Status
             mu[iteration][t][CL]['x_Load'] = {}
             for Load in Circuit.ClustersLoadsNonDispatchable[CL]:
                 if Circuit.LoadType[Load] == 0:
                     mu[iteration][t][CL]['x_Load'][Load] = initial

             # Bus Status
             mu[iteration][t][CL]['x_Bus'] = {}
             for bus in Circuit.ClusterBusVariables[CL]:
                 mu[iteration][t][CL]['x_Bus'][bus] = initial

             # Line Flow Direction
             mu[iteration][t][CL]['Flow'] = {}
             for PDElement in Circuit.ClustersPDElementsReduced[CL]:
                 mu[iteration][t][CL]['Flow'][PDElement] = {}
                 for flow in FlowDirection:
                     mu[iteration][t][CL]['Flow'][PDElement][flow] = initial

     return mu

# *******************************************************************************************************
# * This function is used to build X matrix of the ADMM method
# *******************************************************************************************************
def Build_X(Circuit, RestorationTime):
    X = {}
    initial = 0
    for t in RestorationTime:
        X[t] = {}
        for CL, buses in Circuit.Clusters.items():
            X[t][CL] = {}

            # Load Power
            X[t][CL]['P_Load'] = {}
            X[t][CL]['Q_Load'] = {}
            for load in Circuit.ClustersLoads[CL]:
                X[t][CL]['P_Load'][load] = initial
                X[t][CL]['Q_Load'][load] = initial

            # Voltage
            X[t][CL]['Voltage'] = {}
            for bus in Circuit.ClusterBusVariables[CL]:
                X[t][CL]['Voltage'][bus] = {}
                for phase in Circuit.PhaseSequence:
                    X[t][CL]['Voltage'][bus][phase] = initial

            # Current
            #X[t][CL]['l_Line'] = {}
            #for PDElement in Circuit.ClusterCurrentVariables[CL]:
            #    X[t][CL]['l_Line'][PDElement] = {}
            #    for phase in Circuit.PhaseSequence:
            #        X[t][CL]['l_Line'][PDElement][phase] = initial

            # Line Power
            X[t][CL]['P_Line'] = {}
            X[t][CL]['Q_Line'] = {}
            for PDElement in Circuit.ClustersPDElements[CL]:
                X[t][CL]['P_Line'][PDElement] = {}
                X[t][CL]['Q_Line'][PDElement] = {}
                for phase in Circuit.PhaseSequence:
                    X[t][CL]['P_Line'][PDElement][phase] = initial
                    X[t][CL]['Q_Line'][PDElement][phase] = initial

            # Substation Power
            X[t][CL]['P_sub'] = {}
            X[t][CL]['Q_sub'] = {}
            for SourceBus in Circuit.ClustersSourceBus[CL]:
                X[t][CL]['P_sub'][SourceBus] = {}
                X[t][CL]['Q_sub'][SourceBus] = {}
                for phase in Circuit.PhaseSequence:
                    X[t][CL]['P_sub'][SourceBus][phase] = {}
                    X[t][CL]['Q_sub'][SourceBus][phase] = {}

            # Capacitor bank
            X[t][CL]['Q_cap'] = {}
            for cap in Circuit.ClustersCBs[CL]:
                X[t][CL]['Q_cap'][cap] = {}
                for phase in Circuit.PhaseSequence:
                    X[t][CL]['Q_cap'][cap][phase] = initial

            # PV
            X[t][CL]['P_PV'] = {}
            X[t][CL]['Q_PV'] = {}
            for PVName in Circuit.ClustersPVs[CL]:
                X[t][CL]['P_PV'][PVName] = {}
                X[t][CL]['Q_PV'][PVName] = {}
                for phase in Circuit.PhaseSequence:
                    X[t][CL]['P_PV'][PVName][phase] = initial
                    X[t][CL]['Q_PV'][PVName][phase] = initial

            # Total Generation and Load of the Cluster
            X[t][CL]['P_LoadTotal'] = initial
            X[t][CL]['P_GenTotal'] = initial

    return X
# *******************************************************************************************************
# * This function is used to build Y matrix of the ADMM method
# *******************************************************************************************************
def Build_Y(Circuit, RestorationTime):
    FlowDirection = ['d1', 'd2']
    Y = {}
    initial = 0
    for t in RestorationTime:
        Y[t] = {}
        for CL, buses in Circuit.Clusters.items():
            Y[t][CL] = {}

            # Loads Status
            Y[t][CL]['x_Load'] = {}
            for Load in Circuit.ClustersLoadsNonDispatchable[CL]:
                if Circuit.LoadType[Load] == 0:
                    Y[t][CL]['x_Load'][Load] = initial

            # Bus Status
            Y[t][CL]['x_Bus'] = {}
            for bus in Circuit.ClusterBusVariables[CL]:
                Y[t][CL]['x_Bus'][bus] = initial

            # Switch Status (x_Line)
            Y[t][CL]['x_Switch'] = {}
            for Switch in Circuit.ClustersSwitchableLines[CL]:
                Y[t][CL]['x_Switch'][Switch] = initial

            # Line Flow Direction
            Y[t][CL]['Flow'] = {}
            for PDElement in Circuit.ClustersPDElementsReduced[CL]:
                Y[t][CL]['Flow'][PDElement] = {}
                for flow in FlowDirection:
                    Y[t][CL]['Flow'][PDElement][flow] = initial

    return Y

# *******************************************************************************************************
# * This function is used to build X_bar matrix of the ADMM method
# *******************************************************************************************************
def Build_X_bar(Circuit, mu, RestorationTime, PhaseSequence):
    X_bar = {}
    initial = 0
    iteration = 0

    #Initial iteration
    X_bar[iteration] = {}

    for t in RestorationTime:
        X_bar[iteration][t] = {}
        for CL, buses in Circuit.Clusters.items():
            X_bar[iteration][t][CL] = {}

            # Voltage (Here voltage variable includes one for each boundary buses xcept the one which connects to
            # the parent cluster in that case consider the parent bus in the parent cluster)
            X_bar[iteration][t][CL]['Voltage'] = {}
            for bus in mu[iteration][t][CL]['Voltage'].keys():
                X_bar[iteration][t][CL]['Voltage'][bus] = {}
                for phase in PhaseSequence:
                    X_bar[iteration][t][CL]['Voltage'][bus][phase] = initial

            # Current

            # Line power
            X_bar[iteration][t][CL]['P_Line'] = {}
            X_bar[iteration][t][CL]['Q_Line'] = {}
            for PDElement in mu[iteration][t][CL]['P_Line']: #Circuit.ClustersBoundaryPDElements[CL]:
                X_bar[iteration][t][CL]['P_Line'][PDElement] = {}
                X_bar[iteration][t][CL]['Q_Line'][PDElement] = {}
                for phase in PhaseSequence:
                    X_bar[iteration][t][CL]['P_Line'][PDElement][phase] = initial
                    X_bar[iteration][t][CL]['Q_Line'][PDElement][phase] = initial

    return X_bar

# *******************************************************************************************************
# * This function is used to build Z matrix of the ADMM method
# *******************************************************************************************************
def Build_Z(Circuit, mu, RestorationTime, initial):
    FlowDirection = ['d1', 'd2']
    Z = {}
    iteration = 0
    initial_loads = 0
    # Initial iteration
    Z[iteration] = {}

    for t in RestorationTime:
        Z[iteration][t] = {}
        for CL, buses in Circuit.Clusters.items():
            Z[iteration][t][CL] = {}

            # Loads Status
            Z[iteration][t][CL]['x_Load'] = {}
            for Load in Circuit.ClustersLoadsNonDispatchable[CL]:
                if Circuit.LoadType[Load] == 0:
                    Z[iteration][t][CL]['x_Load'][Load] = initial_loads

            # Bus Status
            Z[iteration][t][CL]['x_Bus'] = {}
            for bus in Circuit.ClusterBusVariables[CL]:
                Z[iteration][t][CL]['x_Bus'][bus] = initial

            # Switchable Lines Status
            Z[iteration][t][CL]['x_Switch'] = {}
            for Switch in mu[iteration][t][CL]['x_Switch']:  # Circuit.ClustersBoundarySwitchableLines[CL]:
                Z[iteration][t][CL]['x_Switch'][Switch] = initial

            # Line Flow Direction
            Z[iteration][t][CL]['Flow'] = {}
            for PDElement in Circuit.ClustersPDElementsReduced[CL]:
                Z[iteration][t][CL]['Flow'][PDElement] = {}
                for flow in FlowDirection:
                    Z[iteration][t][CL]['Flow'][PDElement][flow] = initial


    return Z

# *******************************************************************************************************
# * This function is used to manipulated some values of the Z matrix for the just Polish phase
# *******************************************************************************************************
def Manipulate_Z(Z, Circuit, RestorationTime, Opened_Switches):
    iteration = 0
    if len(Opened_Switches):
        for OpenSwitch in Opened_Switches:
            for CL, Switches in Circuit.ClustersSwitchableLines.items():
                if OpenSwitch in Switches:
                    for t in RestorationTime:
                        Z[iteration][t][CL]['x_Switch'][OpenSwitch] = 0
    return Z

# *******************************************************************************************************
# * This function is used to update X and Y matrix of the ADMM method
# *******************************************************************************************************
def Update_X_and_Y(DSR, X, Y, ClusterIndex):

    for v in DSR.getVars():
        VarFullName = v.VarName
        VarName = VarFullName.split('[')[0]

        # Load Power----------------------
        # Active Power
        if VarName == 'P_Load':
            LoadName = VarFullName.split('[')[1].split(',')[0]
            TimeStep = VarFullName.split('[')[1].split(',')[1][:-1]
            X[TimeStep][ClusterIndex][VarName][LoadName] = v.X
        # Reactive Power
        elif VarName == 'Q_Load':
            LoadName = VarFullName.split('[')[1].split(',')[0]
            TimeStep = VarFullName.split('[')[1].split(',')[1][:-1]
            X[TimeStep][ClusterIndex][VarName][LoadName] = v.X

        # Voltage -------------------------
        elif VarName == 'Voltage':
            phase = VarFullName.split('[')[1].split(',')[0]
            bus = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            X[TimeStep][ClusterIndex][VarName][bus][phase] = v.X

        # Line Power------------------------
        # Active Power
        elif VarName == 'P_Line':
            phase = VarFullName.split('[')[1].split(',')[0]
            PDElement = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            X[TimeStep][ClusterIndex][VarName][PDElement][phase] = v.X
        # Reactive Power
        elif VarName == 'Q_Line':
            phase = VarFullName.split('[')[1].split(',')[0]
            PDElement = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            X[TimeStep][ClusterIndex][VarName][PDElement][phase] = v.X

        # Substation Power-----------------
        # Active Power
        elif VarName == 'P_sub':
            phase = VarFullName.split('[')[1].split(',')[0]
            SourceBus = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            X[TimeStep][ClusterIndex][VarName][SourceBus][phase] = v.X
        # Active Power
        elif VarName == 'Q_sub':
            phase = VarFullName.split('[')[1].split(',')[0]
            SourceBus = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            X[TimeStep][ClusterIndex][VarName][SourceBus][phase] = v.X

        # Capacitor bank-----------------------
        elif VarName == 'Q_cap':
            phase = VarFullName.split('[')[1].split(',')[0]
            cap = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            X[TimeStep][ClusterIndex][VarName][cap][phase] = v.X

        # PV------------------------------------------
        # Active Power
        elif VarName == 'P_PV':
            phase = VarFullName.split('[')[1].split(',')[0]
            PV = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            X[TimeStep][ClusterIndex][VarName][PV][phase] = v.X
        # Reactive Power
        elif VarName == 'Q_PV':
            phase = VarFullName.split('[')[1].split(',')[0]
            PV = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            X[TimeStep][ClusterIndex][VarName][PV][phase] = v.X

        # Total Generation and Load of the Cluster------
        # Total Load
        elif VarName == 'P_LoadTotal':
            TimeStep = VarFullName.split('[')[1][:-1]
            X[TimeStep][ClusterIndex][VarName] = v.X
        # Total Generation
        elif VarName == 'P_GenTotal':
            TimeStep = VarFullName.split('[')[1][:-1]
            X[TimeStep][ClusterIndex][VarName] = v.X

        # Switch Status (x_Line)-----------------------
        elif VarName == 'x_Switch':
            switch = VarFullName.split('[')[1].split(',')[0]
            TimeStep = VarFullName.split('[')[1].split(',')[1][:-1]
            Y[TimeStep][ClusterIndex][VarName][switch] = v.X

        # Loads Status---------------------------
        if VarName == 'x_Load':
            LoadName = VarFullName.split('[')[1].split(',')[0]
            TimeStep = VarFullName.split('[')[1].split(',')[1][:-1]
            Y[TimeStep][ClusterIndex][VarName][LoadName] = v.X

        # Bus Status ----------------------------
        elif VarName == 'x_Bus':
            bus = VarFullName.split('[')[1].split(',')[0]
            TimeStep = VarFullName.split('[')[1].split(',')[1][:-1]
            Y[TimeStep][ClusterIndex][VarName][bus] = v.X

        # Line Flow Direction ---------------------
        elif VarName == 'Flow':
            PDElement = VarFullName.split('[')[1].split(',')[0]
            flow = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            Y[TimeStep][ClusterIndex][VarName][PDElement][flow] = v.X

    return [X, Y]

# *******************************************************************************************************
# * This function is used to update X matrix of the ADMM method during Polish method
# *******************************************************************************************************
def Update_X_Polish(DSR, X, ClusterIndex):

    for v in DSR.getVars():
        VarFullName = v.VarName
        VarName = VarFullName.split('[')[0]

        # Load Power----------------------
        # Active Power
        if VarName == 'P_Load':
            LoadName = VarFullName.split('[')[1].split(',')[0]
            TimeStep = VarFullName.split('[')[1].split(',')[1][:-1]
            X[TimeStep][ClusterIndex][VarName][LoadName] = v.X
        # Reactive Power
        elif VarName == 'Q_Load':
            LoadName = VarFullName.split('[')[1].split(',')[0]
            TimeStep = VarFullName.split('[')[1].split(',')[1][:-1]
            X[TimeStep][ClusterIndex][VarName][LoadName] = v.X

        # Voltage -------------------------
        elif VarName == 'Voltage':
            phase = VarFullName.split('[')[1].split(',')[0]
            bus = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            X[TimeStep][ClusterIndex][VarName][bus][phase] = v.X

        # Line Power------------------------
        # Active Power
        elif VarName == 'P_Line':
            phase = VarFullName.split('[')[1].split(',')[0]
            PDElement = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            X[TimeStep][ClusterIndex][VarName][PDElement][phase] = v.X
        # Reactive Power
        elif VarName == 'Q_Line':
            phase = VarFullName.split('[')[1].split(',')[0]
            PDElement = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            X[TimeStep][ClusterIndex][VarName][PDElement][phase] = v.X

        # Substation Power-----------------
        # Active Power
        elif VarName == 'P_sub':
            phase = VarFullName.split('[')[1].split(',')[0]
            SourceBus = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            X[TimeStep][ClusterIndex][VarName][SourceBus][phase] = v.X
        # Active Power
        elif VarName == 'Q_sub':
            phase = VarFullName.split('[')[1].split(',')[0]
            SourceBus = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            X[TimeStep][ClusterIndex][VarName][SourceBus][phase] = v.X

        # Capacitor bank-----------------------
        elif VarName == 'Q_cap':
            phase = VarFullName.split('[')[1].split(',')[0]
            cap = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            X[TimeStep][ClusterIndex][VarName][cap][phase] = v.X

        # PV------------------------------------------
        # Active Power
        elif VarName == 'P_PV':
            phase = VarFullName.split('[')[1].split(',')[0]
            PV = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            X[TimeStep][ClusterIndex][VarName][PV][phase] = v.X
        # Reactive Power
        elif VarName == 'Q_PV':
            phase = VarFullName.split('[')[1].split(',')[0]
            PV = VarFullName.split('[')[1].split(',')[1]
            TimeStep = VarFullName.split('[')[1].split(',')[2][:-1]
            X[TimeStep][ClusterIndex][VarName][PV][phase] = v.X

        # Total Generation and Load of the Cluster------
        # Total Load
        elif VarName == 'P_LoadTotal':
            TimeStep = VarFullName.split('[')[1][:-1]
            X[TimeStep][ClusterIndex][VarName] = v.X
        # Total Generation
        elif VarName == 'P_GenTotal':
            TimeStep = VarFullName.split('[')[1][:-1]
            X[TimeStep][ClusterIndex][VarName] = v.X

# *******************************************************************************************************
# * This function is used to update loads in Z matrix of the ADMM method during MIP Polish method
# *******************************************************************************************************
def Update_Z_Polish(DSR, Z, ClusterIndex, Last_iter):
    for v in DSR.getVars():
        VarFullName = v.VarName
        VarName = VarFullName.split('[')[0]

        # Load Status----------------------
        if VarName == 'x_Load':
            LoadName = VarFullName.split('[')[1].split(',')[0]
            TimeStep = VarFullName.split('[')[1].split(',')[1][:-1]
            Z[Last_iter][TimeStep][ClusterIndex][VarName][LoadName] = v.X
    return Z

# *******************************************************************************************************
# * This function is used to update X_bar matrix of the ADMM method
# *******************************************************************************************************
def Update_X_bar(X_bar, X, mu, iteration, ClustersNeighborBuses, BusesInClusters, ClustersBoundaryPDElements,
                 ClustersNeighborRelations, PhaseSequence):

    X_bar[iteration] = {}

    for t in X_bar[iteration - 1]:
        X_bar[iteration][t] = {}
        for CL in X_bar[iteration - 1][t]:
            X_bar[iteration][t][CL] = {}

            # Voltage (Update requires for parent bus of the cluster and buses that connects to the children buses)====
            X_bar[iteration][t][CL]['Voltage'] = {}
            for bus in mu[iteration-1][t][CL]['Voltage']:
                if 'ParentNeighborBuses' in ClustersNeighborRelations[CL]:
                    if bus in ClustersNeighborRelations[CL]['ParentNeighborBuses']:
                        X_bar[iteration][t][CL]['Voltage'][bus] = {}
                        ParentCluster = BusesInClusters[bus]
                        ChildClusters = []
                        for child_bus in ClustersNeighborBuses[ParentCluster][bus]:
                            if 'ChildrenNeighborBuses' in ClustersNeighborRelations[ParentCluster]:
                                if child_bus in ClustersNeighborRelations[ParentCluster]['ChildrenNeighborBuses']:
                                    ChildClusters.append(BusesInClusters[child_bus])
                        NumChildCluster = len(ChildClusters)
                        if NumChildCluster != 0:
                            X_bar[iteration][t][CL]['Voltage'][bus] = {}
                            for phase in PhaseSequence:
                                X_bar[iteration][t][CL]['Voltage'][bus][phase] = (
                                                            X[t][ParentCluster]['Voltage'][bus][phase] +
                                                            mu[iteration-1][t][ParentCluster]['Voltage'][bus][phase])
                                for NeighborCluster in ChildClusters:
                                    X_bar[iteration][t][CL]['Voltage'][bus][phase] = (
                                        X_bar[iteration][t][CL]['Voltage'][bus][phase] +
                                            X[t][NeighborCluster]['Voltage'][bus][phase] +
                                            mu[iteration - 1][t][NeighborCluster]['Voltage'][bus][phase])
                                X_bar[iteration][t][CL]['Voltage'][bus][phase] = (
                                                X_bar[iteration][t][CL]['Voltage'][bus][phase])/(NumChildCluster + 1)

                    else:
                        ChildClusters = []
                        for child_bus in ClustersNeighborBuses[CL][bus]:
                            if 'ChildrenNeighborBuses' in ClustersNeighborRelations[CL]:
                                if child_bus in ClustersNeighborRelations[CL]['ChildrenNeighborBuses']:
                                    ChildClusters.append(BusesInClusters[child_bus])
                        NumChildClusters = len(ChildClusters)
                        if NumChildClusters != 0:
                            X_bar[iteration][t][CL]['Voltage'][bus] = {}
                            for phase in PhaseSequence:
                                X_bar[iteration][t][CL]['Voltage'][bus][phase] = (X[t][CL]['Voltage'][bus][phase] +
                                                                    mu[iteration - 1][t][CL]['Voltage'][bus][phase])
                                for Neighborcluster in ChildClusters:
                                    X_bar[iteration][t][CL]['Voltage'][bus][phase] = (
                                        X_bar[iteration][t][CL]['Voltage'][bus][phase] +
                                        X[t][Neighborcluster]['Voltage'][bus][phase] +
                                        mu[iteration - 1][t][Neighborcluster]['Voltage'][bus][phase])
                                X_bar[iteration][t][CL]['Voltage'][bus][phase] = (
                                    X_bar[iteration][t][CL]['Voltage'][bus][phase])/(NumChildClusters+1)
                else:
                    ChildClusters = []
                    for child_bus in ClustersNeighborBuses[CL][bus]:
                        if 'ChildrenNeighborBuses' in ClustersNeighborRelations[CL]:
                            if child_bus in ClustersNeighborRelations[CL]['ChildrenNeighborBuses']:
                                ChildClusters.append(BusesInClusters[child_bus])
                    NumChildClusters = len(ChildClusters)
                    if NumChildClusters != 0:
                        X_bar[iteration][t][CL]['Voltage'][bus] = {}
                        for phase in PhaseSequence:
                            X_bar[iteration][t][CL]['Voltage'][bus][phase] = (X[t][CL]['Voltage'][bus][phase] +
                                                                    mu[iteration - 1][t][CL]['Voltage'][bus][phase])
                            for Neighborcluster in ChildClusters:
                                X_bar[iteration][t][CL]['Voltage'][bus][phase] = (
                                        X_bar[iteration][t][CL]['Voltage'][bus][phase] +
                                        X[t][Neighborcluster]['Voltage'][bus][phase] +
                                        mu[iteration - 1][t][Neighborcluster]['Voltage'][bus][phase])
                            X_bar[iteration][t][CL]['Voltage'][bus][phase] = (
                                        X_bar[iteration][t][CL]['Voltage'][bus][phase]) / (NumChildClusters + 1)

            # Line power=============================================================
            X_bar[iteration][t][CL]['P_Line'] = {}
            X_bar[iteration][t][CL]['Q_Line'] = {}
            for PDElement in X_bar[iteration-1][t][CL]['P_Line']:
                X_bar[iteration][t][CL]['P_Line'][PDElement] = {}
                X_bar[iteration][t][CL]['Q_Line'][PDElement] = {}
                for phase in PhaseSequence:
                    #Active Power ------
                    X_bar[iteration][t][CL]['P_Line'][PDElement][phase] = 0.5*(X[t][CL]['P_Line'][PDElement][phase] +
                                                                mu[iteration-1][t][CL]['P_Line'][PDElement][phase] +
                    X[t][ClustersBoundaryPDElements[CL][PDElement]]['P_Line'][PDElement][phase] +
                    mu[iteration-1][t][ClustersBoundaryPDElements[CL][PDElement]]['P_Line'][PDElement][phase])
                    #Reactive Power -----
                    X_bar[iteration][t][CL]['Q_Line'][PDElement][phase] = 0.5 * (X[t][CL]['Q_Line'][PDElement][phase] +
                                                                mu[iteration-1][t][CL]['Q_Line'][PDElement][phase] +
                                    X[t][ClustersBoundaryPDElements[CL][PDElement]]['Q_Line'][PDElement][phase] +
                            mu[iteration-1][t][ClustersBoundaryPDElements[CL][PDElement]]['Q_Line'][PDElement][phase])

    return X_bar

# *******************************************************************************************************
# * This function is used to update Z matrix of the ADMM method
# *******************************************************************************************************
def Update_Z(Z, Y, mu, iteration, Binary_Projection, ClustersNeighborRelations, BusesInClusters,
             ClustersBoundaryPDElements, ClustersBoundarySwitchableLines, ClustersNeighborBuses, Projected_ADMM_trigger,
             Prox_Oper, t_Prox):
    Z[iteration] = {}
    FlowDirection = ['d1', 'd2']
    for t in Z[iteration-1]:
        Z[iteration][t] = {}
        for CL in Z[iteration-1][t]:
            Z[iteration][t][CL] = {}

            # Loads Status----------------------------------------------------
            Z[iteration][t][CL]['x_Load'] = {}
            if Projected_ADMM_trigger == 1:
                for load in Z[iteration-1][t][CL]['x_Load']:
                    Z[iteration][t][CL]['x_Load'][load] = Binary_Projection(Y[t][CL]['x_Load'][load] +
                                                               mu[iteration-1][t][CL]['x_Load'][load])

            elif Projected_ADMM_trigger == 0:
                for load in Z[iteration-1][t][CL]['x_Load']:
                    W_W = (Y[t][CL]['x_Load'][load] + mu[iteration-1][t][CL]['x_Load'][load])
                    Z[iteration][t][CL]['x_Load'][load] = Prox_Oper(W_W, Binary_Projection(W_W), t_Prox)

            # Bus Status--------------------------------------------------------
            Z[iteration][t][CL]['x_Bus'] = {}
            if Projected_ADMM_trigger == 1:
                for bus in Z[iteration-1][t][CL]['x_Bus']:
                    Z[iteration][t][CL]['x_Bus'][bus] = Binary_Projection(Y[t][CL]['x_Bus'][bus] +
                                                                          mu[iteration - 1][t][CL]['x_Bus'][bus])
            elif Projected_ADMM_trigger == 0:
                for bus in Z[iteration-1][t][CL]['x_Bus']:
                    W_W = (Y[t][CL]['x_Bus'][bus] + mu[iteration - 1][t][CL]['x_Bus'][bus])
                    Z[iteration][t][CL]['x_Bus'][bus] = Prox_Oper(W_W, Binary_Projection(W_W), t_Prox)

            # The x_bus is not considered as consensus variable in the new version.
            # In fact, this augumented term is just used to force this become binary. therefore, all code is commented!!
            '''if Projected_ADMM_trigger == 1:
                for bus in Z[iteration-1][t][CL]['x_Bus']:
                    # if bus is the neighbor buses
                    if bus in ClustersNeighborBuses[CL]:
                        if 'ParentNeighborBuses' in ClustersNeighborRelations[CL]:
                            if bus in ClustersNeighborRelations[CL]['ParentNeighborBuses']:
                                Z[iteration][t][CL]['x_Bus'][bus] = Binary_Projection(0.5*(Y[t][CL]['x_Bus'][bus] +
                                                                            mu[iteration-1][t][CL]['x_Bus'][bus] +
                                                Y[t][BusesInClusters[bus]]['x_Bus'][bus] +
                                                mu[iteration-1][t][BusesInClusters[bus]]['x_Bus'][bus]))

                            else:
                                for childbus in ClustersNeighborRelations[CL]['ChildrenNeighborBuses']:
                                    Z[iteration][t][CL]['x_Bus'][bus] = Binary_Projection(0.5*(Y[t][CL]['x_Bus'][bus] +
                                                                            mu[iteration-1][t][CL]['x_Bus'][bus] +
                                    Y[t][BusesInClusters[ClustersNeighborBuses[CL][bus]]]['x_Bus'][bus] +
                                    mu[iteration-1][t][BusesInClusters[ClustersNeighborBuses[CL][bus]]]['x_Bus'][bus]))

                        else:
                            Z[iteration][t][CL]['x_Bus'][bus] = Binary_Projection(0.5 * (Y[t][CL]['x_Bus'][bus] +
                                                                            mu[iteration - 1][t][CL]['x_Bus'][bus] +
                                    Y[t][BusesInClusters[ClustersNeighborBuses[CL][bus]]]['x_Bus'][bus] +
                                    mu[iteration - 1][t][BusesInClusters[ClustersNeighborBuses[CL][bus]]]['x_Bus'][bus]))

                    else:
                        Z[iteration][t][CL]['x_Bus'][bus] = Binary_Projection(Y[t][CL]['x_Bus'][bus] +
                                                                              mu[iteration - 1][t][CL]['x_Bus'][bus])

            elif Projected_ADMM_trigger == 0:
                for bus in Z[iteration - 1][t][CL]['x_Bus']:
                    # if bus is the neighbor buses
                    if bus in ClustersNeighborBuses[CL]:
                        if 'ParentNeighborBuses' in ClustersNeighborRelations[CL]:
                            if bus in ClustersNeighborRelations[CL]['ParentNeighborBuses']:
                                Z[iteration][t][CL]['x_Bus'][bus] = 0.5*(Y[t][CL]['x_Bus'][bus] +
                                           mu[iteration - 1][t][CL]['x_Bus'][bus] +
                                           Y[t][BusesInClusters[bus]]['x_Bus'][bus] +
                                           mu[iteration - 1][t][BusesInClusters[bus]]['x_Bus'][bus])

                            else:
                                 Z[iteration][t][CL]['x_Bus'][bus] = 0.5*(Y[t][CL]['x_Bus'][bus] +
                                                                 mu[iteration - 1][t][CL]['x_Bus'][bus] +
                                    Y[t][BusesInClusters[ClustersNeighborBuses[CL][bus]]]['x_Bus'][bus] +
                                    mu[iteration - 1][t][BusesInClusters[ClustersNeighborBuses[CL][bus]]]['x_Bus'][bus])

                        else:
                            Z[iteration][t][CL]['x_Bus'][bus] = 0.5 * (Y[t][CL]['x_Bus'][bus] +
                                                                          mu[iteration - 1][t][CL]['x_Bus'][bus] +
                                     Y[t][BusesInClusters[ClustersNeighborBuses[CL][bus]]]['x_Bus'][bus] +
                                     mu[iteration - 1][t][BusesInClusters[ClustersNeighborBuses[CL][bus]]]['x_Bus'][bus])

                    else:
                        Z[iteration][t][CL]['x_Bus'][bus] = Z[iteration-1][t][CL]['x_Bus'][bus]'''


            # Switchable Lines Status=========================
            Z[iteration][t][CL]['x_Switch'] = {}
            if Projected_ADMM_trigger == 1:
                for Switch in Z[iteration - 1][t][CL]['x_Switch']:
                    # if switch is boundary one:
                    if Switch in ClustersBoundarySwitchableLines[CL]:
                        Z[iteration][t][CL]['x_Switch'][Switch] = Binary_Projection(
                        0.5*(Y[t][CL]['x_Switch'][Switch] + mu[iteration - 1][t][CL]['x_Switch'][Switch] +
                             Y[t][ClustersBoundaryPDElements[CL][Switch]]['x_Switch'][Switch] +
                             mu[iteration - 1][t][ClustersBoundaryPDElements[CL][Switch]]['x_Switch'][Switch]))
                    # If switch is inside the cluster
                    else:
                        Z[iteration][t][CL]['x_Switch'][Switch] = Binary_Projection(
                            Y[t][CL]['x_Switch'][Switch] + mu[iteration - 1][t][CL]['x_Switch'][Switch])

            elif Projected_ADMM_trigger == 0:
                for Switch in Z[iteration - 1][t][CL]['x_Switch']:
                    # if switch is boundary one:
                    if Switch in ClustersBoundarySwitchableLines[CL]:
                        W_W = 0.5 * (Y[t][CL]['x_Switch'][Switch] + mu[iteration - 1][t][CL]['x_Switch'][Switch] +
                                   Y[t][ClustersBoundaryPDElements[CL][Switch]]['x_Switch'][Switch] +
                                   mu[iteration - 1][t][ClustersBoundaryPDElements[CL][Switch]]['x_Switch'][Switch])

                        Z[iteration][t][CL]['x_Switch'][Switch] = Prox_Oper(W_W, Binary_Projection(W_W), t_Prox)
                    # If switch is inside the cluster
                    else:
                        W_W = (Y[t][CL]['x_Switch'][Switch] + mu[iteration - 1][t][CL]['x_Switch'][Switch])

                        Z[iteration][t][CL]['x_Switch'][Switch] = Prox_Oper(W_W, Binary_Projection(W_W), t_Prox)

            # Line Flow Direction-------------------
            Z[iteration][t][CL]['Flow'] = {}
            if Projected_ADMM_trigger == 1:
                for PDElement in Z[iteration-1][t][CL]['Flow']:
                    Z[iteration][t][CL]['Flow'][PDElement] = {}
                    # If PDElement is a boundary one
                    if PDElement in ClustersBoundaryPDElements[CL]:
                        for flow in FlowDirection:
                            Z[iteration][t][CL]['Flow'][PDElement][flow] = Binary_Projection(
                            0.5*(Y[t][CL]['Flow'][PDElement][flow] + mu[iteration-1][t][CL]['Flow'][PDElement][flow] +
                            Y[t][ClustersBoundaryPDElements[CL][PDElement]]['Flow'][PDElement][flow] +
                        mu[iteration-1][t][ClustersBoundaryPDElements[CL][PDElement]]['Flow'][PDElement][flow]))
                    # If PDElement is inside the cluster and not the boundary one
                    else:
                        for flow in FlowDirection:
                            Z[iteration][t][CL]['Flow'][PDElement][flow] = Binary_Projection(
                                Y[t][CL]['Flow'][PDElement][flow] + mu[iteration - 1][t][CL]['Flow'][PDElement][flow])
            #-----------
            elif Projected_ADMM_trigger == 0:
                for PDElement in Z[iteration-1][t][CL]['Flow']:
                    Z[iteration][t][CL]['Flow'][PDElement] = {}
                    # If PDElement is a boundary one
                    if PDElement in ClustersBoundaryPDElements[CL]:
                        for flow in FlowDirection:
                            W_W = 0.5 * (Y[t][CL]['Flow'][PDElement][flow] + mu[iteration - 1][t][CL]['Flow'][PDElement][flow] +
                            Y[t][ClustersBoundaryPDElements[CL][PDElement]]['Flow'][PDElement][flow] +
                            mu[iteration - 1][t][ClustersBoundaryPDElements[CL][PDElement]]['Flow'][PDElement][flow])

                            Z[iteration][t][CL]['Flow'][PDElement][flow] = Prox_Oper(W_W, Binary_Projection(W_W), t_Prox)
                    # If PDElement is inside the cluster and not the boundary one
                    else:
                        for flow in FlowDirection:
                            W_W = (Y[t][CL]['Flow'][PDElement][flow] + mu[iteration - 1][t][CL]['Flow'][PDElement][flow])

                            Z[iteration][t][CL]['Flow'][PDElement][flow] = Prox_Oper(W_W, Binary_Projection(W_W), t_Prox)

    return Z

# *******************************************************************************************************
# * This function is used to update Lagrange multipliers or mu matrix of the ADMM method
# *******************************************************************************************************
def Update_mu(mu, X, X_bar, Y, Z, iteration, Polish_status):
    mu[iteration] = {}
    FlowDirection = ['d1', 'd2']
    # time steps cycle
    for t in mu[iteration-1]:
        mu[iteration][t] = {}
        # Cycling over clusters
        for CL in mu[iteration-1][t]:
            mu[iteration][t][CL] = {}

            # Voltage (Lagrange multiplier requires for parent bus of the cluster and buses that connects to the children buses)
            mu[iteration][t][CL]['Voltage'] = {}
            for bus in mu[iteration-1][t][CL]['Voltage']:
                mu[iteration][t][CL]['Voltage'][bus] = {}
                for phase in mu[iteration-1][t][CL]['Voltage'][bus]:
                    mu[iteration][t][CL]['Voltage'][bus][phase] = mu[iteration-1][t][CL]['Voltage'][bus][phase] + (
                    X[t][CL]['Voltage'][bus][phase] - X_bar[iteration][t][CL]['Voltage'][bus][phase])

            # Line power----------------------------------------------------------
            mu[iteration][t][CL]['P_Line'] = {}
            mu[iteration][t][CL]['Q_Line'] = {}
            for PDElement in mu[iteration-1][t][CL]['P_Line']:
                mu[iteration][t][CL]['P_Line'][PDElement] = {}
                mu[iteration][t][CL]['Q_Line'][PDElement] = {}
                for phase in mu[iteration-1][t][CL]['P_Line'][PDElement]:
                    # Active power
                    mu[iteration][t][CL]['P_Line'][PDElement][phase] = mu[iteration-1][t][CL]['P_Line'][PDElement][phase] + (
                    X[t][CL]['P_Line'][PDElement][phase] - X_bar[iteration][t][CL]['P_Line'][PDElement][phase])

                    # Reactive power
                    mu[iteration][t][CL]['Q_Line'][PDElement][phase] = mu[iteration-1][t][CL]['Q_Line'][PDElement][phase] + (
                    X[t][CL]['Q_Line'][PDElement][phase] - X_bar[iteration][t][CL]['Q_Line'][PDElement][phase])

            if not Polish_status:
                # Loads Status-------------------------------------------------------
                mu[iteration][t][CL]['x_Load'] = {}
                for load in mu[iteration-1][t][CL]['x_Load']:
                    mu[iteration][t][CL]['x_Load'][load] = mu[iteration-1][t][CL]['x_Load'][load] + (
                    Y[t][CL]['x_Load'][load] - Z[iteration][t][CL]['x_Load'][load])

                # Bus Status------------------------------------------------------------
                mu[iteration][t][CL]['x_Bus'] = {}
                for bus in mu[iteration-1][t][CL]['x_Bus']:
                    mu[iteration][t][CL]['x_Bus'][bus] = mu[iteration-1][t][CL]['x_Bus'][bus] + (
                        Y[t][CL]['x_Bus'][bus] - Z[iteration][t][CL]['x_Bus'][bus])

                # Switchable Lines Status---------------------------------------------
                mu[iteration][t][CL]['x_Switch'] = {}
                for switch in mu[iteration - 1][t][CL]['x_Switch']:
                    mu[iteration][t][CL]['x_Switch'][switch] = mu[iteration - 1][t][CL]['x_Switch'][switch] + (
                      Y[t][CL]['x_Switch'][switch] - Z[iteration][t][CL]['x_Switch'][switch])

                # Line Flow Direction---------------------------------------------------
                mu[iteration][t][CL]['Flow'] = {}
                for PDElement in mu[iteration-1][t][CL]['Flow']:
                    mu[iteration][t][CL]['Flow'][PDElement] = {}
                    for flow in FlowDirection:
                        mu[iteration][t][CL]['Flow'][PDElement][flow] = mu[iteration-1][t][CL]['Flow'][PDElement][flow] + (
                        Y[t][CL]['Flow'][PDElement][flow] - Z[iteration][t][CL]['Flow'][PDElement][flow])

    return mu
# *******************************************************************************************************
# * This function is used as Proximal operator
# *******************************************************************************************************
def Prox_Oper(W, Proj, t):
    Z = (W + t*Proj) / (1+t)
    return Z

# *******************************************************************************************************
# * This function is used to calculate Primal residual of the ADMM method
# *******************************************************************************************************
def Calculate_PrimalResidual(X, Y, X_bar, Z, iteration, mu, sqrt, Projected_ADMM_trigger,
                             ClustersBoundarySwitchableLines, ClustersBoundaryPDElements, Polish_status,
                             Binary_Residual, Oscillating_Vars):
    residual = 0

    for t in mu[iteration]:
        for CL in mu[iteration][t]:

            if not Binary_Residual:
                # Voltage-------------------------------------------------------------
                for bus in mu[iteration][t][CL]['Voltage']:
                    for phase in mu[iteration][t][CL]['Voltage'][bus]:
                        residual += (X[t][CL]['Voltage'][bus][phase] -
                                     X_bar[iteration][t][CL]['Voltage'][bus][phase])**2

                # Line power----------------------------------------------------------
                for PDElement in mu[iteration][t][CL]['P_Line']:
                    for phase in mu[iteration][t][CL]['P_Line'][PDElement]:
                        # Active power
                        residual += (X[t][CL]['P_Line'][PDElement][phase] -
                                     X_bar[iteration][t][CL]['P_Line'][PDElement][phase])**2

                        # Reactive power
                        residual += (X[t][CL]['Q_Line'][PDElement][phase] -
                                     X_bar[iteration][t][CL]['Q_Line'][PDElement][phase])**2

            # Loads Status--------------------------------------------------------
            if not Polish_status:
                if Projected_ADMM_trigger == 1:
                    for load in mu[iteration][t][CL]['x_Load']:
                        residual += (Y[t][CL]['x_Load'][load] -
                                     Z[iteration][t][CL]['x_Load'][load])**2

            # Bus Status------------------------------------------------------------
            if not Polish_status:
                if Projected_ADMM_trigger == 1:
                    for bus in mu[iteration][t][CL]['x_Bus']:
                        if bus not in Oscillating_Vars:
                            residual += (Y[t][CL]['x_Bus'][bus] -
                                     Z[iteration][t][CL]['x_Bus'][bus])**2

            # Switchable Lines Status---------------------------------------------
            if not Polish_status:
                if Projected_ADMM_trigger == 1:
                    for switch in mu[iteration][t][CL]['x_Switch']:
                        residual += (Y[t][CL]['x_Switch'][switch] -
                                     Z[iteration][t][CL]['x_Switch'][switch]) ** 2
                elif Projected_ADMM_trigger == 0:
                    for switch in mu[iteration][t][CL]['x_Switch']:
                        #if Switch is boundary
                        if switch in ClustersBoundarySwitchableLines[CL]:
                            residual += (Y[t][CL]['x_Switch'][switch] -
                                         Z[iteration][t][CL]['x_Switch'][switch]) ** 2

            # Line Flow Direction---------------------------------------------------
            if not Polish_status:
                if Projected_ADMM_trigger == 1:
                    for PDElement in mu[iteration][t][CL]['Flow']:
                        for flow in mu[iteration][t][CL]['Flow'][PDElement]:
                            residual += (Y[t][CL]['Flow'][PDElement][flow] -
                                         Z[iteration][t][CL]['Flow'][PDElement][flow])**2
                elif Projected_ADMM_trigger == 0:
                    for PDElement in mu[iteration][t][CL]['Flow']:
                        if PDElement in ClustersBoundaryPDElements[CL]:
                            for flow in mu[iteration][t][CL]['Flow'][PDElement]:
                                residual += (Y[t][CL]['Flow'][PDElement][flow] -
                                             Z[iteration][t][CL]['Flow'][PDElement][flow]) ** 2

    return sqrt(residual)

# *******************************************************************************************************
# * This function is used to calculate Dual residual of the ADMM method
# *******************************************************************************************************
def Calculate_DualResidual(X_bar, Z, iteration, mu, sqrt, rho, Projected_ADMM_trigger, ClustersBoundarySwitchableLines,
                           ClustersBoundaryPDElements, Polish_status, Binary_Residual, Oscillating_Vars):
    residual = 0
    Plot_Circuit
    for t in mu[iteration]:
        for CL in mu[iteration][t]:

            if not Binary_Residual:
                # Voltage-------------------------------------------------------------
                for bus in mu[iteration][t][CL]['Voltage']:
                    for phase in mu[iteration][t][CL]['Voltage'][bus]:
                        residual += (rho*(X_bar[iteration][t][CL]['Voltage'][bus][phase] -
                                          X_bar[iteration-1][t][CL]['Voltage'][bus][phase])) ** 2

                # Line power----------------------------------------------------------
                for PDElement in mu[iteration][t][CL]['P_Line']:
                    for phase in mu[iteration][t][CL]['P_Line'][PDElement]:
                        # Active power
                        residual += (rho*(X_bar[iteration][t][CL]['P_Line'][PDElement][phase] -
                                          X_bar[iteration-1][t][CL]['P_Line'][PDElement][phase])) ** 2

                        # Reactive power
                        residual += (rho*(X_bar[iteration][t][CL]['Q_Line'][PDElement][phase] -
                                          X_bar[iteration-1][t][CL]['Q_Line'][PDElement][phase])) ** 2

            # Loads Status--------------------------------------------------------
            if not Polish_status:
                if Projected_ADMM_trigger == 1:
                    for load in mu[iteration][t][CL]['x_Load']:
                        residual += (rho*(Z[iteration][t][CL]['x_Load'][load] -
                                          Z[iteration-1][t][CL]['x_Load'][load])) ** 2

            # Bus Status------------------------------------------------------------
            if not Polish_status:
                if Projected_ADMM_trigger == 1:
                    for bus in mu[iteration][t][CL]['x_Bus']:
                        if bus not in Oscillating_Vars:
                            residual += (rho*(Z[iteration][t][CL]['x_Bus'][bus] -
                                          Z[iteration-1][t][CL]['x_Bus'][bus])) ** 2

            # Switchable Lines Status---------------------------------------------
            if not Polish_status:
                if Projected_ADMM_trigger == 1:
                    for switch in mu[iteration][t][CL]['x_Switch']:
                        residual += (rho * (Z[iteration][t][CL]['x_Switch'][switch] -
                                            Z[iteration - 1][t][CL]['x_Switch'][switch])) ** 2
                elif Projected_ADMM_trigger == 0:
                    for switch in mu[iteration][t][CL]['x_Switch']:
                        #if Switch is boundary
                        if switch in ClustersBoundarySwitchableLines[CL]:
                            residual += (rho * (Z[iteration][t][CL]['x_Switch'][switch] -
                                                Z[iteration - 1][t][CL]['x_Switch'][switch])) ** 2

            # Line Flow Direction---------------------------------------------------
            if not Polish_status:
                if Projected_ADMM_trigger == 1:
                    for PDElement in mu[iteration][t][CL]['Flow']:
                        for flow in mu[iteration][t][CL]['Flow'][PDElement]:
                            residual += (rho*(Z[iteration][t][CL]['Flow'][PDElement][flow] -
                                              Z[iteration-1][t][CL]['Flow'][PDElement][flow])) ** 2
                elif Projected_ADMM_trigger == 0:
                    for PDElement in mu[iteration][t][CL]['Flow']:
                        if PDElement in ClustersBoundaryPDElements[CL]:
                            for flow in mu[iteration][t][CL]['Flow'][PDElement]:
                                residual += (rho * (Z[iteration][t][CL]['Flow'][PDElement][flow] -
                                                    Z[iteration - 1][t][CL]['Flow'][PDElement][flow])) ** 2

    return sqrt(residual)

# *******************************************************************************************************
# * This function is used to plot Primal and Dual residuals of the ADMM method
# *******************************************************************************************************
def Plot_Residuals(PrimalResidual_All, DualResidual_All, iteration, plt):

    iterations_list = list(range(iteration+1))
    plt.figure(2)
    plot1, = plt.plot(iterations_list, PrimalResidual_All, 'k-') # Primal plot
    plot2, = plt.plot(iterations_list, DualResidual_All, 'b-') # Dual plot
    plt.legend([plot1, plot2], ['Primal Residual', 'Dual Residual'])

    plt.grid(True)
    plt.xlabel('Iteration')
    plt.ylabel('Residuals')
    plt.savefig('Residual.pdf')
    plt.ylim([0,0.4])

    plt.show()

# *******************************************************************************************************
# * This function is used to plot proximal operator weight
# *******************************************************************************************************
def Plot_t_Prox(t_Prox, iteration, plt, c1):
    iterations_list = list(range(iteration + 1))
    plt.figure(3)
    plt.plot(iterations_list, t_Prox)

    plt.grid(True)
    plt.xlabel('Iteration')
    plt.ylabel('t_Prox')
    plt.title('t_prox(c1={})'.format(c1))

    plt.show()

# *******************************************************************************************************
# * This function is used to plot total loads
# *******************************************************************************************************
def Plot_TotalLoad(X, Z, iteration, Last_iter, plt, np, RestorationTime, S_base, Loads_max, LoadType):
    P_TotalLoad = []
    y_pos = np.arange(len(RestorationTime))
    for t in RestorationTime:
        P_TotalLoad.append(0)

    for t in RestorationTime:
        for CL in X[t]:
            for load in X[t][CL]['P_Load']:
                if LoadType[load]:
                    P_TotalLoad[RestorationTime.index(t)] += X[t][CL]['P_Load'][load]*S_base
                else:
                    if Z[Last_iter][t][CL]['x_Load'][load]:
                        P_TotalLoad[RestorationTime.index(t)] += Loads_max[load][0]*S_base

            # Previous code
            #P_TotalLoad[RestorationTime.index(t)] += X[t][CL]['P_LoadTotal']*S_base

    plt.figure(4)
    plt.bar(y_pos, P_TotalLoad, align='center', alpha=0.5)
    plt.grid(True)
    plt.xticks(y_pos, RestorationTime)
    plt.ylabel('Total power')
    plt.title('Total restored Loads')

    plt.show()

    return P_TotalLoad

# *******************************************************************************************************
# * This function is used to plot the whole circuit
# *******************************************************************************************************
def Plot_Circuit(X, Y, X_bar, Z, plt, Circuit, TimeStep, iteration, Last_iter, FaultedLines):
    plt.figure(5)
    # Draw Load pickup-----------------
    Result_Loads = {}
    for CL in Z[Last_iter][TimeStep]:
        for load in Z[Last_iter][TimeStep][CL]['x_Load']:
            Result_Loads[load] = Z[Last_iter][TimeStep][CL]['x_Load'][load]

    # Default values :
    count = 0
    X_Coordinates_Previous = 0
    Y_Coordinates_Previous = 0
    # ------------------------------
    for load, pickup in Result_Loads.items():
        LoadBus = Circuit.LoadsConnection[load][0]
        X_Coordinates = Circuit.XY_Coordinates[LoadBus][0]
        Y_Coordinates = Circuit.XY_Coordinates[LoadBus][1]
        if (X_Coordinates != 0) and (Y_Coordinates != 0):
            if (X_Coordinates == X_Coordinates_Previous) and (Y_Coordinates == Y_Coordinates_Previous):
                count += 10
            else:
                count = 0
            if pickup == 0:
                plt.scatter(X_Coordinates+count, Y_Coordinates+count, c='r', s=10)
            elif pickup == 1:
                plt.scatter(X_Coordinates+count, Y_Coordinates+count, c='g', s=10)

        X_Coordinates_Previous = X_Coordinates
        Y_Coordinates_Previous = Y_Coordinates

    # Draw Buses-----------------------
    for bus in Circuit.AllBuses:
        if (Circuit.XY_Coordinates[bus][0] != 0) and (Circuit.XY_Coordinates[bus][1] != 0):
            if Circuit.PVBus[bus] != []:
                plt.scatter(Circuit.XY_Coordinates[bus][0], Circuit.XY_Coordinates[bus][1], c='k', marker='D', s=30)
                plt.text(Circuit.XY_Coordinates[bus][0] - 300, Circuit.XY_Coordinates[bus][1]-1200 , 'PV', fontsize=10)

            if bus in Circuit.SourceBus:
                plt.scatter(Circuit.XY_Coordinates[bus][0], Circuit.XY_Coordinates[bus][1], c='k', marker='*', linewidths=4)
                plt.text(Circuit.XY_Coordinates[bus][0] - 500, Circuit.XY_Coordinates[bus][1] + 600, 'Sub.', fontsize=10)

            # Bus Name:
            #plt.text(Circuit.XY_Coordinates[bus][0]-10, Circuit.XY_Coordinates[bus][1] - 100, '{}'.format(bus), fontsize=5)

            #Important loads **************
            if bus in Circuit.LoadToBus:
                if Circuit.LoadToBus[bus] != []:
                    for load in Circuit.LoadToBus[bus]:
                        if Circuit.LoadsWeight[load] > 1:
                            plt.text(Circuit.XY_Coordinates[bus][0]-10, Circuit.XY_Coordinates[bus][1] - 180, '**', fontsize=8)



    # Draw energized Lines----------------------------------------------------------
    result_line = {}
    for PDElement in Circuit.AllPDElements:
        result_line[PDElement] = 1

    for CL in Z[Last_iter][TimeStep]:
        for Switch in Z[Last_iter][TimeStep][CL]['x_Switch']:
                result_line[Switch] = Z[Last_iter][TimeStep][CL]['x_Switch'][Switch]

    # Faulted Lines
    for Fault_Line in FaultedLines:
        result_line[Fault_Line] = 0

    for PDElement, status in result_line.items():
        Temp_Bus1 = Circuit.PDElementsConnections[PDElement][0]
        Temp_Bus2 = Circuit.PDElementsConnections[PDElement][1]
        if (Circuit.XY_Coordinates[Temp_Bus1][0] != 0) and (Circuit.XY_Coordinates[Temp_Bus1][1] != 0):
            if (Circuit.XY_Coordinates[Temp_Bus2][0] != 0) and (Circuit.XY_Coordinates[Temp_Bus2][1] != 0):
                TempX = [Circuit.XY_Coordinates[Temp_Bus1][0], Circuit.XY_Coordinates[Temp_Bus2][0]]
                TempY = [Circuit.XY_Coordinates[Temp_Bus1][1], Circuit.XY_Coordinates[Temp_Bus2][1]]
                NumPhases = sum(Circuit.PDElementPhaseMatrix[PDElement])
                if status == 0:
                    if Circuit.SwitchableLines[PDElement] == 1:
                        plt.plot(TempX, TempY, 'y-', marker=6, linewidth=1, markersize=8, markevery=100)
                    else:
                        plt.text(Circuit.XY_Coordinates[Temp_Bus1][0]-200, Circuit.XY_Coordinates[Temp_Bus1][1]+1000, 'Fault', fontsize=10)
                        if NumPhases == 3:
                            plt.plot(TempX, TempY, 'r-', marker='x', linewidth=0.8, markersize=4, markevery=100)
                        elif NumPhases == 2:
                             # plt.plot(TempX, TempY, 'r--', linewidth=0.8)
                            plt.plot(TempX, TempY, 'r', marker='x', linewidth=0.8, markersize=4, markevery=100)
                        else:
                             # plt.plot(TempX, TempY, 'r:', linewidth=0.8)
                            plt.plot(TempX, TempY, 'r-', marker='x', linewidth=0.8, markersize=4, markevery=100)
                        #plt.text(TempX, TempY, 'Fault', fontsize=5)
                        
                        
                elif status == 1:
                    if Circuit.SwitchableLines[PDElement] == 1:
                        plt.plot(TempX, TempY, 'b-', marker=7, linewidth=1, markersize=8, markevery=100)
                    else:
                        if NumPhases == 3:
                            plt.plot(TempX, TempY, 'k-', linewidth=1)
                        elif NumPhases == 2:
                            plt.plot(TempX, TempY, 'k--', linewidth=0.9)
                        else:
                            plt.plot(TempX, TempY, 'k--', linewidth=0.9)

    plt.xlabel('X')
    plt.ylabel('Y')

    plt.savefig('IEEE_Network.pdf')
    plt.savefig('IEEE_Network.png')
    plt.show()

# *******************************************************************************************************



# *******************************************************************************************************
# * This function is used to find error of the code in Primal residual
# *******************************************************************************************************
def Find_Points_Primal(X, Y, X_bar, Z, iteration, mu, sqrt, Projected_ADMM_trigger, Last_iter,
                             ClustersBoundarySwitchableLines, ClustersBoundaryPDElements):

    for t in mu[iteration]:
        for CL in mu[iteration][t]:

            # Voltage-------------------------------------------------------------
            for bus in mu[iteration][t][CL]['Voltage']:
                for phase in mu[iteration][t][CL]['Voltage'][bus]:
                    residual = sqrt((X[t][CL]['Voltage'][bus][phase] -
                                 X_bar[iteration][t][CL]['Voltage'][bus][phase])**2)
                    if residual >= 0.5:
                        print('Primal Residual','Voltage', bus)

            # Line power----------------------------------------------------------
            for PDElement in mu[iteration][t][CL]['P_Line']:
                for phase in mu[iteration][t][CL]['P_Line'][PDElement]:
                    # Active power
                    residual = sqrt((X[t][CL]['P_Line'][PDElement][phase] -
                                 X_bar[iteration][t][CL]['P_Line'][PDElement][phase])**2)
                    if residual >= 0.5:
                        print('Primal Residual','P_Line', PDElement)

                    # Reactive power
                    residual = sqrt((X[t][CL]['Q_Line'][PDElement][phase] -
                                 X_bar[iteration][t][CL]['Q_Line'][PDElement][phase])**2)

                    if residual >= 0.5:
                        print('Primal Residual','Q_Line', PDElement)

            # Loads Status--------------------------------------------------------
            if Projected_ADMM_trigger == 1:
                for load in mu[Last_iter][t][CL]['x_Load']:
                    residual = sqrt((Y[t][CL]['x_Load'][load] -
                                 Z[Last_iter][t][CL]['x_Load'][load])**2)

                    if residual >= 0.5:
                        print('Primal Residual','x_Load', load)

            # Bus Status------------------------------------------------------------
            if Projected_ADMM_trigger == 1:
                for bus in mu[Last_iter][t][CL]['x_Bus']:
                    residual = sqrt((Y[t][CL]['x_Bus'][bus] -
                                 Z[Last_iter][t][CL]['x_Bus'][bus])**2)

                    if residual >= 0.5:
                        print('Primal Residual','x_Bus', bus)

            # Switchable Lines Status---------------------------------------------
            if Projected_ADMM_trigger == 1:
                for switch in mu[Last_iter][t][CL]['x_Switch']:
                    residual = sqrt((Y[t][CL]['x_Switch'][switch] -
                                 Z[Last_iter][t][CL]['x_Switch'][switch]) ** 2)

                    if residual >= 0.5:
                        print('Primal Residual','x_Switch', switch)
            elif Projected_ADMM_trigger == 0:
                for switch in mu[Last_iter][t][CL]['x_Switch']:
                    #if Switch is boundary
                    if switch in ClustersBoundarySwitchableLines[CL]:
                        residual = sqrt((Y[t][CL]['x_Switch'][switch] -
                                     Z[Last_iter][t][CL]['x_Switch'][switch]) ** 2)

            # Line Flow Direction---------------------------------------------------
            if Projected_ADMM_trigger == 1:
                for PDElement in mu[Last_iter][t][CL]['Flow']:
                    for flow in mu[Last_iter][t][CL]['Flow'][PDElement]:
                        residual = sqrt((Y[t][CL]['Flow'][PDElement][flow] -
                                     Z[Last_iter][t][CL]['Flow'][PDElement][flow])**2)

                        if residual >= 0.5:
                            print('Primal Residual','Flow', PDElement)

            elif Projected_ADMM_trigger == 0:
                for PDElement in mu[Last_iter][t][CL]['Flow']:
                    if PDElement in ClustersBoundaryPDElements[CL]:
                        for flow in mu[Last_iter][t][CL]['Flow'][PDElement]:
                            residual = sqrt((Y[t][CL]['Flow'][PDElement][flow] -
                                         Z[Last_iter][t][CL]['Flow'][PDElement][flow]) ** 2)


# *******************************************************************************************************
# * This function is used to find error of the code in Dual residual
# *******************************************************************************************************
def Find_Points_Dual(X_bar, Z, iteration, mu, sqrt, rho, Projected_ADMM_trigger, Last_iter,
                     ClustersBoundarySwitchableLines, ClustersBoundaryPDElements):

    for t in mu[iteration]:
        for CL in mu[iteration][t]:

            # Voltage-------------------------------------------------------------
            for bus in mu[iteration][t][CL]['Voltage']:
                for phase in mu[iteration][t][CL]['Voltage'][bus]:
                    residual = sqrt((rho*(X_bar[iteration][t][CL]['Voltage'][bus][phase] -
                                      X_bar[iteration-1][t][CL]['Voltage'][bus][phase])) ** 2)
                    if residual >= 0.5:
                        print('Dual Residual', 'Voltage', bus)

            # Line power----------------------------------------------------------
            for PDElement in mu[iteration][t][CL]['P_Line']:
                for phase in mu[iteration][t][CL]['P_Line'][PDElement]:
                    # Active power
                    residual = sqrt((rho*(X_bar[iteration][t][CL]['P_Line'][PDElement][phase] -
                                      X_bar[iteration-1][t][CL]['P_Line'][PDElement][phase])) ** 2)

                    if residual >= 0.5:
                        print('Dual Residual', 'P_Line', PDElement)

                    # Reactive power
                    residual = sqrt((rho*(X_bar[iteration][t][CL]['Q_Line'][PDElement][phase] -
                                      X_bar[iteration-1][t][CL]['Q_Line'][PDElement][phase])) ** 2)
                    if residual >= 0.5:
                        print('Dual Residual', 'Q_Line', PDElement)

            # Loads Status--------------------------------------------------------
            if Projected_ADMM_trigger == 1:
                for load in mu[Last_iter][t][CL]['x_Load']:
                    residual = sqrt((rho*(Z[Last_iter][t][CL]['x_Load'][load] -
                                      Z[Last_iter-1][t][CL]['x_Load'][load])) ** 2)

                    if residual >= 0.5:
                        print('Dual Residual', 'x_Load', load)

            # Bus Status------------------------------------------------------------
            if Projected_ADMM_trigger == 1:
                for bus in mu[Last_iter][t][CL]['x_Bus']:
                    residual = sqrt((rho*(Z[Last_iter][t][CL]['x_Bus'][bus] -
                                      Z[Last_iter-1][t][CL]['x_Bus'][bus])) ** 2)
                    if residual >= 0.5:
                        print('Dual Residual', 'x_Bus', bus)

            # Switchable Lines Status---------------------------------------------
            if Projected_ADMM_trigger == 1:
                for switch in mu[Last_iter][t][CL]['x_Switch']:
                    residual = sqrt((rho * (Z[Last_iter][t][CL]['x_Switch'][switch] -
                                        Z[Last_iter - 1][t][CL]['x_Switch'][switch])) ** 2)

                    if residual >= 0.5:
                        print('Dual Residual', 'x_Switch', switch)
            elif Projected_ADMM_trigger == 0:
                for switch in mu[Last_iter][t][CL]['x_Switch']:
                    #if Switch is boundary
                    if switch in ClustersBoundarySwitchableLines[CL]:
                        residual = sqrt((rho * (Z[Last_iter][t][CL]['x_Switch'][switch] -
                                            Z[Last_iter - 1][t][CL]['x_Switch'][switch])) ** 2)

            # Line Flow Direction---------------------------------------------------
            if Projected_ADMM_trigger == 1:
                for PDElement in mu[Last_iter][t][CL]['Flow']:
                    for flow in mu[Last_iter][t][CL]['Flow'][PDElement]:
                        residual = sqrt((rho*(Z[Last_iter][t][CL]['Flow'][PDElement][flow] -
                                          Z[Last_iter-1][t][CL]['Flow'][PDElement][flow])) ** 2)

                        if residual >= 0.5:
                            print('Dual Residual', 'Flow', PDElement)
            elif Projected_ADMM_trigger == 0:
                for PDElement in mu[Last_iter][t][CL]['Flow']:
                    if PDElement in ClustersBoundaryPDElements[CL]:
                        for flow in mu[Last_iter][t][CL]['Flow'][PDElement]:
                            residual = sqrt((rho * (Z[Last_iter][t][CL]['Flow'][PDElement][flow] -
                                                Z[Last_iter - 1][t][CL]['Flow'][PDElement][flow])) ** 2)


# *******************************************************************************************************
# * This function is used to save the evolution procedure of some variables during ADMM convergence
# *******************************************************************************************************
def Build_VariableSave(RestorationTime, X, Y, Z):
    SavedItems = {}


    # Total Load
    SavedItems['P_LoadTotal'] = {}
    for t in RestorationTime:
        SavedItems['P_LoadTotal'][t] = [0]

    # Total Load_Z
    SavedItems['P_LoadTotal_Z'] = {}
    for t in RestorationTime:
        SavedItems['P_LoadTotal_Z'][t] = [0]

    # Total Load_Y
    SavedItems['P_LoadTotal_Y'] = {}
    for t in RestorationTime:
        SavedItems['P_LoadTotal_Y'][t] = [0]

    # Total Load by each cluster
    SavedItems['P_LoadCluster_Z'] = {}
    SavedItems['P_LoadCluster_Y'] = {}
    for t in RestorationTime:
        SavedItems['P_LoadCluster_Z'][t] = {}
        SavedItems['P_LoadCluster_Y'][t] = {}
        for CL in X[RestorationTime[0]]:
            SavedItems['P_LoadCluster_Z'][t][CL] = [0]
            SavedItems['P_LoadCluster_Y'][t][CL] = [0]


    # Special Loads saved data
    # s53 load
    # SavedItems['Load_s53a_Z'] = {}
    # SavedItems['Load_s53a_Y'] = {}
    # for t in RestorationTime:
    #     SavedItems['Load_s53a_Z'][t] = [0]
    #     SavedItems['Load_s53a_Y'][t] = [0]

    # Switches Operation
    SavedItems['x_Switch'] = {}
    for CL in X[RestorationTime[0]]:
        for Switch in Y[RestorationTime[0]][CL]['x_Switch']:
            SavedItems['x_Switch'][Switch] = {}
            for t in RestorationTime:
                SavedItems['x_Switch'][Switch][t] = [Y[t][CL]['x_Switch'][Switch]]

    # Switches Operation for Z variables
    SavedItems['x_Switch_Z'] = {}
    for CL in X[RestorationTime[0]]:
        for Switch in Y[RestorationTime[0]][CL]['x_Switch']:
            SavedItems['x_Switch_Z'][Switch] = {}
            for t in RestorationTime:
                SavedItems['x_Switch_Z'][Switch][t] = [Z[0][t][CL]['x_Switch'][Switch]]

    return SavedItems
# **************************************

def ADMM_VariableSave(SavedItems, X, Y, Z, iteration, Last_iter, LoadType, Loads_max, S_base, Polish_status,
                      RestorationTime, LoadsConnection, BusesInClusters):

    for t in RestorationTime:

        # Total Network load
        P_TotalLoad = 0
        for CL in X[t]:
            for load in X[t][CL]['P_Load']:
                if LoadType[load]:
                    P_TotalLoad += X[t][CL]['P_Load'][load]*S_base
                else:
                    if Polish_status:
                        if (Z[Last_iter][t][CL]['x_Load'][load]) == 1:
                            P_TotalLoad += Loads_max[load][0]*S_base
                    else:
                        if (Z[iteration][t][CL]['x_Load'][load]) == 1:
                            P_TotalLoad += Loads_max[load][0]*S_base
        SavedItems['P_LoadTotal'][t].append(P_TotalLoad)

        # Total Network load _ Z relaxed
        P_TotalLoad2 = 0
        for CL in X[t]:
            P_TotalCluster_Z = 0
            for load in X[t][CL]['P_Load']:
                if LoadType[load]:
                    P_TotalLoad2 += X[t][CL]['P_Load'][load] * S_base
                    P_TotalCluster_Z += X[t][CL]['P_Load'][load] * S_base
                else:
                    if Polish_status:
                        P_TotalLoad2 += Z[Last_iter][t][CL]['x_Load'][load]*Loads_max[load][0] * S_base
                        P_TotalCluster_Z += Z[Last_iter][t][CL]['x_Load'][load]*Loads_max[load][0] * S_base
                    else:
                        P_TotalLoad2 += Z[iteration][t][CL]['x_Load'][load]*Loads_max[load][0] * S_base
                        P_TotalCluster_Z += Z[iteration][t][CL]['x_Load'][load]*Loads_max[load][0] * S_base
            SavedItems['P_LoadCluster_Z'][t][CL].append(P_TotalCluster_Z)
        SavedItems['P_LoadTotal_Z'][t].append(P_TotalLoad2)

        # Total Network load_Y and each cluster total load
        P_TotalLoad1 = 0
        for CL in X[t]:
            P_TotalCluster_Y = 0
            for load in X[t][CL]['P_Load']:
                if LoadType[load]:
                    P_TotalLoad1 += X[t][CL]['P_Load'][load] * S_base
                    P_TotalCluster_Y += X[t][CL]['P_Load'][load] * S_base
                else:
                    if Polish_status:
                        P_TotalLoad1 += Y[t][CL]['x_Load'][load]*Loads_max[load][0] * S_base
                        P_TotalCluster_Y += Y[t][CL]['x_Load'][load]*Loads_max[load][0] * S_base
                    else:
                        P_TotalLoad1 += Y[t][CL]['x_Load'][load]*Loads_max[load][0] * S_base
                        P_TotalCluster_Y += Y[t][CL]['x_Load'][load]*Loads_max[load][0] * S_base
            SavedItems['P_LoadCluster_Y'][t][CL].append(P_TotalCluster_Y)
        SavedItems['P_LoadTotal_Y'][t].append(P_TotalLoad1)

        # Special Loads saved data
        # # s53 load
        # LoadName1 = 's53a'
        # Load1_Cluster = BusesInClusters[LoadsConnection[LoadName1][0]]
        # # Z relaxed
        # if Polish_status:
        #     SavedItems['Load_s53a_Z'][t].append(Z[Last_iter][t][Load1_Cluster]['x_Load'][LoadName1]*Loads_max[LoadName1][0] * S_base)
        # else:
        #     SavedItems['Load_s53a_Z'][t].append(
        #         Z[iteration][t][Load1_Cluster]['x_Load'][LoadName1] * Loads_max[LoadName1][0] * S_base)
        #
        # # Y
        # SavedItems['Load_s53a_Y'][t].append(Y[t][Load1_Cluster]['x_Load'][LoadName1]*Loads_max[LoadName1][0] * S_base)


        # Switches Operation
        for CL in X[t]:
            for Switch in Y[t][CL]['x_Switch']:
                SavedItems['x_Switch'][Switch][t].append(Y[t][CL]['x_Switch'][Switch])

        # Switches Operation_Z
        for CL in X[t]:
            for Switch in Y[t][CL]['x_Switch']:
                if Polish_status:
                    SavedItems['x_Switch_Z'][Switch][t].append(Z[Last_iter][t][CL]['x_Switch'][Switch])
                else:
                    SavedItems['x_Switch_Z'][Switch][t].append(Z[iteration][t][CL]['x_Switch'][Switch])
    return SavedItems



