import win32com.client

class DSS():

    def __init__(self):
        self.dssObj = win32com.client.Dispatch("OpenDSSEngine.DSS")

        # Initialize DSS object
        if self.dssObj.start(0) == False:
            print("Unable to start the OpenDSS Engine")

        else:
            self.dssText = self.dssObj.Text
            self.dssCircuit = self.dssObj.ActiveCircuit
            self.dssSolution = self.dssCircuit.Solution
            self.dssElement = self.dssCircuit.ActiveCktElement
            self.dssBus = self.dssCircuit.ActiveBus

    # ***********************************************************************
    def Version(self):
        return self.dssObj.Version

    # ***********************************************************************
    def Compile_DSS(self):
        self.dssObj.ClearAll()

    # ***********************************************************************
    def Solve(self):
        self.dssSolution.Solve()

    # ***********************************************************************
    def IncMatrix(self):
        self.dssText.command = 'CalcIncMatrix'
        return self.dssCircuit.Solution.IncMatrix

    # ***********************************************************************
    def AllBuses(self):
        return self.dssCircuit.AllBusNames

    # ***********************************************************************
    def FindAllPDElements(self):
        return self.dssCircuit.Solution.IncMatrixrows

    # ***********************************************************************
    # This method is used to find base power of the circuit in terms of the first transformer
    def FindS_Base(self):
        self.dssCircuit.Transformers.First
        return self.dssCircuit.Transformers.kva

    # ***********************************************************************
    # This method is used to find base voltage of the primary feeder
    def FindV_Base(self):
        self.dssCircuit.Transformers.First
        return self.dssCircuit.Transformers.kV

    # ***********************************************************************
    # This method is used to find source bus of the substation
    def FindSourceBus(self):
        self.__Source = []
        self.dssCircuit.Vsources.First
        self.dssCircuit.SetActiveElement(self.dssCircuit.Vsources.Name)
        self.__Source.append(self.dssCircuit.ActiveCktElement.BusNames[0].split('.')[0])
        return self.__Source


