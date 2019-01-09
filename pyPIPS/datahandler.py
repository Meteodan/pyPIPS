

def getDataHandler(model_name, base_dir, times, microphys=None, multitime=True):
    from wrfmodule import WRFDataHandler
    from commasmodule import COMMASDataHandler
    from arpsmodule import ARPSDataHandler

    if model_name == "WRF":
        return WRFDataHandler(base_dir, times)
    elif model_name == "COMMAS":
        return COMMASDataHandler(base_dir, times, multitime=multitime)
    elif model_name == "ARPS":
        return ARPSDataHandler(base_dir, times, microphys)

class DataHandler(object):
    def __init__(self, model_name):
        self._model_name = model_name
        self._fields = {}
        self._constants = {}
        return

    def setRun(self):
        self._abstract('setRun')
        return

    def setTime(self):
        self._abstract('setTime')
        return

    def loadGrid(self):
        self._abstract('loadGrid')
        return

    def loadMicrophysics(self):
        self._abstract('loadMicrophysics')
        return

    def loadMeanWind(self):
        self._abstract('loadMeanWind')
        return

    def loadHorizWind(self):
        self._abstract('loadHorizWind')
        return

    def loadVertWind(self):
        self._abstract('loadVertWind')
        return

    def fileNameBuilder(self):
        self._abstract('fileNameBuilder')
        return

    def dumpToFile(self, format='npz'):
        if format not in [ 'npz', 'nc' ]:
            raise ValueError("Keyword format ('%s') must be one of 'npz' or 'nc'" % format)

        if format == 'npz':
            self._dumpNPZ()
        elif format == 'nc':
            self._dumpNC()
        return

    def slice(self, coord, axis='x'):
        if axis not in [ 'x', 'y', 'z', 'e' ]:
            raise ValueError("Keyword axis ('%s') must be one of 'x', 'y', 'z', or 'e'" % axis)

        return

    def __getitem__(self, key):
        if key in self._fields:
            return self._fields[key]
        elif key in self._constants:
            return self._constants[key]
        else:
            raise KeyError("Variable '%s' not in DataHandler object." % key)
        return

    def _dumpNPZ(self):
        return

    def _dumpNC(self):
        return

    def _abstract(self, func_name):
        raise NotImplementedError("'%s()' is an abstract function that needs to be overridden in a subclass!" % func_name)
