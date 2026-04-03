import slophep.Predictions.FormFactorsBToDstst as FFBToDstst

class ISGW2(FFBToDstst.ISGW2_BToD2st):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("Bs", "Ds2*+", par, scale, *ffargs)
        self._name = "BsToDs2st_ISGW2"
        internalparams = {
            "msb" : 5.2                    ,
            "msd" : 0.55                   ,
            "bb2" : 0.54*0.54              ,
            "mbb" : 5.38                   ,
            "msq" : 1.82                   ,
            "bx2" : 0.41*0.41              ,
            "mbx" : (5.0*2.61+3.0*2.54)/8.0,
            "mqm" : 0.1                    ,
            "nfp" : 3.0
        }
        self._internalparams.update(internalparams)

class BLR(FFBToDstst.BLR_BToD2st):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("Bs", "Ds2*+", par, scale, *ffargs)
        self._name = "BsToDs2st_BLR"


class LLSW(FFBToDstst.LLSW_BToD2st):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("Bs", "Ds2*+", par, scale, *ffargs)
        self._name = "BsToDs2st_LLSW"