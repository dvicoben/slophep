import slophep.FormFactors.FormFactorsBToDstst as FFBToDstst

class ISGW2(FFBToDstst.FFBToD1st_ISGW2):
    _name = "BsToDs1st@ISGW2"
    def __init__(self):
        super().__init__("Bs", "Ds1*+")

    def define_userparams(self):
        ffpar = {
            "msb"     : 5.2                ,
            "msd"     : 0.55               ,
            "bb2"     : 0.54*0.54          ,
            "mbb"     : 5.38               ,
            "msq"     : 1.82               ,
            "bx2"     : 0.41*0.41          ,
            "mbx"     : (3.0*2.54+2.46)/4.0,
            "mqm"     : 0.1                ,
            "nfp"     : 3.0                ,
            "SmearQ2" : True
        }
        return ffpar


class LLSW(FFBToDstst.FFBToD1st_LLSW):
    _name = "BsToDs1st@LLSW"
    def __init__(self):
        super().__init__("Bs", "Ds1*+")


class BLR(FFBToDstst.FFBToD1st_BLR):
    _name = "BsToDs1st@BLR"
    def __init__(self):
        super().__init__("Bs", "Ds1*+")