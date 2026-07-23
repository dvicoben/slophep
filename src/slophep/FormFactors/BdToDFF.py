import slophep.FormFactors.FormFactorsBToP as FFBToP

class BSZ(FFBToP.FFBToP_BSZ):
    _name = "FFBdToD@BSZ"
    def __init__(self):
        super().__init__("B0", "D+")


class BLPR(FFBToP.FFBToP_BLPR):
    _name = "FFBdToD@BLPR"
    def __init__(self):
        super().__init__("B0", "D+")


class BLPRXP(FFBToP.FFBToP_BLPRXP):
    _name = "FFBdToD@BLPRXP"
    def __init__(self):
        super().__init__("B0", "D+")


class BGL(FFBToP.FFBToP_BGL):
    _name = "FFBdToD@BGL"
    def __init__(self):
        super().__init__("B0", "D+")


class CLN(FFBToP.FFBToP_CLN):
    _name = "FFBdToD@CLN"
    def __init__(self):
        super().__init__("B0", "D+")


# class BGL_Hammer(FFBToP.FFBToP_BGL):
# class CLN_Hammer(FFBToP.FFBToP_CLN):
# class BLPR_Hammer(FFBToP.FFBToP_BLPR):
