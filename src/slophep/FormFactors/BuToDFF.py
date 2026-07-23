import slophep.FormFactors.FormFactorsBToP as FFBToP

class BSZ(FFBToP.FFBToP_BSZ):
    _name = "FFBuToD@BSZ"
    def __init__(self):
        super().__init__("B+", "D0")


class BLPR(FFBToP.FFBToP_BLPR):
    _name = "FFBuToD@BLPR"
    def __init__(self):
        super().__init__("B+", "D0")


class BLPRXP(FFBToP.FFBToP_BLPRXP):
    _name = "FFBuToD@BLPRXP"
    def __init__(self):
        super().__init__("B+", "D0")


class BGL(FFBToP.FFBToP_BGL):
    _name = "FFBuToD@BGL"
    def __init__(self):
        super().__init__("B+", "D0")


class CLN(FFBToP.FFBToP_CLN):
    _name = "FFBuToD@CLN"
    def __init__(self):
        super().__init__("B+", "D0")


# class BGL_Hammer(FFBToP.FFBToP_BGL):
# class CLN_Hammer(FFBToP.FFBToP_CLN):
# class BLPR_Hammer(FFBToP.FFBToP_BLPR):
