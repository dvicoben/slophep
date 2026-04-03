from flavio import default_parameters

def get_default_params() -> dict:
    par = default_parameters.get_central_all()
    # Excited D states, from hammer
    extrapar = {
        "m_D10"   : 2.421 ,
        "m_D1+"   : 2.423 ,
        "m_D1*0"  : 2.427 ,
        "m_D1*+"  : 2.427 ,
        "m_D0*0"  : 2.300 ,
        "m_D0*+"  : 2.349 ,
        "m_D2*0"  : 2.461 ,
        "m_D2*+"  : 2.465 ,
        "m_Ds0*+" : 2.3178,
        "m_Ds1*+" : 2.4595,
        "m_Ds1+"  : 2.5351,
        "m_Ds2*+" : 2.5691
    }
    par.update(extrapar)
    return par