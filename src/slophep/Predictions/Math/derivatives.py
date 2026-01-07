import numpy as np
import sys

_DER_EPS = 1.e-4

def derivative_left(f, val: float, eps: float = _DER_EPS) -> float:
    return (4.*f(val-eps)-f(val-2.*eps)-3.*f(val))/(-2.*eps)

def derivative_right(f, val: float, eps: float = _DER_EPS) -> float:
    return (4.*f(val+eps)-f(val+2.*eps)-3.*f(val))/(2.*eps)

def derivative_centre(f, val: float, eps: float = _DER_EPS) -> float:
    return (f(val + eps) - f(val - eps))/(2.*eps)

def derivative(f, val: float, eps: float = _DER_EPS, 
               bound_lo: float = -sys.float_info.max, 
               bound_hi: float = sys.float_info.max):
    if val <= bound_lo:
        return derivative_right(f, bound_lo, eps)
    if val >= bound_hi:
        return derivative_right(f, bound_hi, eps)
    
    return derivative_centre(f, val, eps)