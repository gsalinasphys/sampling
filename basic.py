import numpy as np


def mag(vs: np.ndarray) -> np.ndarray:
    return np.sqrt(np.sum(vs**2, axis = 1))

def cases(cond: np.ndarray, inside: np.ndarray, outside: np.ndarray):
    return inside*np.heaviside(-cond, 0.) + outside*np.heaviside(cond, 1.)