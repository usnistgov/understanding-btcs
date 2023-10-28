import numpy as np


class Linear:
    def __init__(self):
        self.name = "linear"
    
    def f(self, c):
        return c
    
    def df_dc(self, c):
        return np.ones_like(c, dtype=float)


def f_langmuir_dimensional(c, k, H):
    return H*c/(1 + k*c)


def df_langmuir_dimensional_dH(c, k, H):
    return c/(1 + k*c)


def df_langmuir_dimensional_dk(c, k, H):
    return -H*c*c/(1 + k*c)/(1 + k*c)


def f_langmuir_dimensional_prime(c, k, H):
    return H/(1 + k*c)/(1 + k*c)


def f_langmuir(c, k):
    return f_langmuir_dimensional(c, k, 1 + k)


def df_langmuir_dk(c, k):
    return c*(1 - c)/(1 + k*c)/(1 + k*c)


def f_langmuir_prime(c, k):
    return (k + 1)/(1 + k*c)/(1 + k*c)


class Langmuir:
    def __init__(self, k):
        self.k = k
        self.name = "Langmuir"

    def f(self, c):
        return f_langmuir(c, self.k)

    def df_dc(self, c):
        return f_langmuir_prime(c, self.k)

    def s(self, c):
        var = (1 + c*self.k)
        return var*var/(1 + self.k + var * var)

    def get_c_rarefaction(self, x, t):
        if x <= self.s(1.)*t:
            return 1.

        if x >= self.s(0.)*t:
            return 0.

        return 1./self.k*(np.sqrt((1 + self.k)* x / (t - x)) - 1.)
    
    def get_c_shock(self, x, t):
        lamda = self.shock_lambda()
        if x <= lamda*t:
            return 1.

        return 0.

    def is_rarefaction_wave(self):
        return -1. <= self.k < 0.

    def is_shock_wave(self):
        return self.k > 0.

    def shock_lambda(self):
        return 1/2
    
    def get_c_local_equilibrium(self, x, t):
        if self.is_shock_wave():
            return self.get_c_shock(x, t)
        
        return self.get_c_rarefaction(x, t)

