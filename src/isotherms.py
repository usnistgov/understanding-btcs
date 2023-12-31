"""! @brief defines the Langmuir isotherm class
"""

from numpy import sqrt


class Langmuir:
    """! Class for Langmuir isotherm

    @param k The value of \f$\kappa\f$ chosen
    @param name The name of the isotherm

    """
    def __init__(self, k):
        """! Initialize the class
        
        @param k The value of \f$\kappa\f$ chosen
        @param name The name of the isotherm

        @return An instance of the class with specified k
        """
        self.k = k
        self.name = "Langmuir"

    def F(self, c):
        """!
        @return \f$L(c; k) = (1 + \kappa)c/(1 + \kappa c)\f$
        """
        return (1. + self.k)*c/(1 + self.k*c)

    def F_prime(self, c):
        """!
        @return \f$L^\prime(c; k) = (1 + \kappa)/(1 + \kappa c)^2\f$
        """
        return (1. + self.k)/(1. + self.k*c)/(1. + self.k*c)

    def s(self, c):
        """! 
        @return speed at given concentration in x,t coordinate system
        """
        var = (1 + c*self.k)
        return var*var/(1 + self.k + var * var)

    def get_c_rarefaction(self, x, t):
        """!
        @return concentration \f$c(x,t)\f$ associated with rarefaction wave
        """
        if x <= self.s(1.)*t:
            return 1.

        if x >= self.s(0.)*t:
            return 0.

        return 1./self.k*(sqrt((1 + self.k)* x / (t - x)) - 1.)

    def shock_lambda(self):
        """! speed of shock wave
        @return float representing speed of shock wave"""
        return 0.5
    
    def get_c_shock(self, x, t):
        """!
        @return concentration \f$c(x,t)\f$ associated with rarefaction wave
        """
        lamda = self.shock_lambda()
        if x <= lamda*t:
            return 1.

        return 0.

    def is_rarefaction_wave(self):
        """!
        @return whether rarefaction wave or not
        """
        return -1. <= self.k < 0.

    def is_shock_wave(self):
        """!
        @return whether shock wave or not
        """
        return self.k > 0.
    
    def get_c_local_equilibrium(self, x, t):
        """!
        @return \f$c(x,t)\f$ associated with local equilibrium
        """
        if self.is_shock_wave():
            return self.get_c_shock(x, t)
        elif self.is_rarefaction_wave():
            return self.get_c_rarefaction(x, t)
        
        raise Exception("Langmuir isotherm not defined for kappa = %5.2f" % self.k)

