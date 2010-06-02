

to_array = lambda S: [S['Rho'], S['Pre'],
                      S['v'][0], S['v'][1], S['v'][2],
                      S['B'][0], S['B'][1], S['B'][2]]



class CylindricalProblem:

    def __init__(self, I={ }, O={ }, gamma=1.4):

        self._setup()

        self.I_state.update(I)
        self.O_state.update(O)
        self.adiabatic_gamma = gamma


    def initial_model(self, P):

        assert len(P.shape) is 3
        from numpy import sqrt, array, zeros_like, linspace, where, newaxis

        Nx, Ny = P.shape[0:2]

        X = linspace(-1,1,Nx)[:,newaxis]
        Y = linspace(-1,1,Ny)[newaxis,:]

        PI = zeros_like(P)
        PO = zeros_like(P)

        for i in range(8):
            PI[:,:,i] = to_array(self.I_state)[i]
            PO[:,:,i] = to_array(self.O_state)[i]

        for i in range(8):
            P[:,:,i] = where(sqrt(X**2 + Y**2) < 0.16, PI[:,:,i], PO[:,:,i])




class RMHDCylindricalA(CylindricalProblem):

    def __init__(self, I={ }, O={ }, pre=1.0, gamma=1.4):

        self.pre = pre
        CylindricalProblem.__init__(self, I=I, O=O, gamma=gamma)

    def _setup(self):

        self.I_state = { 'Rho':1.0, 'Pre':self.pre, 'v': [0,0,0], 'B': [4,0,0] }
        self.O_state = { 'Rho':1.0, 'Pre':    0.01, 'v': [0,0,0], 'B': [4,0,0] }



class ShockTubeProblem:

    def __init__(self, L={ }, R={ }, gamma=1.4, orientation='x'):

        self._setup()

        self.L_state.update(L)
        self.R_state.update(R)
        self.adiabatic_gamma = gamma
        self.orientation = orientation

    def initial_model(self, P):

        if len(P.shape) is 2:

            Nx = P.shape[0]

            for i in range(8):
                P[:Nx/2,i] = to_array(self.L_state)[i]
                P[Nx/2:,i] = to_array(self.R_state)[i]

        elif len(P.shape) is 3:

            """
            Offer the option to transpose shocktube direction. This may be
            useful in checking the dimensional equivalence of 2d algorithms.
            """
            Nx, Ny = P.shape[0:2]

            if self.orientation == 'x':

                for i in range(8):
                    P[:Nx/2 ,:,i] = to_array(self.L_state)[i]
                    P[ Nx/2:,:,i] = to_array(self.R_state)[i]

            elif self.orientation == 'y':

                for S in [self.L_state, self.R_state]:
                    S['v'][0], S['v'][1] = S['v'][1], S['v'][0]
                    S['B'][0], S['B'][1] = S['B'][1], S['B'][0]

                for i in range(8):
                    P[:,:Ny/2 ,i] = to_array(self.L_state)[i]
                    P[:, Ny/2:,i] = to_array(self.R_state)[i]


    def get_states(self):

        return to_array(self.L_state), to_array(self.R_state)



class SRShockTube1(ShockTubeProblem):

    """
    This is Problem 1 of Marti & Muller's article, except that the
    pressure of the left state has been set apart from zero.

    See: http://relativity.livingreviews.org/Articles/lrr-2003-7/
    """

    def __init__(self, **kwargs):

        ShockTubeProblem.__init__(self, **kwargs)

    def _setup(self):

        self.L_state = { 'Rho':10.0, 'Pre':13.33, 'v': [0,0,0], 'B': [0,0,0] }
        self.R_state = { 'Rho': 1.0, 'Pre': 0.01, 'v': [0,0,0], 'B': [0,0,0] }


class SRShockTube2(ShockTubeProblem):

    """
    This is Problem 2 of Marti & Muller's article. It is very
    challenging and still fails many modern SR codes.

    See: http://relativity.livingreviews.org/Articles/lrr-2003-7/
    """

    def __init__(self, **kwargs):

        ShockTubeProblem.__init__(self, **kwargs)

    def _setup(self):

        self.L_state = { 'Rho': 1.0, 'Pre':1000.00, 'v': [0,0,0], 'B': [0,0,0] }
        self.R_state = { 'Rho': 1.0, 'Pre':   0.01, 'v': [0,0,0], 'B': [0,0,0] }


class RMHDShockTube1(ShockTubeProblem):

    def __init__(self, **kwargs):

        ShockTubeProblem.__init__(self, **kwargs)

    def _setup(self):

        self.L_state = { 'Rho': 1.000, 'Pre':1.0, 'v': [0,0,0], 'B': [0.5, 1.0, 0.0] }
        self.R_state = { 'Rho': 0.125, 'Pre':0.1, 'v': [0,0,0], 'B': [0.5,-1.0, 0.0] }


class RMHDShockTube2(ShockTubeProblem):

    def __init__(self, **kwargs):

        ShockTubeProblem.__init__(self, **kwargs)

    def _setup(self):

        self.L_state = { 'Rho': 1.08, 'Pre': 0.95, 'v': [ 0.40, 0.3, 0.2], 'B': [2.0, 0.3, 0.3] }
        self.R_state = { 'Rho': 0.95, 'Pre': 1.00, 'v': [-0.45,-0.2, 0.2], 'B': [2.0,-0.7, 0.5] }


class RMHDShockTube3(ShockTubeProblem):

    def __init__(self, **kwargs):

        ShockTubeProblem.__init__(self, **kwargs)

    def _setup(self):

        self.L_state = { 'Rho': 1.00, 'Pre': 0.1, 'v': [ 0.999, 0.0, 0.0], 'B': [10.0, 7.0, 7.0] }
        self.R_state = { 'Rho': 1.00, 'Pre': 0.1, 'v': [-0.999, 0.0, 0.0], 'B': [10.0,-7.0,-7.0] }


class RMHDShockTube4(ShockTubeProblem):

    def __init__(self, **kwargs):

        ShockTubeProblem.__init__(self, **kwargs)

    def _setup(self):

        self.L_state = { 'Rho': 1.0, 'Pre': 5.0, 'v': [0.0, 0.3, 0.4], 'B': [1.0, 6.0, 2.0] }
        self.R_state = { 'Rho': 0.9, 'Pre': 5.3, 'v': [0.0, 0.0, 0.0], 'B': [1.0, 5.0, 2.0] }


class RMHDContactWave(ShockTubeProblem):

    def __init__(self, **kwargs):

        ShockTubeProblem.__init__(self, **kwargs)

    def _setup(self):

        self.L_state = { 'Rho': 1.0, 'Pre': 1.0, 'v': [0.0, 0.7, 0.2], 'B': [5.0, 1.0, 0.5] }
        self.R_state = { 'Rho': 0.1, 'Pre': 1.0, 'v': [0.0, 0.7, 0.2], 'B': [5.0, 1.0, 0.5] }


class RMHDRotationalWave(ShockTubeProblem):

    def __init__(self, **kwargs):

        ShockTubeProblem.__init__(self, **kwargs)

    def _setup(self):

        vR = [0.400000,-0.300000, 0.500000]
        vL = [0.377347,-0.482389, 0.424190]

        self.L_state = { 'Rho': 1.0, 'Pre': 1.0, 'v': vL, 'B': [2.4, 1.0,-1.600000] }
        self.R_state = { 'Rho': 1.0, 'Pre': 1.0, 'v': vR, 'B': [2.4,-0.1,-2.178213] }
