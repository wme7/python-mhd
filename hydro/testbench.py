

to_array = lambda S: [S['Rho'], S['Pre'],
                      S['v'][0], S['v'][1], S['v'][2],
                      S['B'][0], S['B'][1], S['B'][2]]


class TestProblemBase:

    def __init__(self):
        raise NotImplementedError

    def initial_model(self, domain, Ng, Nq):
        from numpy import zeros, linspace, vectorize, zeros_like, newaxis as n
        x0, x1, dx, N = domain.x0, domain.x1, domain.dx, domain.N
        vfunc = vectorize(self.prim_at_point)

        if len(domain.N) is 1:
            X = linspace(x0[0]-Ng*dx[0], x1[0]+Ng*dx[0], N[0]+2*Ng)
            Y = zeros_like(X)
            Z = zeros_like(Y)

        if len(domain.N) is 2:
            X = linspace(x0[0]-Ng*dx[0], x1[0]+Ng*dx[0], N[0]+2*Ng)[:,n]
            Y = linspace(x0[1]-Ng*dx[1], x1[1]+Ng*dx[1], N[1]+2*Ng)[n,:]
            Z = zeros_like(Y)

        if len(domain.N) is 3:
            X = linspace(x0[0]-Ng*dx[0], x1[0]+Ng*dx[0], N[0]+2*Ng)[:,n,n]
            Y = linspace(x0[1]-Ng*dx[1], x1[1]+Ng*dx[1], N[1]+2*Ng)[n,:,n]
            Z = linspace(x0[2]-Ng*dx[2], x1[2]+Ng*dx[2], N[2]+2*Ng)[n,n,:]

        P = zeros(tuple([n+2*Ng for n in domain.N])+(Nq,))
        for i in range(Nq):
            P[...,i] = vfunc(X,Y,Z,i)

        return P


class CylindricalProblem(TestProblemBase):

    def __init__(self, I={ }, O={ }, gamma=1.4):
        self._setup_()
        self.I_state.update(I)
        self.O_state.update(O)
        self.adiabatic_gamma = gamma

    def prim_at_point(self, x,y,z,i):
        I = self.I_state
        O = self.O_state
        r = (x**2 + y**2 + z**2)**0.5
        return to_array(I)[i] if r<0.16 else to_array(O)[i]


class RMHDCylindricalA(CylindricalProblem):

    def __init__(self, pre=1.0, **kwargs):
        self.pre = pre
        CylindricalProblem.__init__(self, **kwargs)

    def _setup_(self):
        self.I_state = { 'Rho':1.0, 'Pre':self.pre, 'v': [0,0,0], 'B': [4,0,0] }
        self.O_state = { 'Rho':1.0, 'Pre':    0.01, 'v': [0,0,0], 'B': [4,0,0] }


class QuadrantProblem(TestProblemBase):

    def __init__(self, NE={ }, NW={ }, SE={ }, SW={ }, gamma=1.4):
        self._setup_()
        self.NE_state.update(NE)
        self.NW_state.update(NW)
        self.SE_state.update(SE)
        self.SW_state.update(SW)
        self.adiabatic_gamma = gamma

    def prim_at_point(self, x,y,z,i):
        if x<=0.0 and y<=0.0:
            return to_array(self.SW_state)[i]
        if x<=0.0 and y >0.0:
            return to_array(self.NW_state)[i]
        if x >0.0 and y<=0.0:
            return to_array(self.SE_state)[i]
        if x >0.0 and y >0.0:
            return to_array(self.NE_state)[i]


class SRQuadrantA(QuadrantProblem):

    def __init__(self, **kwargs):
        QuadrantProblem.__init__(self, **kwargs)

    def _setup_(self):
        self.SW_state = {'Rho':1.0, 'Pre':1.0, 'v': [0.00, 0.00, 0.00], 'B': [0,0,0]}
        self.NW_state = {'Rho':1.0, 'Pre':1.0, 'v': [0.99, 0.00, 0.00], 'B': [0,0,0]}
        self.SE_state = {'Rho':1.0, 'Pre':1.0, 'v': [0.00, 0.99, 0.00], 'B': [0,0,0]}
        self.NE_state = {'Rho':1.0, 'Pre':1.0, 'v': [0.00, 0.00, 0.00], 'B': [0,0,0]}


class SRQuadrantB(QuadrantProblem):

    def __init__(self, **kwargs):
        QuadrantProblem.__init__(self, **kwargs)

    def _setup_(self):
        self.SW_state = {'Rho':1.0, 'Pre':1.0, 'v': [0.00, 0.00, 0.00], 'B': [0,0,0]}
        self.NW_state = {'Rho':1.0, 'Pre':1.0, 'v': [0.80, 0.00, 0.00], 'B': [0,0,0]}
        self.SE_state = {'Rho':1.0, 'Pre':1.0, 'v': [0.00, 0.80, 0.00], 'B': [0,0,0]}
        self.NE_state = {'Rho':1.0, 'Pre':1.0, 'v': [0.00, 0.00, 0.00], 'B': [0,0,0]}


class ShockTubeProblem(TestProblemBase):

    def __init__(self, L={ }, R={ }, gamma=1.4, orientation='x'):
        self._setup_()
        self.L_state.update(L)
        self.R_state.update(R)
        self.adiabatic_gamma = gamma
        self.orientation = orientation

    def prim_at_point(self, x,y,z,i):
        L = self.L_state
        R = self.R_state
        return to_array(L)[i] if x<0.5 else to_array(R)[i]


class SRShockTube1(ShockTubeProblem):
    """
    This is Problem 1 of Marti & Muller's article, except that the
    pressure of the left state has been set apart from zero.

    See: http://relativity.livingreviews.org/Articles/lrr-2003-7/
    """
    def __init__(self, **kwargs):
        ShockTubeProblem.__init__(self, **kwargs)

    def _setup_(self):
        self.L_state = {'Rho':10.0, 'Pre':13.33, 'v': [0,0,0], 'B': [0,0,0]}
        self.R_state = {'Rho': 1.0, 'Pre': 0.01, 'v': [0,0,0], 'B': [0,0,0]}


class SRShockTube2(ShockTubeProblem):
    """
    This is Problem 2 of Marti & Muller's article. It is very challenging
    and still fails many modern SR codes.

    See: http://relativity.livingreviews.org/Articles/lrr-2003-7/
    """
    def __init__(self, **kwargs):
        ShockTubeProblem.__init__(self, **kwargs)

    def _setup_(self):
        self.L_state = {'Rho': 1.0, 'Pre':1000.00, 'v': [0,0,0], 'B': [0,0,0]}
        self.R_state = {'Rho': 1.0, 'Pre':   0.01, 'v': [0,0,0], 'B': [0,0,0]}


class SRShockTube3(ShockTubeProblem):
    """
    This is the transverse velocity problem presented in section 6-1
    of Zhang & MacFadyen (2005). It is identical to SRShockTube2, except
    for the presence of transverse velocity on both sides of the domain.
    """
    def __init__(self, **kwargs):
        ShockTubeProblem.__init__(self, **kwargs)

    def _setup_(self):
        self.L_state = {'Rho': 1.0, 'Pre':1000.00, 'v': [0,0.9,0], 'B': [0,0,0]}
        self.R_state = {'Rho': 1.0, 'Pre':   0.01, 'v': [0,0.9,0], 'B': [0,0,0]}


class RMHDShockTube1(ShockTubeProblem):

    def __init__(self, **kwargs):
        ShockTubeProblem.__init__(self, **kwargs)

    def _setup_(self):
        self.L_state = {'Rho': 1.000, 'Pre':1.0, 'v': [0,0,0], 'B': [0.5, 1.0, 0.0]}
        self.R_state = {'Rho': 0.125, 'Pre':0.1, 'v': [0,0,0], 'B': [0.5,-1.0, 0.0]}


class RMHDShockTube2(ShockTubeProblem):

    def __init__(self, **kwargs):
        ShockTubeProblem.__init__(self, **kwargs)

    def _setup_(self):
        self.L_state = {'Rho': 1.08, 'Pre': 0.95, 'v': [ 0.40, 0.3, 0.2], 'B': [2.0, 0.3, 0.3]}
        self.R_state = {'Rho': 0.95, 'Pre': 1.00, 'v': [-0.45,-0.2, 0.2], 'B': [2.0,-0.7, 0.5]}


class RMHDShockTube3(ShockTubeProblem):

    def __init__(self, **kwargs):
        ShockTubeProblem.__init__(self, **kwargs)

    def _setup_(self):
        self.L_state = {'Rho': 1.00, 'Pre': 0.1, 'v': [ 0.999, 0.0, 0.0], 'B': [10.0, 7.0, 7.0]}
        self.R_state = {'Rho': 1.00, 'Pre': 0.1, 'v': [-0.999, 0.0, 0.0], 'B': [10.0,-7.0,-7.0]}


class RMHDShockTube4(ShockTubeProblem):

    def __init__(self, **kwargs):
        ShockTubeProblem.__init__(self, **kwargs)

    def _setup_(self):
        self.L_state = {'Rho': 1.0, 'Pre': 5.0, 'v': [0.0, 0.3, 0.4], 'B': [1.0, 6.0, 2.0]}
        self.R_state = {'Rho': 0.9, 'Pre': 5.3, 'v': [0.0, 0.0, 0.0], 'B': [1.0, 5.0, 2.0]}


class RMHDContactWave(ShockTubeProblem):

    def __init__(self, **kwargs):
        ShockTubeProblem.__init__(self, **kwargs)

    def _setup_(self):
        self.L_state = {'Rho': 1.0, 'Pre': 1.0, 'v': [0.0, 0.7, 0.2], 'B': [5.0, 1.0, 0.5]}
        self.R_state = {'Rho': 0.1, 'Pre': 1.0, 'v': [0.0, 0.7, 0.2], 'B': [5.0, 1.0, 0.5]}


class RMHDRotationalWave(ShockTubeProblem):

    def __init__(self, **kwargs):
        ShockTubeProblem.__init__(self, **kwargs)

    def _setup_(self):
        vR = [0.400000,-0.300000, 0.500000]
        vL = [0.377347,-0.482389, 0.424190]

        self.L_state = {'Rho': 1.0, 'Pre': 1.0, 'v': vL, 'B': [2.4, 1.0,-1.600000]}
        self.R_state = {'Rho': 1.0, 'Pre': 1.0, 'v': vR, 'B': [2.4,-0.1,-2.178213]}


class KelvinHelmholtzProblem(TestProblemBase):

    def __init__(self, I={ }, O={ }, gamma=1.4):
        self._setup_()
        self.I_state.update(I)
        self.O_state.update(O)
        self.adiabatic_gamma = gamma

    def prim_at_point(self, x,y,z,i):
        inside = -0.25 < y and y < 0.25
        I = self.I_state
        O = self.O_state
        p = to_array(I)[i] if inside else to_array(O)[i]
        return p + self.perturbation(x,y,z,i)


class AthenaKelvinHelmholtz(KelvinHelmholtzProblem):

    def __init__(self, **kwargs):
        KelvinHelmholtzProblem.__init__(self, **kwargs)

    def _setup_(self):        
        self.I_state = {'Rho': 2.0, 'Pre': 2.5, 'v': [-0.5,0.0,0.0], 'B': [0.5,0.0,0.0]}
        self.O_state = {'Rho': 1.0, 'Pre': 2.5, 'v': [ 0.5,0.0,0.0], 'B': [0.5,0.0,0.0]}

    def perturbation(self, x,y,z,i):
        from math import sin, pi
        if i == 2:
            return 0.01 * sin(2*pi*x)
        elif i == 3:
            return 0.01 * sin(2*pi*x)
        else:
            return 0.0
