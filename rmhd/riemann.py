
import sr_riemann

def exact_sr(problem, tfinal=0.4):

    from os import system
    from numpy import zeros

    if problem.L_state['B'] != [0,0,0] or problem.R_state['B'] != [0,0,0]:
        print """
              Warning! exact riemann solver for SR only, B = (%s / %s) being ignored
              Your solution will be incorrect.
              """ % (problem.L_state['B'], problem.R_state['B'])

    if problem.L_state['v'][1:3] != [0,0] or problem.R_state['v'][1:3] != [0,0]:
        print """
              Warning! exact riemann solver for normal velocity components only.
              Your solution will be incorrect. Use exact_sr_vt for riemann problems
              with tangential velocities.
              """

    sr_riemann.riemann(tfinal, problem.adiabatic_gamma,
                       problem.L_state['Pre'],
                       problem.L_state['Rho'],
                       problem.L_state['v'][0],
                       problem.R_state['Pre'],
                       problem.R_state['Rho'],
                       problem.R_state['v'][0])

    soln_lines = open('solution.dat').readlines()
    system('rm -f solution.dat')

    x,pre,rho,vel,u = [[float(m.split()[i]) for m in soln_lines] for i in range(5)]
    Nx = len(x)

    P = zeros((Nx,8))
    P[:,0] = rho
    P[:,1] = pre
    P[:,2] = vel

    return P


def exact_sr_vt(problem, tfinal=0.4):

    from os import system
    from numpy import zeros

    if problem.L_state['B'] != [0,0,0] or problem.R_state['B'] != [0,0,0]:
        print """
              Warning! exact riemann solver for SR only, B = (%s / %s) being ignored
              Your solution will be incorrect.
              """ % (problem.L_state['B'], problem.R_state['B'])

    if problem.L_state['v'][2] != 0 or problem.R_state['v'][2] != 0:
        print """
              Warning! exact riemann solver configured to have transverse velocities
              stored in y-component only. Ignoring z velocity value.
              Your solution will be incorrect.
              """

    sr_riemann.riemann_vt(tfinal, problem.adiabatic_gamma,
                          problem.L_state['Pre'],
                          problem.L_state['Rho'],
                          problem.L_state['v'][0],
                          problem.L_state['v'][1],
                          problem.R_state['Pre'],
                          problem.R_state['Rho'],
                          problem.R_state['v'][0],
                          problem.R_state['v'][1])

    soln_lines = open('solution.dat').readlines()
    system('rm -f solution.dat')

    x,pre,rho,velx,vely = [[float(m.split()[i]) for m in soln_lines] for i in range(5)]
    Nx = len(x)

    P = zeros((Nx,8))
    P[:,0] = rho
    P[:,1] = pre
    P[:,2] = velx
    P[:,3] = vely

    return P
