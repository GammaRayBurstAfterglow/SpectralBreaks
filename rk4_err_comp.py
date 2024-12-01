from math import exp, sqrt, sin, cos, log, pi
import numpy as np
import matplotlib.pyplot as plt

# Define the derivative part of the ODE.  Since this is an ODE,
#   we might need both x and y -- so pass them both regardless
#   of whether or not they're used in the computation.
def f(x,y):
    return 0.5*sin(x)-cos(pi*x)


def euler1(f,a,b,N,y0):
    """ Implementation of Euler's method for solving 1D 1st-order ODEs.
        Input arguments:
        (1) f: function to compute derivative at given x, y
        (2) a: starting point of interval over which we want to solve
                the ODE
        (3) b: ending point of the interval
        (4) N: number of steps to take
        (5) y0: intial value of y, i.e. y(a)
        Returns:
        (1) y: value of y after N steps of Euler method, i.e. y(b)
    """
    h = (b-a)/float(N)
    xs = a + np.arange(N)*h
    y = y0

    for x in xs[:-1]:
        y += h*f(x,y)

    return y


def heun1(f,a,b,N,y0):
    """ Implementation of Heun's method for solving 1D 1st-order ODEs.
        Input arguments:
        (1) f: function to compute derivative at given x, y
        (2) a: starting point of interval over which we want to solve
                the ODE
        (3) b: ending point of the interval
        (4) N: number of steps to take
        (5) y0: intial value of y, i.e. y(a)
        Returns:
        (1) y: value of y after N steps of Heun's method, i.e. y(b)
    """
    h = (b-a)/float(N)
    xs = a + np.arange(N)*h
    y = y0

    for x in xs:
        y += 0.5*h*f(x,y) + 0.5*h*f( x+h, y+h*f(x,y) )

    return y


def rk4_1(f,a,b,N,y0):
    """ Implementation of 4th-order Runge-Kutta method for solving 1D 1st-
        order ODEs.
        Input arguments:
        (1) f: function to compute derivative at given x, y
        (2) a: starting point of interval over which we want to solve
                the ODE
        (3) b: ending point of the interval
        (4) N: number of steps to take
        (5) y0: intial value of y, i.e. y(a)
        Returns:
        (1) y: value of y after N steps of RK4 method, i.e. y(b)
    """
    h = (b-a)/float(N)
    xs = a + np.arange(N)*h
    y = y0

    for x in xs:
        k0 = h*f(x,y)
        k1 = h*f(x+0.5*h, y+0.5*k0)
        k2 = h*f(x+0.5*h, y+0.5*k1)
        k3 = h*f(x+h, y+k2)
        y += (k0 + 2*k1 + 2*k2 + k3)/6.0

    return y


# By putting our main function inside this if statement, we can safely import
#   the module from other scripts without having this code execute every time
if __name__ == '__main__':

    # Set initial quantities
    a = 0.0
    b = 20.0
    y0 = 0.0

    # Create the arrays to hold x and y values
    ys_euler = []
    ys_heun  = []
    ys_rk4   = []

    # Fill out the array of N values to use
    Ns = [int(10**(1.0+0.25*i)) for i in range(13)]

    # Run through the Ns and fill out the ys array for each method
    for N in Ns:
        ys_euler.append( euler1(f,a,b,N,y0) )
        ys_heun.append( heun1(f,a,b,N,y0) )
        ys_rk4.append( rk4_1(f,a,b,N,y0) )

    # Compute the analytical solution
    ans = y0 - 0.5*cos(b) - sin(pi*b)/pi + 0.5*cos(a) + sin(pi*a)/pi
    # Compute the errors for each method
    euler_errs = [abs(y-ans)/abs(y) for y in ys_euler]
    heun_errs = [abs(y-ans)/abs(y) for y in ys_heun]
    rk4_errs = [abs(y-ans)/abs(y) for y in ys_rk4]

    # Compute the error scaling.  Since err = k*h^sigma, err = k*N^-sigma,
    #   which means err2/err1 = (N2/N1)^-sigma.  This can be solved for
    #   sigma, which is what appears below
    sig_euler = [-log(euler_errs[i+1]/euler_errs[i])/log(Ns[i+1]/Ns[i]) for i in range(len(Ns)-1)]
    sig_heun = [-log(heun_errs[i+1]/heun_errs[i])/log(Ns[i+1]/Ns[i]) for i in range(len(Ns)-1)]
    sig_rk4 = [-log(rk4_errs[i+1]/rk4_errs[i])/log(Ns[i+1]/Ns[i]) for i in range(len(Ns)-1)]


    # Plot the results
    # Note the use of multiple subplots now, and the different ax variables
    #   for each
    fig, (ax1,ax2) = plt.subplots(nrows=2, gridspec_kw={'height_ratios':[2,1]}, sharex=True)
    plt.subplots_adjust(wspace=0.4,hspace=0.0)
    ax1.plot( Ns, euler_errs, 'ko-', linewidth=2, markersize=3, label="Euler" )
    ax1.plot( Ns, heun_errs, 'ro-', linewidth=2, markersize=3, label="Heun" )
    ax1.plot( Ns, rk4_errs, 'bo-', linewidth=2, markersize=3, label="RK4" )
    ax2.plot( Ns[1:], sig_euler, 'ko-', linewidth=2, markersize=3, label="Euler")
    ax2.plot( Ns[1:], sig_heun, 'ro-', linewidth=2, markersize=3, label="Heun")
    ax2.plot( Ns[1:], sig_rk4, 'bo-', linewidth=2, markersize=3, label="RK4")
    plt.xscale('log')
    ax1.set_yscale('log')
    ax2.set_yscale('linear')
    ax2.set_ylim(-1.0,5.5)
    ax1.tick_params(right=True)
    ax1.tick_params(labelright=True)
    ax2.tick_params(right=True)
    
    ax1.legend(loc='lower left') # Make sure to include the "label" argument for your curves
    ax2.set_xlabel("N")
    ax1.set_ylabel("relative error")
    ax2.set_ylabel("error scaling")
    
    plt.show()
