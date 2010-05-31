
def shocktube(P, extent=(-1,1), **kwargs):

    from pylab import sqrt, linspace, subplot, plot, text, xlabel, figure, show
    from pylab import subplots_adjust, setp, gca, LinearLocator, rcParams, legend

    rho, pre   = P[:,0], P[:,1]
    vx, vy, vz = P[:,2], P[:,3], P[:,4]
    Bx, By, Bz = P[:,5], P[:,6], P[:,7]

    plot_args = { }
    plot_args['marker'] = kwargs.get('marker', 'o')
    plot_args['c'     ] = kwargs.get('c'     , 'k')
    plot_args['mfc'   ] = kwargs.get('mfc'   , 'None')

    plot_args.update(kwargs)
    rcParams.update({'axes.labelsize':16, 'ytick.major.pad':8})

    x = linspace(extent[0],extent[1],P.shape[0])
    g = 1 / sqrt(1-(vx**2+vy**2+vz**2))

    ax = subplot(2,3,1)
    plot(x,rho, **plot_args)
    text(0.9,0.9, r"$\rho$", transform = ax.transAxes, fontsize=20)
    setp(ax.get_xticklabels(), visible=False)
    if 'label' in plot_args: legend(loc='lower left')

    ax = subplot(2,3,2)
    plot(x,pre, **plot_args)
    text(0.9,0.9, r"$P$", transform = ax.transAxes, fontsize=20)
    setp(ax.get_xticklabels(), visible=False)

    ax = subplot(2,3,3)
    plot(x, g, **plot_args)
    text(0.9,0.9, r"$\gamma$", transform = ax.transAxes, fontsize=20)
    setp(ax.get_xticklabels(), visible=False)

    ax = subplot(2,3,4)
    plot(x, vx, **plot_args)
    text(0.9,0.9, r"$v_x$", transform = ax.transAxes, fontsize=20)
    xlabel(r"$x$")

    ax = subplot(2,3,5)
    plot(x, vy, **plot_args)
    text(0.9,0.9, r"$v_y$", transform = ax.transAxes, fontsize=20)
    xlabel(r"$x$")

    ax = subplot(2,3,6)
    plot(x, By, **plot_args)
    text(0.9,0.9, r"$B_y$", transform = ax.transAxes, fontsize=20)
    xlabel(r"$x$")

    subplots_adjust(hspace=0.05)

