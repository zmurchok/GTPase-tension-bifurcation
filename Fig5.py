import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.integrate import solve_ivp
from scipy import optimize
import jax.numpy as jnp
from jax import jacfwd
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.transforms as mtransforms






def plot_setup():
    plt.rcParams['text.usetex'] = False
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['mathtext.fontset'] = 'cm'
    # plt.rcParams['figure.autolayout'] = 'True'

    sns.set_style('ticks',{'axes.edgecolor': '[0,0,0]',
                           'xtick.direction':'in',
                           'ytick.direction':'in',
                           'ytick.right':'True',
                           'xtick.top':'True',
                           'xtick.color':'k',
                           'ytick.color':'k'
                            })






def add_arrow_to_line2D(
    axes, line, arrow_locs=[0.2, 0.4, 0.6, 0.8],
    arrowstyle='-|>', arrowsize=1, transform=None):
    """
    Add arrows to a matplotlib.lines.Line2D at selected locations.

    Parameters:
    -----------
    axes:
    line: Line2D object as returned by plot command
    arrow_locs: list of locations where to insert arrows, % of total length
    arrowstyle: style of the arrow
    arrowsize: size of the arrow
    transform: a matplotlib transform instance, default to data coordinates

    Returns:
    --------
    arrows: list of arrows
    """
    if not isinstance(line, mlines.Line2D):
        raise ValueError("expected a matplotlib.lines.Line2D object")
    x, y = line.get_xdata(), line.get_ydata()

    arrow_kw = {
        "arrowstyle": arrowstyle,
        "mutation_scale": 10 * arrowsize,
    }

    color = line.get_color()
    use_multicolor_lines = isinstance(color, np.ndarray)
    if use_multicolor_lines:
        raise NotImplementedError("multicolor lines not supported")
    else:
        arrow_kw['color'] = color

    linewidth = line.get_linewidth()
    if isinstance(linewidth, np.ndarray):
        raise NotImplementedError("multiwidth lines not supported")
    else:
        arrow_kw['linewidth'] = linewidth

    if transform is None:
        transform = axes.transData

    arrows = []
    for loc in arrow_locs:
        s = np.cumsum(np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2))
        n = np.searchsorted(s, s[-1] * loc)
        arrow_tail = (x[n], y[n])
        arrow_head = (np.mean(x[n:n + 2]), np.mean(y[n:n + 2]))
        p = mpatches.FancyArrowPatch(
            arrow_tail, arrow_head, transform=transform,
            **arrow_kw)
        axes.add_patch(p)
        arrows.append(p)
    return arrows




def f(t, y, b, beta):
    """
    ODE system

    Parameters:
    -----------
    t: time
    y: dimensions of system
    b: basal activation rate b
    beta: feedback from tension parameter

    Returns:
    --------
    [dGdt, dLdt]: differential equations
    """

    gamma = 1.5
    G_T = 2
    ell0 = 1
    phi = 0.75
    G_h = 0.3
    epsilon = 0.1
    alpha = 10
    n = 4

    G, L = y
    restlength = (ell0 - phi * G**n / (G_h**n + G**n))
    squashing = beta / (1 + np.e**(-alpha * (L - restlength)))

    dGdt = (b + squashing + gamma * G**n / (1+G**n) ) * (G_T - G) - G
    dLdt = -epsilon*(L - restlength)

    return [dGdt, dLdt]



if __name__ == '__main__':
    plot_setup()

    b = 0.14088
    beta_pars = [0.16508040,0.16508060,0.16508100]

    #plot set up
    fig, axes = plt.subplots(1,3, figsize=(6.5,2), sharex=True, sharey=True, constrained_layout=True)

    #loop over two parameters of interest
    for i in range(len(beta_pars)):
        #set parameter dictionary
        p = {'beta':beta_pars[i],'b':b}
        ax = axes[i]

        #set up functions to integrate
        fpar = lambda t,y: f(t,y,**p)
        frhs = lambda y: fpar(0,y)

        #find steady-states using root finding and initial guesses
        node = optimize.root(frhs, [1,0.2])
        saddle = optimize.root(frhs, [0.7,0.25])

        #plot steady-states with a dot
        fill='k'
        nodeeq = ax.scatter(node.x[0],node.x[1],c=fill,zorder=20,edgecolors='k',label='Node')
        fill='w'
        saddleeq = ax.scatter(saddle.x[0],saddle.x[1],c=fill,zorder=20,edgecolors='k',label='Saddle')

        #numerically calculate jacobian
        #automatic numerical differentiation (we could to this by hand or symbolically but I'm lazy)
        fgrad = jacfwd(frhs)

        #get eigenvalues w and eigenvectors v of the saddle point (columns of v are eigenvectors)
        w, v = np.linalg.eig(fgrad(saddle.x))

        #solution on unstable manifold (GREY CURVE)
        T = 2000 #need large time
        tol = 1e-4 #distance from saddle node along eigenvector for initial condition
        y0 = saddle.x - tol*v[:,0] #initial condition
        sol = solve_ivp(fpar, [0,T], y0, method='Radau', rtol=1e-12, atol=1e-12)
        G, L = sol.y
        het, = ax.plot(G, L, ls='-', color='grey', zorder=10, label='Heteroclinic')
        if i == 0:
            arrowpos = [0.0003]
        else:
            arrowpos = [0.003,0.997]
        add_arrow_to_line2D(ax, het, arrow_locs=arrowpos,
                        arrowstyle='->')

        #plot saddle manifolds
        h = 3e-3 #length of manifold to plot
        fastidx = np.argmin(w) # negative eigenvalue is STABLE manifold
        slowidx = np.argmax(w) # positive eigenvalue is UNSTABLE manifold
        stablemanifold, = ax.plot([saddle.x[0],saddle.x[0]+(h*v[:,fastidx])[0]],[saddle.x[1],saddle.x[1]+(h*v[:,fastidx])[1]],'C1',lw=3,label='Stable Manifold')
        ax.plot([saddle.x[0],saddle.x[0]-(h*v[:,fastidx])[0]],[saddle.x[1],saddle.x[1]-(h*v[:,fastidx])[1]],'C1',lw=3)
        #plot slow eigenvector
        unstablemanifold, = ax.plot([saddle.x[0],saddle.x[0]+(h*v[:,slowidx])[0]],[saddle.x[1],saddle.x[1]+(h*v[:,slowidx])[1]],'C1',label='Unstable Manifold')
        ax.plot([saddle.x[0],saddle.x[0]-(h*v[:,slowidx])[0]],[saddle.x[1],saddle.x[1]-(h*v[:,slowidx])[1]],'C1')

        #get eigenvalues and eigenvectors of stable node
        w, v = np.linalg.eig(fgrad(node.x))

        #check if not complex eigevalues
        if (np.imag(w[0]) == 0):
            fastidx = np.argmin(w) #most negative eigenvalue is fast manifold
            slowidx = np.argmax(w) #least negative eigenvalue is slow manifold

            #plot fast eigenvector
            fastmanifold, = ax.plot([node.x[0],node.x[0]+(h*v[:,fastidx])[0]],[node.x[1],node.x[1]+(h*v[:,fastidx])[1]],'C2',lw=3,label='Fast Manifold')
            ax.plot([node.x[0],node.x[0]-(h*v[:,fastidx])[0]],[node.x[1],node.x[1]-(h*v[:,fastidx])[1]],'C2',lw=3)
            #plot slow eigenvector
            slowmanifold, = ax.plot([node.x[0],node.x[0]+(h*v[:,slowidx])[0]],[node.x[1],node.x[1]+(h*v[:,slowidx])[1]],'C2',label='Slow Manifold')
            ax.plot([node.x[0],node.x[0]-(h*v[:,slowidx])[0]],[node.x[1],node.x[1]-(h*v[:,slowidx])[1]],'C2')
        else:
            ax.text(x=0.5,y=0.1,s='complex evals',transform=ax.transAxes)

        #zoom in
        ax.set_xlim([0.87,0.89])
        ax.set_ylim([0.259,0.262])

    axes[1].set_xlabel('$G$')
    axes[0].set_ylabel('$L$')

    #legend
    axes[2].legend(handles=[het,saddleeq,stablemanifold,unstablemanifold,nodeeq,fastmanifold,slowmanifold],
                   labels=['Trajectory','Saddle','Stable Dir.','Unstable Dir.','Node','Fast Dir.','Slow Dir.'],
                   bbox_to_anchor=(1.02,0.5), loc=6)

    axes[0].text(-34/72, 1.0+2/72, 'A', transform=axes[0].transAxes, va='bottom')
    axes[1].text(-9/72, 1.0+2/72, 'B', transform=axes[1].transAxes, va='bottom')
    axes[2].text(-9/72, 1.0+2/72, 'C', transform=axes[2].transAxes, va='bottom')

    #savefig
    fig.savefig('pplane_beta={:0.12f}_b={:0.12f}.pdf'.format(p['beta'],p['b']),format='pdf')
