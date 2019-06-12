import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import msprime

try:
    import stdpopsim
except:
    stdpopsim = None

#matplotlib.use('Agg')

def plot_migration(input='result.out', output='result.png', g = 30, ymax=0.00025):
    """
    Plot migration between two groups.

    Give the path to the result file and the path to the output file. The ymax can
    be adjusted for larger or smaller ranges."""
    df = pd.read_csv(input, sep = "\s+")
    df = df[df['Iter'] == max(df['Iter'])][df['Clump'] == -1]
    zto = df[df['Type'] == 'Migr'][df['From'] == 0]
    otz = df[df['Type'] == 'Migr'][df['From'] == 1]

    plt.plot(zto['Start']*30, zto['Rate'], drawstyle = 'steps-pre')
    plt.plot(otz['Start']*30, otz['Rate'], drawstyle = 'steps-pre')

    plt.xscale('log')
    plt.ylim((0, ymax))
    plt.xlim((0, 5e5))
    plt.legend(['Population 0 to 1', 'Population 1 to 0'])
    plt.title('Migration')
    plt.xlabel(f"Log Years (g = {g})") 
    plt.ylabel(f"Migration Rate (4N0m)")
    plt.tight_layout()

    plt.show()

    plt.savefig(output) 
    
def plot_with_guide(input, guide, output, g = 30, ymax = 0.00025, N0=14312):
    """
    Plots migration along with a guide from a simulations."""
    df = pd.read_csv(input, sep = "\s+")
    df = df[df['Iter'] == max(df['Iter'])][df['Clump'] == -1]
    zto_m = df[df['Type'] == 'Migr'][df['From'] == 0]
    otz_m = df[df['Type'] == 'Migr'][df['From'] == 1]
    zto_ne = df[df['Type'] == 'Coal'][df['From'] == 0]
    otz_ne = df[df['Type'] == 'Coal'][df['From'] == 1]

    truth = pd.read_csv(guide)

    fig, (ax1, ax2) = plt.subplots(2,1,sharex=True,figsize=(10, 8))

    # Plot Ne
    #ax1.subplot(2,1,1)
    ax1.set_yscale('log')
    ax1.set_ylim((1e3, 1e6))
    ax1.plot(zto_ne['Start']*g, zto_ne['Ne'], drawstyle = 'steps-pre')
    ax1.plot(otz_ne['Start']*g, otz_ne['Ne'], drawstyle = 'steps-pre')
    ax1.plot(truth['time']*g, truth['ceu_ne']*N0, drawstyle = 'default', c = "black")
    ax1.plot(truth['time']*g, truth['yri_ne']*N0, drawstyle = 'default', c = "black")

    ax1.set_xscale('log')
    ax1.set_xlim((0, 1e6))
    #ax1.legend(['Population 0 to 1', 'Population 1 to 0'])
    ax1.set_title('Population Size')
    ax1.set_xlabel(f"Log Years (g = {g})") 
    ax1.set_ylabel(f"Migration Rate (4N0m)")
   # ax1.tight_layout()

    #ax2.subplot(2,1,2)
    ax2.plot(zto_m['Start']*30, zto_m['Rate'], drawstyle = 'steps-pre')
    ax2.plot(otz_m['Start']*30, otz_m['Rate'], drawstyle = 'steps-pre')
    ax2.plot(truth['time']*g, truth['ceu_m']/(4*N0), drawstyle = 'default', c = "black")
    ax2.plot(truth['time']*g, truth['yri_m']/(4*N0), drawstyle = 'default', c = "black")

    ax2.set_xscale('log')
    ax2.set_ylim((0, ymax))
    ax2.set_xlim((0, 1e6))
  #  ax2.set_legend(['Population 0 to 1', 'Population 1 to 0'])
    ax2.set_title('Migration')
    ax2.set_xlabel(f"Log Years (g = {g})") 
    ax2.set_ylabel(f"Migration Rate (4N0m)")
  #  ax2.tight_layout()

    fig.savefig(output)





def plot_rainbow(input, output, g = 30, model = None, steps = None, pop_id = 1):
    """ 
    Creates a plot of all iterations by colour.

    Give the path to the full result file and each 
    iteration will be displayed in a different colour.
    Useful for assessing convergence."""

    f, ax = plt.subplots(figsize=(7, 7))
    ax.set(xscale="log", yscale="log")
    ax.set_ylim([10e2, 10e4])
    ax.set_xlim([1e3, 1e6])

    if model == "ooa":
        if stdpopsim is not None:
            model = getattr(stdpopsim.homo_sapiens,"GutenkunstThreePopOutOfAfrica")()
        else:
            print ("Module stdpopsim not available - cannot go on.")
            sys.exit(1)

    if model is not None:
        ddb = msprime.DemographyDebugger(**model.asdict())
        if steps is None:
            end_time = ddb.epochs[-2].end_time + 10000
            steps = np.exp(np.linspace(1,np.log(end_time),31))
        num_samples = [0 for _ in range(ddb.num_populations)]
        num_samples[pop_id] = 20
        coal_rate, P = ddb.coalescence_rate_trajectory(steps=steps,
            num_samples=num_samples, double_step_validation=False)
        steps = steps * g
        ax.plot(steps, 1/(2*coal_rate), c="black", drawstyle = 'steps-pre')

    nt = pd.read_csv(input, sep = '\s+')
    nt['Start'] *= g
    for k, g in nt.groupby(['Iter']):
        g.plot(x='Start',y='Ne', ax = ax, drawstyle = 'steps-pre', label = k)
    ax.legend(title='EM Iterations', ncol = 5)
    f.savefig(output)
    
