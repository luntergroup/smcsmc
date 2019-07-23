import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import msprime

try:
    import stdpopsim
except:
    stdpopsim = None

matplotlib.use('Agg')

def plot_migration(input='result.out', output='result.png', g = 30, ymax=0.00025):
    """
    Plot the migration between *two* groups. 

    .. admonition::
        It would make a lot more sense if this could plot any number of migrations, but at present it can only do two.

    The function pulls from the input, which must be some sort of aggregated output from :code:`smcsmc`, and plots the migration rate over time for a given 
    generational time :code:`g`, which scales the number of generations to years. The Y axis (which can be modified with the :code:`ymax` parameter) represents 
    :math:`m_{ij}`, or the proportion of population j being placed by migration from population i *backwards in time*. This is important, as most intution about migration
    is understood forward in time. 

    This is a barebones function and there is no option to provide a "truth" bar. See the other plotting functions for various definitions of "truth".

    As an example, here is the resulting (symmetric) migration from a simulation:

    .. figure:: ../img/plot_migration.png

    .. todo::

        Typo above, it is not log years, rather years and the axis is logged.

    :param str input: This must be the file path to some sort of aggregated output from :code:`smcsmc`. 
                        This means that it can be either a :code:`chunkfinal.out` or :code:`result.out` but *not* an
                        individual chunk. We need all the data here.
    :param str output: Filepath to save the plot.
    :param int g: The length of one generation in years. This is used to scale the x axis.
    :param float ymax: The maximum y value to plot. This is used to scale the plots up or down. 
    """
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
    This function is very similar to :meth:`smcsmc.plot.plot_migration`  except that it includes the ability to add a bar of "truth". In this case, the function uses a specific form of "truth"
    generated from recording all epochs output by :code:`SCRM`. Additionally, we provide both the effective population size and the migration rates.

    The structure of the truth guide is like so:
    
    +------------+----------+----------+---------+---------+
    | Start Time | Pop_1 Ne | Pop_2 Ne | Pop_1 M | Pop_2 M |
    +------------+----------+----------+---------+---------+
    | 0          | 3        | 3        | 6       | 2       |
    +------------+----------+----------+---------+---------+
    | 100        | 0.65     | 0.3      | 14      | 0       |
    +------------+----------+----------+---------+---------+
    | ...        | .        | .        | .       | .       |
    +------------+----------+----------+---------+---------+
    | 10000      | 0.4      | 0.4      | 8       | 6       |
    +------------+----------+----------+---------+---------+
   
    Saved as a CSV file for the :code:`guide` argument.

    .. todo::
        This function is part of a WIP tutorial on simulating with :code:`SCRM`. More details and convenience functions to come.

    :param str input: The full file path to an aggregated output file from :code:`smcsmc`.
    :param str guide: The full file path to a CSV formated as above with the truth of a simulation. 
    :param str output: Filepath to save the plot.
    :param int g: The length of one generation in years. This is used to scale the x axis.
    :param float ymax: The maximum y value to plot. This is used to scale the plots up or down.  
    """
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
    Creates a plot of all iterations by colour. This plot is useful for assessing convergence.

    :param str input: The full file path to an aggregated output file from :code:`smcsmc`.
    :param str output: Filepath to save the plot.
    :param int g: The length of one generation in years. This is used to scale the x axis.
    :param stdpopsim.Model model: Model for plotting.
    :param float ymax: The maximum y value to plot. This is used to scale the plots up or down.  
    :param int steps: Don't worry about this.
    :param int pop_id: If your model includes multiple populations, which one do you want to plot?

    .. figure:: ../img/rainbow.png
        :align: center
 """

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
    
