import math
import sys
import subprocess
import smcsmc.populationmodels
import pandas as pd
 
class Simulation:
    def __init__(self, L, midpoint, duration, proportion, flatarg = False):

        if isinstance(midpoint, set):
            midpoint = int(*midpoint)
            duration = int(*duration)
            proportion = float(*proportion)

        g0 = 133.0
        g1 = 133016.0
        g_split = 200000/29
        epochs = 27
        N0 = 14312
        self.samples = 4
        eps = 0.99  # make sure we set the new pop size / migration rate just before the change time
        mid = math.exp( math.log(g1/g0)/(2*(epochs-1)) ) / eps
        set_times = [ g0 * eps * math.exp( (math.log(g1/g0)*i) / (epochs-1) ) for i in range(epochs) ]
        eval_times = [ time * mid for time in set_times ]
        eval_times = [ set_times[0] / mid ] + eval_times
        set_times = [0] + set_times

        mu = 1.25e-8
        rho = 3e-9
        #L = 13223520
        g=29

        self.L = L
        self.midpoint = midpoint / g
        self.duration = duration / g
        self.proportion = proportion
       
        #if type == "scrm2seg":
         #   p = populationmodels.Population( num_samples = self.samples, sequence_length = L )
         #   infilename = "/dev/stdin"
         #   outfilename = "/dev/stdout"
         #   missing_leaves = []
         #   phased = True
         #   p.convert_scrm_to_seg(infilename, outfilename, missing_leaves, phased )
         #   sys.exit(0)

        self.model = ["scrm {} 1".format(self.samples),
                 "-l 100000 -p 10",  ## do not seed!
                 "-t {} -r {} {}".format(self.L * mu * 4 *N0, self.L * rho * 4 * N0, self.L),
                 "-I 2 {} {}".format(self.samples//2,self.samples//2)]

        split = False

        self.params = {}
        self.params['time'] = []
        self.params['ceu_ne'] = []
        self.params['yri_ne'] = []
        self.params['ceu_m'] = []
        self.params['yri_m'] = []

        for i,g_set in enumerate(set_times):
            g_eval = eval_times[i]

            # set model parameters for standard (P) model
            unscaled_time_set = g_set / (4*N0)
            ceu_popsize_unscaled = self.ceu(g_eval) / N0
            yri_popsize_unscaled = self.yri(g_eval) / N0

            ceu_migr_unscaled = self.migr_ceu(g_eval) * 4 * N0
            yri_migr_unscaled = self.new_migr_yri(g_eval,  midpoint = self.midpoint, years = self.duration, total = self.proportion, flat=flatarg) * 4 * N0

            self.model.append("-en {} 1 {}".format(unscaled_time_set, ceu_popsize_unscaled))
            self.model.append("-en {} 2 {}".format(unscaled_time_set, yri_popsize_unscaled))
            if g_set < g_split:
                self.model.append("-ema {} 0.0 {} {} 0.0".format(unscaled_time_set, ceu_migr_unscaled, yri_migr_unscaled))
            elif not split:
                self.model.append("-ej {} 2 1".format(unscaled_time_set))
                split = True

            self.params['time'].append(unscaled_time_set * 4 * N0)
            self.params['ceu_ne'].append(ceu_popsize_unscaled)
            self.params['yri_ne'].append(yri_popsize_unscaled)
            self.params['ceu_m'].append(ceu_migr_unscaled)
            self.params['yri_m'].append(yri_migr_unscaled)

        self.df = pd.DataFrame(self.params)

        self.model = " ".join(self.model)

    def write_df(self, output):
        self.df.to_csv(output, index = False)

    def F(self, x,a,b,c,d,e,f,g):
        x = math.log(x * 29)  # transform to log years
        return math.exp((a*x*x + b*x + c)*(d*x + e - math.cos( (x - f) / g ) ))

    # model for yri population size
    # x = time in generations
    def yri(self, x):
        return self.ceu( max(5555, x) )

    
    # model for ceu population size
    # x = time in generations
    def ceu(self, x):
        x = min(4100000/29, x)
        a = 0.02869
        b = -0.5555
        c = 3.016
        d = -3.651
        e = 63.3
        f = 10.8
        g = 0.51
        return self.F(x,a,b,c,d,e,f,g)

    # model for migration from yri into ceu, forward in time
    # x = time in generations
    # returns migration rate per generation
    # (total expected number of migrations: 0.435)
    # (new: total expected number of migrations: 0.498)
    def migr_ceu(self, x):
        ## version of 9/4/2018
        if False:
            if x < 2600: return 0         # 75 kya
            if x > 6900: return 0.0002  # 200 kya
            return 0.0002 * (1.0 - (math.log(x)-math.log(6900))/(math.log(2600)-math.log(6900)))
        ## version of 10/4/2018 (model 2): push back migration into ceu a little.  Looks good now.
        if x < 3448: return 0         # 100 kya
        if x > 6900: return 0.00025   # 200 kya
        return 0.00025 * (1.0 - (math.log(x)-math.log(6900))/(math.log(3448)-math.log(6900)))
        

    # model for migration from ceu into yri, forward in time
    def migr_yri(self, x, variant):
        start, mid, end = 1379, 2069, 3103   # 40k, 60k, 90k; "model 3"
        strength = 0.001
        if variant == -1:
            start, mid, end = 1025, 1700, 2900  # 30k, 50k, 85k; "model 1", 9/4/2018
            strength = 0.001
        if variant == 1:
            start, mid, end = 1724, 2586, 3448  # 50k, 75k, 100k; "model 2", 10/4/2018, moved peak to around 75kya
            strength = 0.001
        if x < start: return 0
        if x > end: return 0
        if x < mid:
            return strength * (1.0 - (math.log(x)-math.log(mid)) / (math.log(start)-math.log(mid)))
        return strength * (1.0 - (math.log(x)-math.log(mid)) / (math.log(end)-math.log(mid)))

    """
    This is a new version of the migration rate which is 
            a) square, rather than peaked
            b) consistent in terms of the proportion of migrating individuals
            c) relatively systematic
    Instead of distinct cases, we have a continum on three axis
            a) The proportion of the population migrating (this is the integral under the 
                    migration curve)
            b) The midpoint of the migration
            c) The length of the migration. This functions as as continuous transition between a "pulse"
                    and "continuous" migration.

    Note that input times here must be in terms of generations.
    """
    def new_migr_yri(self, x, midpoint, years, total, flat = False):
            start, end = (midpoint - (years / 2)), (midpoint + (years / 2) ) 
            if flat: 
                    proportion = 0.00025
            else:
                    proportion = total / years 

            # If its before or after the migration return 0
            if x < start: return 0
            if x > end: return 0
            return proportion

    def run_scrm(self, output):
        self.model += f" > {output}"
        subprocess.run(self.model, shell = True)

    def scrm_to_seg(self,input, output):
        p = smcsmc.populationmodels.Population( num_samples = self.samples, sequence_length = self.L )
        infilename = input
        outfilename = output
        missing_leaves = []
        phased = True
        p.convert_scrm_to_seg(infilename, outfilename, missing_leaves, phased )



