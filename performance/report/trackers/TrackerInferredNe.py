from CGATReport.Tracker import TrackerSQL

import matplotlib
import matplotlib.cm
import matplotlib.pyplot as plt
import math


def median( x ) :
        if len(x) == 0:
                return 0
        if len(x) % 2 == 1:
                return x[ len(x) // 2 ]
        return (x[ len(x) // 2 ] + x[ len(x) // 2 + 1 ]) / 2.0


class PlotNe(TrackerSQL):
        
        def __init__(self, **kwargs):
                self.name =      kwargs.get("name")         # experiment name
                self.minne =     kwargs.get("minne", 100)   # minimum Ne value to plot
                self.maxne =     kwargs.get("maxne", 1e5)   # maximum Ne value to plot
                self.mint  =     kwargs.get("mint", 30 )    # minimum t (generations)
                self.maxt  =     kwargs.get("maxt", 80000 )
                self.bias =      kwargs.get("bias", "1.0" )
                self.guide =     kwargs.get("guide", "0.0" )
                self.truth =     kwargs.get("truth",None)
                self.show_conv = kwargs.get("showconv", "False" )
                self.cache = False
                TrackerSQL.__init__(self)
                
        # tracks and slices
        def getTracks(self):
                # use the number of particles as tracks
                statement = "SELECT DISTINCT experiment.np FROM experiment WHERE name = '{}' AND str_parameter = 'guide{}_bias{}_mstepTrue'".format(self.name,self.guide,self.bias)
                times = sorted(list(set( self.getValues(statement) )))
                return [ "P"+str(int(t)) for t in times ]
        
        def getSlices(self):
                return ["slice1"]
    
        def __call__(self, track, slice):
                np = int( track[1:] )

                if self.truth != None:
                        vals = map(float,self.truth.split(';'))
                        if min(vals) > 0:
                                truth_ne = [vals[i]*vals[0] for i in range(1,len(vals),2)]
                                truth_t = [0] + [4*vals[0]*vals[i] for i in range(2,len(vals),2)]
                        else:
                                # a series parameters for ms/scrm -eG statements; exponential growth for zigzag model
                                ne = vals[0] * vals[1]   # value at start of exponential growth period
                                truth_ne = [ne]          # initial linear part
                                truth_t = [0]
                                for i in range(2,len(vals)-1,2):
                                        a = vals[i+1]
                                        t0 = vals[i]
                                        t1 = vals[i+2]
                                        for j in range(10):
                                                tj = t0 + (j/10.0)*(t1-t0)
                                                truth_t.append( tj * 4 * vals[0] )
                                                truth_ne.append( ne * math.exp( -a * (tj - t0) ) )
                                        ne = ne * math.exp( -a * (t1 - t0) )
                
                # get tuples (start_time, iteration, ne, experiment_id)
                statement = "SELECT result.start, result.iter, result.ne, result.exp_id FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                            "WHERE experiment.name = '{}' AND str_parameter = 'guide{}_bias{}_mstepTrue' AND type = 'Coal' AND experiment.np = {}".format(
                                    self.name,self.guide,self.bias,np)
                values = self.get(statement)
                
                # calculate median across experiments and add as exp_id = 0
                starts = sorted(list(set( v[0] for v in values )))
                iters = sorted(list(set( v[1] for v in values )))
                for s in starts:
                        for i in iters:
                                nes = [ v[2] for v in values if v[0] == s and v[1] == i ]
                                if len(nes) > 0:
                                        values.append( [s,i,median(nes),0] )

                # select data to plot: either data + median at last iteration, or median across all iterations
                main_curve = [ (v[0],v[2]) for v in values if v[1]==max(iters) and v[3] == 0 ]
                if self.show_conv == "False":
                        experiments = sorted(list(set( v[3] for v in values )))
                        bg_curves = [ [(v[0],v[2]) for v in values if v[1]==max(iters) and v[3] == experiment]
                                      for experiment in experiments ]
                else:
                        bg_curves = [ [(v[0],v[2]) for v in values if v[1]==iteration and v[3] == 0]
                                      for iteration in iters ]

                # make plot
                colormap = matplotlib.cm.inferno
                fig=plt.figure()
                ax=plt.subplot(111)
                if self.truth != None:
                        plt.plot( truth_t, truth_ne,
                                  drawstyle='steps-post', color='black', lw=1 )
                for idx, c in enumerate(bg_curves):
                        if self.show_conv == "False":
                                col = 0.1
                        else:
                                col = 0.75 * (idx / float(len(bg_curves)))
                        plt.plot( [p[0] for p in c], [p[1] for p in c],
                                  drawstyle='steps-post', color=colormap(col), lw=1, linestyle=":" )
                plt.plot( [p[0] for p in main_curve], [p[1] for p in main_curve],
                          drawstyle='steps-post', color=colormap(0.8), lw=3 )
                plt.ylim( self.minne, self.maxne )
                plt.xlim( self.mint, self.maxt )
                plt.title("Inferred pop sizes")
                ax.set_xscale("log", nonposx='clip')
                ax.set_yscale("log", nonposy='clip')
                return {"text": "#$mpl %i$#" % fig.number}
                
                
