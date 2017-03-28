from CGATReport.Tracker import TrackerSQL
import sys
import inspect, os

sys.path.extend( [ "../experiments"] )
from twopops_unidirmigr_ancestralaware_experiment import experiment_name as experiment_ancestralaware


class TwopopsUnidirmigrAncestralawaredependence_bynp(TrackerSQL):

    # tracks and slices
    def getTracks(self):
        return ["Np100","Np500","Np1000","Np5000","Np100_8s","Np500_8s","Np1000_8s","Np5000_8s"]

    def getSlices(self):
        # use the epoch start times as slices, so that the ancestral aware boolean is represented as columns
        statement = "SELECT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id WHERE name = '{}'".format(experiment_ancestralaware)
        times = sorted(list(set( self.getValues(statement) )))
        return [ "T"+str(int(t)) for t in times ]

    def __call__(self, track, slice):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:]) # The 800 in e.g. T800
        where = "experiment.name = '{}' AND result.start = {} AND (type = 'Coal' or type = 'Migr' )".format(experiment_ancestralaware, time)
        if   track[:6] == "Np5000":  where += " AND experiment.np = 5000"
        elif track[:6] == "Np1000":  where += " AND experiment.np = 1000"
        elif track[:5] == "Np500":   where += " AND experiment.np = 500"
        elif track[:5] == "Np100":   where += " AND experiment.np = 100"
        if track[-3:] == "_8s":      where += " AND experiment.num_samples=8"
        else:                        where += " AND experiment.num_samples=4"

        # get ancestral_aware and Ne estimates at the first iteration
        statement = "SELECT experiment.ancestral_aware, result.frm, result.type, result.ne, result.rate " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE result.iter = 0 AND {}".format(where)
        values = self.get(statement)

        # extract the aa boolean to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( 'AA'+str(v[0])+'Type'+v[2]+'Pop'+str(v[1]) for v in values )))

        results = {}
        for key in keys:
            if key[7:11] == "Coal":
                results[ key ] = [ne for (aa, pop, type, ne, rate) in values if str(aa) == key[2] and str(pop) == key[-1] and type == key[7:11]]
            elif key[7:11] == "Migr":
                results[ key ] = [rate for (aa, pop, type, ne, rate) in values if str(aa) == key[2] and str(pop) == key[-1] and type == key[7:11]]

        return results

