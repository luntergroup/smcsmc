from CGATReport.Tracker import TrackerSQL
import sys
import inspect, os

sys.path.extend( [ "../experiments"] )
from constpopsize_bylength_experiment import experiment_name as experiment_bylength
from constpopsize_byparticle_experiment import experiment_name as experiment_byparticle


class ConstpopsizeLengthdependence(TrackerSQL):

    # table to query (all -- but this is not used I think)
    #pattern = "(.*)$"

    # tracks and slices
    def getTracks(self):
        return ["missingData","fixedData","variableData"]

    def getSlices(self):
        # use the epoch start times as slices, so that the lengths are represented as columns
        statement = "SELECT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id WHERE name = '{}'".format(experiment_bylength)
        times = sorted(list(set( self.getValues(statement) )))
        return [ "T"+str(int(t)) for t in times ]
    
    def __call__(self, track, slice):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:])
        where = "experiment.name = '{}' AND result.start = {} AND type = 'Coal'".format(experiment_bylength, time)
        if track == "missingData": where += " AND experiment.missing_leaves = '[0, 1]'"
        elif track == "fixedData": where += " AND experiment.dataseed = 1"
        else:                      where += " AND experiment.missing_leaves = '[]' AND experiment.dataseed = infseed"

        # get last iteration
        statement = "SELECT result.iter " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE {}".format(where)
        lastiter = sorted( self.getValues(statement) )[-1]

        # get sequence length and Ne estimates at the last iteration
        statement = "SELECT experiment.sequence_length, result.ne " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE result.iter = {} AND {}".format(lastiter, where)
        values = self.get(statement)

        # extract the unique sequence lengths to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( v[0] for v in values )))
        results = {}
        for sl in keys:
            results[ int(sl) ] = [ne for (sl0, ne) in values if sl0 == sl]
        return results
        

class ConstpopsizeParticledependence(TrackerSQL):

    # tracks and slices
    def getTracks(self):
        return ["variableData"]

    def getSlices(self):
        # use the epoch start times as slices, so that the lengths are represented as columns
        statement = "SELECT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id WHERE name = '{}'".format(experiment_byparticle)
        times = sorted(list(set( self.getValues(statement) )))
        return [ "T"+str(int(t)) for t in times ]
    
    def __call__(self, track, slice):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:])
        where = "experiment.name = '{}' AND result.start = {} AND type = 'Coal'".format(experiment_byparticle, time)
        if track == "fixedData": where += " AND experiment.dataseed = 1"
        else:                    where += " AND experiment.dataseed = infseed"

        # get last iteration
        statement = "SELECT result.iter " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE {}".format(where)
        lastiter = sorted( self.getValues(statement) )[-1]
        
        # get sequence length and Ne estimates at the last iteration
        statement = "SELECT experiment.np, result.ne " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE result.iter = {} AND {}".format(lastiter, where)
        values = self.get(statement)

        # extract the unique particle counts to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( v[0] for v in values )))
        results = {}
        for np in keys:
            results[ int(np) ] = [ne for (np0, ne) in values if np0 == np]
        return results