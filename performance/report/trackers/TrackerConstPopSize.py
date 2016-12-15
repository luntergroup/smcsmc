from CGATReport.Tracker import TrackerSQL
import sys
import inspect, os

sys.path.extend( [ "../experiments"] )
#import constpopsize_bylength_cfg
import constpopsize_bylength_experiment


class ConstpopsizeLengthdependence(TrackerSQL):

    # table to query (all -- but this is not used I think)
    pattern = "(.*)$"

    # tracks and slices
    def getTracks(self):
        return ["missingData","fixedData","variableData"]

    def getSlices(self):
        # use the sequence lengths as slices; get them directly from module that generate them
        return map(int, constpopsize_bylength_experiment.seqlens)
    
    def __call__(self, track, slice):
        # find the experiments (smcsmc runs) for this track - take name directly from module that generated them
        where = "name = '{}' AND sequence_length = {}".format(constpopsize_bylength_experiment.experiment_name,
                                                              slice)
        if track == "missingData": where += " AND missing_leaves = '[0, 1]'"
        elif track == "fixedData": where += " AND dataseed = 1"
        else:                      where += " AND missing_leaves = '[]' AND dataseed = infseed"
        statement = "SELECT id FROM experiment WHERE {}".format(where)
        expids = self.getValues(statement)

        # get last iteration
        statement = "SELECT iter FROM result WHERE exp_id = {} AND type = 'Coal'".format(expids[0])
        lastiter = sorted( self.getValues(statement) )[-1]

        # get Ne estimates, at the last iteration, for all epochs (and ALL experiments)
        statement = "SELECT start, ne, exp_id FROM result WHERE type = 'Coal' AND iter = {}".format(lastiter)
        values = self.get(statement)

        # extract the unique start times, to serve as keys in the results
        keys = sorted(list(set( v[0] for v in values if v[2] in expids )))

        # build a dictionary { key: [ne values] }
        results = {}
        for t in keys:
            results[ "T" + str(int(t)) ] = [ne
                                            for (t0, ne, expid) in values
                                            if t0 == t and expid in expids]

        return results
        
    
