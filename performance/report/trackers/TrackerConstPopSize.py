from CGATReport.Tracker import TrackerSQL
import sys
import inspect, os

sys.path.extend( [ "../experiments"] )
from constpopsize_bylength_experiment     import experiment_name as experiment_bylength
from constpopsize_byparticle_experiment   import experiment_name as experiment_byparticle
from constpopsize_bylag_experiment        import experiment_name as experiment_bylag
from constpopsize_byparticles2_experiment import experiment_name as experiment_byparticles2
from constpopsize_byparticles3_experiment import experiment_name as experiment_byparticles3

class ConstpopsizeLengthdependence(TrackerSQL):

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
        return ["variableData","fixedData"]

    def getSlices(self):
        # use the epoch start times as slices, so that the lengths are represented as columns
        statement = "SELECT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id WHERE name = '{}'".format(experiment_byparticle)
        times = sorted(list(set( self.getValues(statement) )))
        return [ "T"+str(int(t)) for t in times ]
    
    def __call__(self, track, slice):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:])
        where = "experiment.name = '{}' AND result.start = {} AND type = 'Coal'".format(experiment_byparticle, time)
        if track == "fixedData": where += " AND experiment.dataseed = 100"
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

        # failsafe: if just one result, return none
        for key in results:
            if len(results[key]) == 1:
                results[key] = []
                
        return results



class ConstpopsizeLagdependence(TrackerSQL):

    # tracks and slices
    def getTracks(self):
        return ["fixedData","variableData"]

    def getSlices(self):
        # use the epoch start times as slices, so that the lengths are represented as columns
        statement = "SELECT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id WHERE name = '{}'".format(experiment_bylag)
        times = sorted(list(set( self.getValues(statement) )))
        return [ "T"+str(int(t)) for t in times ]
    
    def __call__(self, track, slice):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:])
        where = "experiment.name = '{}' AND result.start = {} AND type = 'Coal'".format(experiment_bylag, time)
        if track == "fixedData": where += " AND experiment.dataseed = 100"
        else:                    where += " AND experiment.dataseed = infseed"

        # get last iteration
        statement = "SELECT result.iter " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE {}".format(where)
        lastiter = sorted( self.getValues(statement) )[-1]
        
        # get lag and Ne estimates at the last iteration
        statement = "SELECT experiment.lag, result.ne " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE result.iter = {} AND {}".format(lastiter, where)
        values = self.get(statement)

        # extract the lags to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( v[0] for v in values )))
        results = {}
        for lag in keys:
            results[ lag ] = [ne for (lag0, ne) in values if lag0 == lag]
        return results



class ConstpopsizeEMConvergence(TrackerSQL):

    options = "name = '{}'".format(experiment_bylag)
    trackfield = "lag"
    
    def getTracks(self):
        # use particle counts as tracks
        statement = "SELECT DISTINCT {} FROM experiment WHERE {}".format(self.trackfield,self.options)
        return sorted( self.getValues(statement) )

    def getSlices(self):
        # use the epoch start times as slices, so that the EM iterations are represented as columns
        statement = "SELECT DISTINCT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id " \
                    "WHERE {}".format(self.options)
        times = sorted( self.getValues(statement) )
        return [ "T{}".format(int(t)) for t in times ]
    
    def __call__(self, track, slice, **kwargs):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:])
        lag = float(track)
        where = "result.start = {} "\
                "AND {} = {} "\
                "AND type = 'Coal' AND {}".format(time, self.trackfield, lag, self.options)

        # get iteration and Ne estimates
        statement = "SELECT result.iter, result.ne " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE {}".format(where)
        values = self.get(statement)

        # extract the lags to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( v[0] for v in values )))
        return { key : [ne for (key0, ne) in values if key0 == key]
                 for key in keys }



class ConstpopsizeEMConvergence2(TrackerSQL):

    options = "name = '{}'".format(experiment_byparticles2)
    trackfield = "np"
    
    def getTracks(self):
        # use particle counts as tracks
        statement = "SELECT DISTINCT {} FROM experiment WHERE {}".format(self.trackfield,self.options)
        return sorted( self.getValues(statement) )

    def getSlices(self):
        # use the epoch start times as slices, so that the EM iterations are represented as columns
        statement = "SELECT DISTINCT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id " \
                    "WHERE {}".format(self.options)
        times = sorted( self.getValues(statement) )
        return [ "T{}".format(int(t)) for t in times ]
    
    def __call__(self, track, slice, **kwargs):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:])
        np = int(track)
        where = "result.start = {} "\
                "AND {} = {} "\
                "AND type = 'Coal' AND {}".format(time, self.trackfield, np, self.options)

        # get iteration and Ne estimates
        statement = "SELECT result.iter, result.ne " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE {}".format(where)
        values = self.get(statement)

        # extract the lags to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( v[0] for v in values )))
        return { key : [ne for (key0, ne) in values if key0 == key]
                 for key in keys }
    
    


class ConstpopsizeParticles2dependence(TrackerSQL):

    # tracks and slices
    def getTracks(self):
        return ["variableData"]

    def getSlices(self):
        # use the epoch start times as slices, leaving numbers of particles to be represented as columns
        statement = "SELECT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id " \
                    "WHERE name = '{}'".format(experiment_byparticles2)
        times = sorted(list(set( self.getValues(statement) )))
        return [ "T"+str(int(t)) for t in times ]
    
    def __call__(self, track, slice, options = None):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:])
        where = "experiment.name = '{}' AND result.start = {} AND type = 'Coal'".format(experiment_byparticles2, time)

        # get last iteration
        statement = "SELECT result.iter FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE {}".format(where)
        lastiter = sorted( self.getValues(statement) )[-1]
        
        # get lag and Ne estimates at the last iteration
        statement = "SELECT experiment.np, result.ne FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE result.iter = {} AND {}".format(lastiter, where)
        values = self.get(statement)

        # extract the particle counts to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( v[0] for v in values )))
        results = {}
        for np in keys:
            results[ np ] = [ne for (np0, ne) in values if np0 == np]
        return results


class ConstpopsizeParticles3dependence(TrackerSQL):

    # tracks and slices
    def getTracks(self):
        return ["variableData"]

    def getSlices(self):
        # use the epoch start times as slices, leaving numbers of particles to be represented as columns
        statement = "SELECT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id " \
                    "WHERE name = '{}'".format(experiment_byparticles3)
        times = sorted(list(set( self.getValues(statement) )))
        return [ "T"+str(int(t)) for t in times ]
    
    def __call__(self, track, slice, options = None):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:])
        where = "experiment.name = '{}' AND result.start = {} AND type = 'Coal'".format(experiment_byparticles3, time)

        # get last iteration
        statement = "SELECT result.iter FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE {}".format(where)
        lastiter = sorted( self.getValues(statement) )[-1]
        
        # get lag and Ne estimates at the last iteration
        statement = "SELECT experiment.np, result.ne FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE result.iter = {} AND {}".format(lastiter, where)
        values = self.get(statement)

        # extract the particle counts to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( v[0] for v in values )))
        results = {}
        for np in keys:
            results[ np ] = [ne for (np0, ne) in values if np0 == np]
        return results

    


