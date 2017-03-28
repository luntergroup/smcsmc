from CGATReport.Tracker import TrackerSQL
import sys
import inspect, os

sys.path.extend( [ "../experiments"] )
from constpopsize_bylag_truestart_experiment import experiment_name as experiment_bylag
from constpopsize_bylag_experiment           import experiment_name as experiment_bylag_falsestart
from constpopsize_ancestralaware_experiment  import experiment_name as experiment_ancestralaware
from constpopsize_phasedependence_experiment import experiment_name as experiment_phase
from constpopsize_bylength_experiment        import experiment_name as experiment_bylength
from constpopsize_byparticle_experiment      import experiment_name as experiment_byparticle
#from constpopsize_bylag_experiment           import experiment_name as experiment_bylag
from constpopsize_byparticles2_experiment    import experiment_name as experiment_byparticles2
from constpopsize_byparticles3_experiment    import experiment_name as experiment_byparticles3



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


class ConstpopsizeLagdependence(TrackerSQL):

    # tracks and slices
    def getTracks(self):
        return ["2s","4s","8s"]

    def getSlices(self):
        # use the epoch start times as slices, so that the lengths are represented as columns
        statement = "SELECT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id WHERE name = '{}'".format(experiment_bylag)
        times = sorted(list(set( self.getValues(statement) )))
        return [ "T"+str(int(t)) for t in times ]

    def __call__(self, track, slice):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:])
        where = "experiment.name = '{}' AND result.start = {} AND type = 'Coal'".format(experiment_bylag, time)
        if track == "2s":    where += " AND experiment.num_samples = 2"
        elif track == "4s":  where += " AND experiment.num_samples = 4"
        elif track == "8s":  where += " AND experiment.num_samples = 8"
        else: raise ValueError("The track is undefined; should be '2s' '4s' or '8s'; but the track is set to "+track)

        # get sequence length and Ne estimates at the first iteration
        statement = "SELECT experiment.lag, result.ne " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE result.iter = 0 AND {}".format( where)
        values = self.get(statement)

        # extract the unique sequence lengths to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( v[0] for v in values )))

        results = {}
        for lag in keys:
            results[ lag ] = [ne for (lag0, ne) in values if lag0 == lag]
        return results

class ConstpopsizeLagdependence_lastiter(TrackerSQL):

    # tracks and slices
    def getTracks(self):
        return ["2s","4s","8s"]

    def getSlices(self):
        # use the epoch start times as slices, so that the lengths are represented as columns
        statement = "SELECT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id WHERE name = '{}'".format(experiment_bylag)
        times = sorted(list(set( self.getValues(statement) )))
        return [ "T"+str(int(t)) for t in times ]

    def __call__(self, track, slice):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:])
        where = "experiment.name = '{}' AND result.start = {} AND type = 'Coal'".format(experiment_bylag, time)
        if track == "2s":    where += " AND experiment.num_samples = 2"
        elif track == "4s":  where += " AND experiment.num_samples = 4"
        elif track == "8s":  where += " AND experiment.num_samples = 8"
        else: raise ValueError("The track is undefined; should be '2s' '4s' or '8s'; but the track is set to "+track)

        # get last iteration
        statement = "SELECT result.iter " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE {}".format(where)
        lastiter = sorted( self.getValues(statement) )[-1]

        # get sequence length and Ne estimates at the first iteration
        statement = "SELECT experiment.lag, result.ne " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE result.iter = {} AND {}".format(lastiter, where)
        values = self.get(statement)

        # extract the unique sequence lengths to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( v[0] for v in values )))

        results = {}
        for lag in keys:
            results[ lag ] = [ne for (lag0, ne) in values if lag0 == lag]
        return results

class ConstpopsizeLagdependenceFalsestart(TrackerSQL):

    # tracks and slices
    def getTracks(self):
        return ["2s","4s","8s"]

    def getSlices(self):
        # use the epoch start times as slices, so that the lengths are represented as columns
        statement = "SELECT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id WHERE name = '{}'".format(experiment_bylag_falsestart)
        times = sorted(list(set( self.getValues(statement) )))
        return [ "T"+str(int(t)) for t in times ]

    def __call__(self, track, slice):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:])
        where = "experiment.name = '{}' AND result.start = {} AND type = 'Coal'".format(experiment_bylag_falsestart, time)
        if track == "2s":    where += " AND experiment.num_samples = 2 AND np=500"
        elif track == "4s":  where += " AND experiment.num_samples = 4 AND np=500"
        elif track == "8s":  where += " AND experiment.num_samples = 8 AND np=500"
        else: raise ValueError("The track is undefined; should be '2s' '4s' or '8s'; but the track is set to "+track)

        # get sequence length and Ne estimates at the first iteration
        statement = "SELECT experiment.lag, result.ne " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE result.iter = 0 AND {}".format( where)
        values = self.get(statement)

        # extract the unique sequence lengths to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( v[0] for v in values )))

        results = {}
        for lag in keys:
            results[ lag ] = [ne for (lag0, ne) in values if lag0 == lag]
        return results

class ConstpopsizeLagdependenceFalsestart_lastiter(TrackerSQL):

    # tracks and slices
    def getTracks(self):
        return ["2s","4s","8s"]

    def getSlices(self):
        # use the epoch start times as slices, so that the lengths are represented as columns
        statement = "SELECT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id WHERE name = '{}'".format(experiment_bylag_falsestart)
        times = sorted(list(set( self.getValues(statement) )))
        return [ "T"+str(int(t)) for t in times ]

    def __call__(self, track, slice):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:])
        where = "experiment.name = '{}' AND result.start = {} AND type = 'Coal'".format(experiment_bylag_falsestart, time)
        if track == "2s":    where += " AND experiment.num_samples = 2 AND np=500"
        elif track == "4s":  where += " AND experiment.num_samples = 4 AND np=500"
        elif track == "8s":  where += " AND experiment.num_samples = 8 AND np=500"
        else: raise ValueError("The track is undefined; should be '2s' '4s' or '8s'; but the track is set to "+track)

        # get last iteration
        statement = "SELECT result.iter " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE {}".format(where)
        lastiter = sorted( self.getValues(statement) )[-1]

        # get sequence length and Ne estimates at the first iteration
        statement = "SELECT experiment.lag, result.ne " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE result.iter = {} AND {}".format(lastiter, where)
        values = self.get(statement)

        # extract the unique sequence lengths to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( v[0] for v in values )))

        results = {}
        for lag in keys:
            results[ lag ] = [ne for (lag0, ne) in values if lag0 == lag]
        return results


class ConstpopsizeAncestralawaredependence_bynp_multiters(TrackerSQL):

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
        where = "experiment.name = '{}' AND result.start = {} AND type = 'Coal'".format(experiment_ancestralaware, time)
	if   track[:6] == "Np5000":  where += " AND experiment.np=5000"
	elif track[:6] == "Np1000":  where += " AND experiment.np=1000"
	elif track[:5] == "Np500":   where += " AND experiment.np=500"
	elif track[:5] == "Np100":   where += " AND experiment.np=100"
	if track[-3:] == "_8s":      where += " AND experiment.num_samples=8"
	else:                        where += " AND experiment.num_samples=4"

	# I would like to plot the result from the first iteration and the result from the last iteration...
	#   values returns aa, ne; could return aa, iter, ne and then create a key that defines both aa and iter
	#   self.getValues(statement) returns a list of the first value of each row
	#   self.get(statement) returns a list of tuples
	#   e.g. self.get("SELECT column1, column2 FROM table")
        #         returns [(1,2),(2,4),(3,2)]

        # get ancestral_aware and Ne estimates at the last iteration
        statement = "SELECT experiment.ancestral_aware, result.iter, result.ne " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE (result.iter = 0 OR result.iter = 4) AND {}".format(where)
        values = self.get(statement)

        # extract the aa boolean to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( "Iter"+str(v[1])+"_AA"+str(v[0]) for v in values )))

        results = {}
        for key in keys:
        	results[ key ] = [ ne for (aa, iter, ne) in values if (str(aa) == key[-1] and str(iter) == key[4]) ]
        return results

class ConstpopsizeAncestralawaredependence_bynp(TrackerSQL):

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
        where = "experiment.name = '{}' AND result.start = {} AND type = 'Coal'".format(experiment_ancestralaware, time)
	if   track[:6] == "Np5000":  where += " AND experiment.np = 5000"
        elif track[:6] == "Np1000":  where += " AND experiment.np = 1000"
        elif track[:5] == "Np500":   where += " AND experiment.np = 500"
        elif track[:5] == "Np100":   where += " AND experiment.np = 100"
        if track[-3:] == "_8s":      where += " AND experiment.num_samples=8"
        else:                        where += " AND experiment.num_samples=4"

        # get ancestral_aware and Ne estimates at the first iteration
        statement = "SELECT experiment.ancestral_aware, result.ne " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE result.iter = 0 AND {}".format(where)
        values = self.get(statement)

	# extract the aa boolean to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( v[0] for v in values )))

        results = {}
        for sl in keys:
            results[ sl ] = [ne for (sl0, ne) in values if sl0 == sl]
        return results


class ConstpopsizePhasedependence_bynp(TrackerSQL):

    # tracks and slices
    def getTracks(self):
        return ["Np100","Np500","Np1000","Np5000","Np100_8s","Np500_8s","Np1000_8s","Np5000_8s"]

    def getSlices(self):
        # use the epoch start times as slices, so that the ancestral aware boolean is represented as columns
        statement = "SELECT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id WHERE name = '{}'".format(experiment_phase)
        times = sorted(list(set( self.getValues(statement) )))
        return [ "T"+str(int(t)) for t in times ]

    def __call__(self, track, slice):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:]) # The 800 in e.g. T800
        where = "experiment.name = '{}' AND result.start = {} AND type = 'Coal'".format(experiment_phase, time)
        if   track[:6] == "Np5000":  where += " AND experiment.np = 5000"
        elif track[:6] == "Np1000":  where += " AND experiment.np = 1000"
        elif track[:5] == "Np500":   where += " AND experiment.np = 500"
        elif track[:5] == "Np100":   where += " AND experiment.np = 100"
        if track[-3:] == "_8s":      where += " AND experiment.num_samples=8"
        else:                        where += " AND experiment.num_samples=4"

        # get ancestral_aware and Ne estimates at the first iteration
        statement = "SELECT experiment.phased, result.ne " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE result.iter = 0 AND {}".format(where)
        values = self.get(statement)

        # extract the aa boolean to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( v[0] for v in values )))

        results = {}
        for sl in keys:
            results[ sl ] = [ne for (sl0, ne) in values if sl0 == sl]
        return results

class ConstpopsizePhasedependence_bynp_multiters(TrackerSQL):

    # tracks and slices
    def getTracks(self):
        return ["Np100","Np500","Np1000","Np5000","Np100_8s","Np500_8s","Np1000_8s","Np5000_8s"]

    def getSlices(self):
        # use the epoch start times as slices, so that the ancestral aware boolean is represented as columns
        statement = "SELECT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id WHERE name = '{}'".format(experiment_phase)
        times = sorted(list(set( self.getValues(statement) )))
        return [ "T"+str(int(t)) for t in times ]

    def __call__(self, track, slice):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:]) # The 800 in e.g. T800
        where = "experiment.name = '{}' AND result.start = {} AND type = 'Coal'".format(experiment_phase, time)
        if   track[:6] == "Np5000":  where += " AND experiment.np = 5000"
        elif track[:6] == "Np1000":  where += " AND experiment.np = 1000"
        elif track[:5] == "Np500":   where += " AND experiment.np = 500"
        elif track[:5] == "Np100":   where += " AND experiment.np = 100"
        if track[-3:] == "_8s":      where += " AND experiment.num_samples=8"
        else:                        where += " AND experiment.num_samples=4"

        # get ancestral_aware and Ne estimates at the last iteration
        statement = "SELECT experiment.phased, result.iter, result.ne " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE (result.iter = 0 OR result.iter = 4) AND {}".format(where)
        values = self.get(statement)

        # extract the aa boolean to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( "Iter"+str(v[1])+"_Ph"+str(v[0]) for v in values )))

        results = {}
        for key in keys:
                results[ key ] = [ ne for (ph, iter, ne) in values if (str(ph) == key[-1] and str(iter) == key[4]) ]
        return results

class ConstpopsizeEMConvergence(TrackerSQL):

    def __init__(self, **kwargs):
        self.name =       kwargs.get("name",experiment_bylag)
        self.trackfield = kwargs.get("track","lag")
        self.cache =      False                    # this doesn't seem to work...
        TrackerSQL.__init__(self)
    
    def getTracks(self):
        # use lags or particle counts as tracks
        statement = "SELECT DISTINCT {} FROM experiment WHERE name = '{}'".format(self.trackfield,self.name)
        return sorted( self.getValues(statement) )

    def getSlices(self):
        # use the epoch start times as slices, so that the EM iterations are represented as columns
        statement = "SELECT DISTINCT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id " \
                    "WHERE name = '{}'".format(self.name)
        times = sorted( self.getValues(statement) )
        return [ "T{}".format(int(t)) for t in times ]
    
    def __call__(self, track, slice, **kwargs):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:])
        trackvalue = float(track)
        where = "result.start = {} "\
                "AND {} = {} "\
                "AND type = 'Coal' AND name = '{}'".format(time, self.trackfield, trackvalue, self.name)

        # get iteration and Ne estimates
        statement = "SELECT result.iter, result.ne " \
                    "FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE {}".format(where)
        values = self.get(statement)

        # extract the lags to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( v[0] for v in values )))
        return { key : [ne for (key0, ne) in values if key0 == key]
                 for key in keys }




class ConstpopsizeNe(TrackerSQL):

    def __init__(self, **kwargs):
        self.name =       kwargs.get("name",experiment_byparticle)
        self.field =      kwargs.get("column","np")
        self.cache =      False                    # this doesn't seem to work...
        TrackerSQL.__init__(self)

    # tracks and slices
    def getTracks(self):
        return ["variableData","fixedData"]

    def getSlices(self):
        # use the epoch start times as slices, leaving numbers of particles to be represented as columns
        statement = "SELECT DISTINCT result.start FROM result INNER JOIN experiment ON experiment.id = result.exp_id " \
                    "WHERE name = '{}'".format(self.name)
        times = sorted( self.getValues(statement) )
        return [ "T"+str(int(t)) for t in times ]
    
    def __call__(self, track, slice):
        # generate the selector ('where') clause for this experiment, track and slice
        time = float(slice[1:])
        where = "experiment.name = '{}' AND result.start = {} AND type = 'Coal'".format(self.name, time)
        if track == "fixedData": where += " AND experiment.dataseed = 100"
        else:                    where += " AND experiment.dataseed = infseed"

        # get last iteration
        statement = "SELECT result.iter FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE {}".format(where)
        lastiter = sorted( self.getValues(statement) )[-1]
        
        # get lag and Ne estimates at the last iteration
        statement = "SELECT experiment.{}, result.ne FROM experiment INNER JOIN result ON experiment.id = result.exp_id " \
                    "WHERE result.iter = {} AND {}".format(self.field, lastiter, where)
        values = self.get(statement)

        # extract the particle counts to serve as keys in the results, and build dictionary { key : [ne values] }
        keys = sorted(list(set( v[0] for v in values )))
        return { key : [ne for (key0, ne) in values if key0 == key]
                 for key in keys }



