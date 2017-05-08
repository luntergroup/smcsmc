from CGATReport.Tracker import TrackerSQL
import sys
import inspect, os


class Experiment(TrackerSQL):

    # tracks and slices
    def getTracks(self):
        experiments = sorted(list(set( self.getValues( "SELECT name FROM experiment" ))))
        return experiments

    fields = [ ['simulate_command','simulate',-1],     # -1 = show 'typical' example
               ['inference_command','inference',-1],
               ['np','particles',0],                   # 0 = (numerical or string); show min/max or fixed value
               ['sequence_length','length',0],
               ['lag','lag',0],
               ['num_samples','samples',0],
               ['ancestral_aware','ancestral aware',0],
               ['phased','phased',0],
               ['infer_recombination','infer recombination',0],
               ['bias_heights','bias heights',-1],
               ['bias_strengths','bias strengths',-1],
               ['guided_recomb','guided recombination',0],
               ['dataseed','sequence_seed',0],
               ['infseed','inference_seed',0],
               ['missing_leaves','missing_leaves',0],
               ['smcsmc_version','smcsmc_version',0],
               ['smcsmc_runtime','runtime',1] ]        # 1 = show min/max (seconds) and total (hours)

    def getSlices(self):
        return [ field[1] for field in self.fields ]

    def __call__(self, track, slice):
        # generate the selector ('where') clause for this experiment, track and slice

        sqlfield, _, code = [ field for field in self.fields if field[1] == slice ][0]
        statement = "SELECT {} FROM experiment WHERE name = '{}'".format(sqlfield, track)
        values = self.getValues(statement)

        if code == -1:
            return { 'value(s)' : "(Typ:) " + str(values[0]) }

        if code == 1:
            return { 'value(s)' : '{:1.2f} hr  ({:1.1f} - {:1.1f} sec)'.format( sum(values)/3600, min(values), max(values) ) }
        
        numvalues = len(set(values))

        if numvalues == 1:
            return { 'value(s)' : values[0] }

        if numvalues < 4:
            return { 'value(s)' : ','.join( map( str, sorted(list(set(values))) )) }

        return { 'value(s)' : '{} - {}'.format( min(values), max(values) ) }
