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
               ['dataseed','sequence_seed',0],
               ['infseed','inference_seed',0],
               ['missing_leaves','missing_leaves',0] ]
               
    def getSlices(self):
        return [ field[1] for field in self.fields ]

    def __call__(self, track, slice):
        # generate the selector ('where') clause for this experiment, track and slice

        sqlfield, _, code = [ field for field in self.fields if field[1] == slice ][0]
        statement = "SELECT {} FROM experiment WHERE name = '{}'".format(sqlfield, track)
        values = self.getValues(statement)

        if code == -1:
            return { 'value(s)' : "(Typ:) " + str(values[0]) }

        numvalues = len(set(values))

        if numvalues == 1:
            return { 'value(s)' : values[0] }

        if numvalues < 4:
            return { 'value(s)' : ','.join( map( str, sorted(list(set(values))) )) }
        
        return { 'value(s)' : '{} - {}'.format( min(values), max(values) ) }
