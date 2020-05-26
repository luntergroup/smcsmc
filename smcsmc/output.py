import pandas as pd # this is slow but I'm not sure what I'll need. 
import pdb

class Event:
   def __init__(self, line, header):
        self.line = line.strip().split() 
        self.records = dict(zip(header, self.line))
        

class Output:
    """Contains all the information from the output of SMC2, as well as convenience functions for interacting with it."""
    def __init__(self, path, N0 = 14312, g=29):
        self.df = pd.read_csv(path, sep='\s+')
        self.N0 = N0 # for downstream stuff
        self.g = g  # for downstream stuff
        #self.records = []
        #with open(path, 'r') as f:
        #    self.header = f.readline().strip().split()
        #    for line in f:
        #        self.records.append(Event(line, self.header))

    def subset_time(self, start, finish, bleed = False, final_iter = True):
        """
        Take a subset of the output file

        Parameters
        ----------
        start: float
            Begining of the interval you wish to subset
        end: float
            End of the interval you wish to subset
        bleed: bool
            Should additional records be included, encompasing the full range of the interval?
        final_iter: bool
            Include only the last iteration. 
        """

        f = self.df

        if final_iter:
            f = f.loc[f['Iter'] == max(f['Iter'])]

        f2 = f[f['Start']*self.g > start]
        f2 = f2[f2['End']*self.g < finish]
        
        if bleed:
            lower_end = f[f['Epoch'] == min(f2['Epoch'])-1]
            lower_end.loc[:, 'Start'] = start/self.g
            upper_end = f[f['Epoch'] == max(f2['Epoch'])+1]
            upper_end.loc[:, 'End'] = finish/self.g
            f = pd.concat([lower_end, f2, upper_end])

        self.filtered = f
        #self.filtered.reset_index()

    def integrate_migration(self, start, finish, bleed = True, From = 1):        

        prop = 0
        self.subset_time(start = start, finish = finish, bleed = bleed)
        current = self.filtered[self.filtered['Type'] == 'Migr']
        current = current.loc[self.filtered['From']==From]
        self.current = current

        self.current.sort_values(by = ['Epoch'], ascending = False, inplace=True)

        #pdb.set_trace()
        current = current.sort_values(by = ['Start'], ascending = False)
        rates = sum([[c[1]['Rate']] * round(c[1]['End']-c[1]['Start']) for c in current.iterrows()], [])

        for rate in rates: 
            prop += rate*(1-prop)

        return prop

 
        #for mig in self.current.iterrows():
        #    gen = mig[1]['End'] - mig[1]['Start']
        #    prop = prop + mig[1]['Rate'] * gen
            #pdb.set_trace()

        #return(prop)
