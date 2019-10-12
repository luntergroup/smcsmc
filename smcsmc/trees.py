#import smcsmc
import gzip
import pandas as pd
from collections import OrderedDict
from numpy import unique

class Event:
    def __init__(self, tokens):
        self.right = tokens[1]
        self.time = tokens[2]
        self.lineages = self.parse_lineages(tokens[5])
        self.migration_events = []
        self.nmig = 0

        self.nhaps = 4 
        self.pops = [0,0,1,1]
        self.ids = [0,0,1,1]

    def parse_lineages(self, lineages):
        out = [c for c in range(len(lineages)) if lineages[c] == '1']
        if len(out) == 0:
            out = 'R'

        return(out)

    def append(self, tokens):
        if tokens[0] == "C":
            self.coal = tokens[2]
            self.population = tokens[3]
            self.recipient = self.parse_lineages(tokens[5])
        
        elif tokens[0] == "M":
            self.migration_events.append(Migration(tokens))
            self.nmig += 1

    def print(self):
        print("Event at " + self.right + " in " + ''.join([str(int(s)) for s in self.lineages])+"\n" +
               "-\tTime: " + self.time + '\n' +
               '-\tCoal Time: ' + self.coal + '\n' +
               '-\tN Migration: ' + str(self.nmig))
        for e in self.migration_events:
            e.print()

    def initialiseTables(self):
        self.NodeTable = NodeTable()
        self.IndividualTable = IndividualTable()
        self.EdgeTable = EdgeTable()
        self.PopulationTable = PopulationTable()
        
        # Base entries
        for i in range(self.nhaps):
            self.NodeTable.table.loc[i] = ['1', '0', self.pops[i], self.ids[i], None] 
            self.IndividualTable.table.loc[self.ids[i]] = ['1', None, None]

        for i in range(len(unique(self.pops))):
            self.PopulationTable.add(None)

        # Have some entries for initialising 
        # internal nodes
        time1 = 50000
        time2 = 50000
        self.NodeTable.add(time1, 0, 0)
        self.NodeTable.add(time2, 1, 1)
        self.NodeTable.add(time1+time2, 0, -1)

        self.EdgeTable.add(self.right, self.left, self.nhaps, 0)
        self.EdgeTable.add(self.right, self.left, self.nhaps, 1) 
        self.EdgeTable.add(self.right, self.left, self.nhaps+1, 2)
        self.EdgeTable.add(self.right, self.left, self.nhaps+1, 3)
        self.EdgeTable.add(self.right, self.left, self.nhaps+2, self.nhaps)
        self.EdgeTable.add(self.right, self.left, self.nhaps+2, self.nhaps+1)

        self.tables = TableCollection({'Node': self.NodeTable,
                 'Individual': self.IndividualTable,
                 'Edge': self.EdgeTable,
                 'Population': self.PopulationTable})

    def printTables(self):   
        for v in self.tables.dict.values():
            v.print()

class Migration(Event):
    def __init__(self, tokens):
        assert(tokens[0] == "M")
        super().__init__(tokens) 
        self.To = tokens[3]
        self.From = tokens[4]
        self.migration_events = None

    def print(self):
        print("*\t[" + ''.join([str(s) for s in self.lineages]) + "] \tMigration " + self.From + " -> " + self.To + " @ " + str(self.time) + '\n')

        
class Record:
    def __init__(self, tree):
        self.events = OrderedDict()
        # These are some hardcoded parameters that should
        # be inferred from the record when I have a moment.
        self.nhaps = 4 
        self.pops = [0,0,1,1]
        self.ids = [0,0,1,1]

        with gzip.open(tree, 'rb') as f:
            current_event = Event([s.decode("utf-8") for s in f.readline().strip().split()])
            for line in f:
                tokens = [s.decode("utf-8") for s in line.strip().split(b'\t')]
                if tokens[0] == "R": 
                    self.addEvent(current_event)
                    current_event = Event(tokens)
                    # This update is based on the fact
                    # that the tree files decrease in position
                    # Should the opposite occur,
                    # this would be backwards.
                    self.events[list(self.events.keys())[-1]].left = current_event.right
                    #assert(self.events[list(self.events.keys())[-1]].right < current_event.left) 
                else:
                    current_event.append(tokens)

    def addEvent(self, event):
        self.events[event.right] =event

 
class TableCollection:
    def __init__(self, dict):
        self.dict = dict

    def update(self, entry):
        # 1. Add a coalescent node
        entry.tables['Node'].add()
        # 2. Attach the proper nodes to the new one
        # 3. Find the entries that were cut off (if any)
        # 4. Destroy their parents 
        # 5. Check that it still forms a connected tree

    def return_parents(self,node):
        return(self.parents(node, table = self.dict['Edge'].table))

    def parents(self,node, table): 
        if int(self.dict['Node'].table.loc[node]['individual']) == -1: 
            return ["Root"]
        else: 
            #p2.extend([int(table.loc[table['child'] == node]['parent'])])
            new_node = int(table.loc[table['child'] == node]['parent'])
            return [new_node] + self.parents(new_node,  table = table)

class Table: 
    def __init__(self, columns):
        self.table = pd.DataFrame(columns = columns) 
        self.name = "Generic"

    def print(self):
        print("\n*** " + self.name )
        print(self.table) 


class NodeTable(Table):
    def __init__(self): 
        super().__init__(columns = ['flags', 'time', 'population', 'individual', 'metadata'])
        self.name = "Node Table"

    def add(self, time, population, individual):
        row = {'flags': 0, 'time': time, 'population': population, 'individual': individual, 'metadata': None}
        self.table = self.table.append(pd.DataFrame(row, index = [len(self.table.index)]))
        self.table.reset_index()
        #return(len(self.table.index))


class IndividualTable(Table):
    def __init__(self):
        super().__init__(columns = ['flags', 'location', 'metadata'])
        self.name = "Individual Table"

class EdgeTable(Table): 
    def __init__(self):
        self.columns = ['left', 'right', 'parent', 'child']
        super().__init__(columns = self.columns)
        self.name = "Edge Table"

    def add(self, left, right, parent, child):
        row = {'left': left, 'right': right, 'parent': parent, 'child': child}
        self.table = self.table.append(pd.DataFrame(row, index = [len(self.table.index)]))
        #self.table.reset_index()

class PopulationTable(Table):
    def __init__(self):
        self.columns = ['metadata']
        super().__init__(columns = self.columns)
        self.name = "Population Table"

    def add(self, metadata):
        row = {'metadata': metadata}
        self.table = self.table.append(pd.DataFrame(row, index = [len(self.table.index)]))
        #self.table.reset_index()



