import gzip
from pandas import DataFrame
from collections import OrderedDict
from numpy import unique, argmin, argmax
from copy import deepcopy
import pdb

class NoEdgeFound(Exception):
    pass

class RedundantNode(Exception):
    pass

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
            if self.recipient == "R":
                self.recipient = [0,1,2,3]
            #for x in self.lineages:
            #    try: 
            #        self.recipient.remove(x) # need to get rid of the donating branch.
            #    except AttributeError:
            #        self.recipient = [0,1,2,3] # This reroots the tree
        
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
        time1 = 1000
        time2 = 1000
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

        self.tables.root = 6

    def printTree(self):   
        print("[" + self.right + ","+  self.left + ")")
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
                else:
                    current_event.append(tokens)

                self.addEvent(current_event)
                current_event.left = '1' # in preperation for scaling

    def addEvent(self, event):
        self.events[event.right] =event

    def updateTrees(self):
        first = list(self.events.keys())[0]
        self.events[first].initialiseTables()
        tree = self.events[first].tables
        for k, v in self.events.items():
            self.events[k].tables = tree
            self.events[k].tables.update(self.events[k])
            tree = deepcopy(self.events[k].tables)



     
class TableCollection:
    def __init__(self, dict):
        self.dict = dict

    def update(self, entry):
        # 1. Add a coalescent node 
        recombining_node = self.find_mrc_node_older_than(entry.lineages, entry.time)
        mrc_recipient = self.find_mrc_node_older_than(entry.recipient, entry.coal)

        coal_id = entry.tables.dict['Node'].add(time = entry.coal, population = entry.population, individual = entry.tables.dict['Node'].table.loc[recombining_node]['individual'])
        if coal_id == 13:
            pdb.set_trace()
        # 2. Attach the proper nodes to the new one
        #   Add an edge to represent the R to the C node
                #   Remove the edge from the recombining node to any previous nodes
        #   This potentially leaves redundant nodes... I'm hoping these will be removed by tskit
        try:
            parent = entry.tables.dict['Edge'].remove_parents(entry.right, entry.left, recombining_node)
            if self.root == parent:
                self.root = recombining_node
                entry.tables.dict['Edge'].add(entry.right, entry.left, coal_id, parent)

            #if len(entry.tables.dict['Edge'].table.loc[entry.tables.dict['Edge'].table['parent'] == parent]) > 0:
            #    entry.tables.dict['Edge'].add(entry.right, entry.left, coal_id, parent)
        except NoEdgeFound:
            print("No edge found...?")

        #   Thread into the middle of the recipient node and its parent
        if coal_id != recombining_node:
            entry.tables.dict['Edge'].add(entry.right, entry.left, coal_id, recombining_node)

        entry.tables.dict['Edge'].thread(self.dict['Node'].table, entry.right, entry.left, coal_id, mrc_recipient)

        # Deal with the case where the threaded node is older than the upper join 

        entry.tables.evaluate_reroot(entry.tables.dict['Edge'].table, coal_id, self.root, '0') 
        entry.tables.dict['Edge'].prune()
        
    def parents(self,node): 
        """Find all parents of a node on the path to the root"""
        try: 
            table = self.dict['Edge'].table
        except KeyError:
            raise("You must update with the filled tables before finding parents.") 

        if int(self.dict['Node'].table.loc[node]['individual']) == -1: 
            return ["Root"]
        else: 
            new_node = int(table.loc[table['child'] == node]['parent'])
            return [new_node] + self.parents(new_node)

    def intersect_parents(self, nodes): 
        """Find the parents of a group of nodes, and find their intersection"""
        list_of_parents = [set(self.parents(node)) for node in nodes]
        assert(len(list_of_parents) > 1)
        common_parents = list(list_of_parents[0].intersection(*list_of_parents[1:]))
        common_parents = [node for node in common_parents if node != "Root"]
        return common_parents

    def find_times(self, nodes):
        """Find the times of a group of nodes, minus the root""" 
        return [self.dict['Node'].table.loc[node]['time'] for node in nodes]

    def find_mrc_node(self, nodes):
        """Identify the most recent common shared node for a group of samples."""
        assert(type(nodes) == list)
        if len(nodes) == 1:
            return nodes[0]
        else:
            common_parents = self.intersect_parents(nodes)
            ages = self.find_times(common_parents)
            assert(float(min(ages)) > 0)
            mrc = common_parents[argmin(ages)] 
            return mrc

    def find_mrc_node_older_than(self, nodes, time):
        """Identify the most recent common shared node for a group of samples."""
        assert(type(nodes) == list)
        if len(nodes) == 1:
            return nodes[0]
        else:
            common_parents = self.intersect_parents(nodes)
            ages = self.find_times(common_parents)
            if len(ages) > 1: 
                ages = [float(age) for age in ages if float(age) < float(time)]
            assert(float(max(ages)) > 0)
            mrc = common_parents[argmax(ages)] 
            return mrc


    def evaluate_reroot(self, table, new_node, old_root, id):
        # Do we have a connection between the old root and the new node?
        if len(table.loc[(table['child'] == old_root) & (table['parent'] == new_node)]) != 0:
            if float(self.dict['Node'].table.loc[old_root]['time']) < float(self.dict['Node'].table.loc[new_node]['time']):
                # Then we have a new root.
                self.dict['Node'].table.loc[old_root]['individual'] = id
                self.dict['Node'].table.loc[new_node]['individual'] = -1
                self.root = new_node

class Table: 
    def __init__(self, columns):
        self.table = DataFrame(columns = columns) 
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
        self.table = self.table.append(DataFrame(row, index = [len(self.table.index)]))
        return( len(self.table.index)-1)
        

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
        if parent != child: 
            row = {'left': left, 'right': right, 'parent': parent, 'child': child}
            self.table = self.table.append(DataFrame(row, index = [len(self.table.index)])) 

    def remove(self, left, right, parent, child):
        # I've made the indices a bit confusing here.
        id = self.table.loc[(self.table['parent'] == parent) & (self.table['child'] == child)].index
        if len(id) == 0:
            raise NoEdgeFound()
    
        parent = self.table.loc[id]['parent']       
        self.table.drop([self.table.index[int(id)]], inplace = True)
        return(int(parent))

    def remove_parents(self, left, right, child):
        self.table.reset_index(drop = True, inplace = True)
        id = self.table.loc[self.table['child'] == child].index
        if len(id) == 0:
            raise NoEdgeFound()

        parent = self.table.loc[id]['parent']       
        self.table.drop([self.table.index[id[0]]], inplace = True)
        return(int(parent)) 
        
    def thread(self, table, left, right, new, recipient):
        try:
            #pdb.set_trace()
            parent = self.remove_parents(right, left, recipient)
            if float(table.loc[new]['time']) < float(table.loc[parent]['time']):
                self.add(left, right, parent, new)
        except NoEdgeFound:
            #self.add(right, left, recipient, new)
            print("This must connect directly to the Root")

        self.add(left, right, new, recipient)

    def prune(self):
        # Find internal nodes and reroute them.
        for child in self.table['child'].tolist():
            try:
                parent = self.table.loc[self.table['child'] == child]['parent'].tolist()[0]
                grandparent = self.table.loc[self.table['child'] == parent]['parent'].tolist()[0]

                if len( self.table.loc[self.table['parent'] == parent] ) == 1:
                    if len( self.table.loc[self.table['parent'] == grandparent] ) == 1:
                            parent_idx = self.table.loc[self.table['parent'] == parent].index[0]
                            grandparent_idx = self.table.loc[self.table['parent'] == grandparent].index[0]
                            self.table.drop([parent_idx], inplace = True)
                            self.table.loc[grandparent_idx]['child'] = child
            except IndexError: 
                continue




class PopulationTable(Table):
    def __init__(self):
        self.columns = ['metadata']
        super().__init__(columns = self.columns)
        self.name = "Population Table"

    def add(self, metadata):
        row = {'metadata': metadata}
        self.table = self.table.append(DataFrame(row, index = [len(self.table.index)]))
