from __future__ import print_function

import sys
import gzip
from collections import namedtuple
import pdb
import logging
import pandas as pd

##
## To do:
## - end position of edges of initial tree do not always get updated correctly
## - smcsmc: initial events are not always recorded, leading sometimes to inconsistent migration histories at beginning
## - include samples (at height 0)?  Need to know start heights and populations
## - check that migration table conforms to tskit requirements
##


#
# define objects, constants, and globals
#
Event = namedtuple('Event', 'height type pos frm to descendants nodeid')
Node = namedtuple('Node', 'flags time population individual metadata')
Edge = namedtuple('Edge', 'left right parent child')
Migration = namedtuple('Migration', 'left right node source dest time descendants')
Population = namedtuple('Population', 'metadata')

class migrationSegment:
    def __init__(self, event):
        self.event = event
        self.start = event.pos

    def stop(self, pos):
        self.stoppos = pos
        self.complete = True

    def to_event(self):
        assert(self.complete)
        return(Migration(self.start, self.stoppos, self.event.nodeid, self.event.frm, self.event.to, self.event.height, self.event.descendants))

RECOMBINATION_ANCESTOR = 2 ** 62

STORE_RECOMBINATION = True

nodelist = []

#
# build an event from a .trees line
#
def makeEvent( datum ):
    global nodelist
    type, pos, height, frm, to, desc = datum
    event = Event( float(height), type, float(pos), int(frm), int(to), int(''.join(reversed(desc)),2), len(nodelist) )
    metadata = type if type != "R" else float(pos)
    nodelist.append( Node( 0, float(height), int(frm), -1, metadata) )
    return event

#
# modify descendants
#
def setDescendants( event, descendants ):
    return Event( event.height, event.type, event.pos, event.frm, event.to, descendants, event.nodeid )

#
# modify right endpoint
#
def setRight( edge, right ):
    return Edge( edge.left, right, edge.parent, edge.child )

#
# modify population - for recombination
#
def setPopulation( node, population ):
    return Node( node.flags, node.time, population, node.individual, node.metadata )

#
# for printing
#
def eventRepr( event ):
    return "{:7.1f} {} {:10.1f} {: 1}->{: 1}  {}".format(event.height, event.type, float(event.pos), event.frm, event.to, ''.join(reversed("{:b}".format(event.descendants))))

def treeRepr( tree ):
    return "\n  ".join( [""] + [eventRepr(event) for event in tree ] )

#
# update a current tree with a new event
#
def update( tree, event, start_height, lineage ):
    # 'tree' is a list of Events defining a tree, ordered by height
    #print("Update: start=",start_height," lineage={:b}".format(lineage))
    
    # find first event beyond start_height
    idx = 0
    while idx < len(tree) and tree[idx].height <= start_height:
        idx += 1
    
    # trace ancestors of recombining lineage, and remove descendants from that line.
    # keep nodes, which now may have no descendants, and mark them as ancestors of recombination
    # if a tree currently has a loose branch, no ancestors will be present
    #pdb.set_trace()
    for anc_idx in range( idx, len(tree) ):
        if lineage & tree[anc_idx].descendants:
            check( lineage & tree[anc_idx].descendants == lineage, "Inconsistent tree structure" )
            tree[anc_idx] = setDescendants( tree[anc_idx], (tree[anc_idx].descendants & ~lineage) | RECOMBINATION_ANCESTOR )
            
    # find first event beyond event height
    while idx < len(tree) and tree[idx].height < event.height:
        idx += 1
    check( idx >= len(tree) or tree[idx].height > event.height, "duplicate heights found" )

    # insert event
    tree = tree[:idx] + [event] + tree[idx:]
    if event.type != 'C':
        return tree

    # for a coalescence, add descendants to ancestors of lineage that we coalesced with.
    # if coalescence was into own lineage, use recombination_ancestor mark to identify those nodes
    global back_coals
    coalescence_lineage = event.descendants & ~lineage
    if coalescence_lineage == 0:
        coalescence_lineage = RECOMBINATION_ANCESTOR
        back_coals['current'] = True
        back_coals['nBC'] += 1
        back_coals['time'] = event.height
        back_coals['descendants'] = event.descendants
    else:
        back_coals['current'] = False
        back_coals['nFC'] += 1

    idx += 1

    #if coalescence_lineage & RECOMBINATION_ANCESTOR:
      #  pdb.set_trace()

    while idx < len(tree):
        if tree[idx].descendants & coalescence_lineage:
            tree[idx] = setDescendants( tree[idx], tree[idx].descendants | lineage )
        tree[idx] = setDescendants( tree[idx], tree[idx].descendants & ~RECOMBINATION_ANCESTOR )
        idx += 1
    return tree


#
# normalize tree - remove recombination ancestor mark and superfluous "coalescences"
#
def normalize( tree ):
    norm_tree = []
    maxlineage = max( event.descendants for event in tree )
    lineages = [ 2**i for i in range(len(bin(maxlineage))-2) ]  ## all tips
    for event in tree:
        event = setDescendants( event, event.descendants & ~RECOMBINATION_ANCESTOR )
        if event.descendants == 0: 
            continue
        if event.type == "C":
            null_coalescence = [True
                                for lineage in lineages
                                if lineage & event.descendants == event.descendants ]
            if null_coalescence != []:
                continue
            lineages.append( event.descendants )
        norm_tree.append(event)
    return norm_tree

#
# update edges - enter rightmost sequence position into disappeared edges, and initialize new ones
#       now also updates a migration table in the same way as the edgelist, and outputs a list 
#       of "completed" (i.e. having both start and end coordinates) migration events.
#
def update_edges_and_migrations( current_edges, tree, edgelist, current_pos, current_migrations, segments, hap, d = False):
    if tree == []:
        maxlineage = 1
    else:
        maxlineage = max( event.descendants for event in tree )
    new_edges = {}
    new_migrations = {}
    ## map from lineage (encoded by descendants) to current event on that lineage
    ## tips correspond to null events; should represent explicitly
    lineages = {2**i : None for i in range(len(bin(maxlineage))-2) }
    for event in tree:
        parent = event.descendants
        children = [ lineage for lineage in lineages.keys() if lineage & parent ]
        check( children != [], "No children found for parent" )
        for child in children:
            edge = (event.nodeid, lineages[child])
            del lineages[child]
            if edge[1] is None: 
                continue
            if edge not in current_edges:
                new_edges[ edge ] = len(edgelist)
                edgelist.append( Edge( current_pos, current_pos, edge[0], edge[1] ) )
            else:
                new_edges[ edge ] = current_edges[ edge ]
                del current_edges[ edge ]
        lineages[parent] = event.nodeid
    ## enter end sequence position for edges that are no longer current: 
    for edge in current_edges:          
        edgelist[ current_edges[ edge ] ] = setRight( edgelist[ current_edges[ edge ] ], current_pos )
    ## done
     
    ## Add any new migration events to the list
    migrations = [e for e in tree if e.type == 'M' and (e.descendants & hap)] 


    global r_count
    global back_coals
    all_migrations = [ e for e in tree if e.type == 'M' ]
    should_be_deleted = []
    ## find migrations which have changed descendants and check if they are 
    ## removed from the list. 
    for e in all_migrations:
        key = (e.pos, e.height)
        if key in current_migrations:
            if current_migrations[key].event.descendants & hap:
                if e.descendants & hap != hap:
                    #print("descendants present in previous tree but not in the one at " + str(e.pos) )
                    r_count += 1
                    should_be_deleted.append((e.pos, e.height))



    ## This is all just debug
    for m in migrations: 
        key = (m.pos, m.height)
        #pdb.set_trace()
        if key in current_migrations:
            if not event_in_BC(m):
                new_migrations[key] = migrationSegment(m)
                del current_migrations[key]
            #else:
            #    pdb.set_trace()
            #else:
                #print("killed a semgnet")
        else: 
            ## This line is changing the behaviour to examine the behavior if BM segments should be includded:
            ## remove all the if statements if you want to revert to the original behaviour
            new_migrations[key] = migrationSegment(m)
 
   

    ## Now we have a continueing list of the segments in the 
    ## new_migrations object.
    ## 
    ## Now deal with the ones that are ending 
    global back_coals
    if d: print("ending " + str(len(current_migrations.keys())) + " segments")
    for m in current_migrations.values():
        m.stop(current_pos)
        mig = m.to_event()
        key = (mig.left, mig.time)
        if key not in segments:
            segments[key] = m.to_event()
            if back_coals['current']:
                back_coals['BC'] += 1
            else:
                back_coals['FC'] += 1
           
    return (new_edges, new_migrations, segments)


def event_in_BC(m):
    global back_coals
    if back_coals['current']:
        if m.height < back_coals['time'] and m.height > back_coals['rtime']:
            if m.descendants & back_coals['descendants']:
                return True
    return False

def is_migration_node(node):
    ## Assuming data is global, which it seems to be
    events = { datum.nodeid : datum for datum in data }
    return (events[node].type == 'M')

#
# check various consistencies
#
def check( isok, message ):
    'turn off the printing for now'
    #if not isok:
        #print(message,"\n")


Output = namedtuple("Output", 'nodelist, edgelist, migrationlist')


# Globals for recording various things
r_count = 0
back_coals = {'nFC': 0, 'nBC': 0,'BC': 0, 'FC': 0, 'current': False, 'time': 0, 'rtime': 0, 'descendants':0}


#
# process data from bottom (leftmost) up
#
def trees2tskit(infile, hap = 2, suffix = '.trees.gz', write_all = False, d = False):
    prefix = infile[ :-len(suffix) ]
    data = gzip.open(infile).readlines()
    data = [ makeEvent(line.decode('utf-8').strip().split()) for line in reversed(data) ]


    idx = 0
    tree = []
    current_edges = {}
    edgelist = []
    segments = {}
    migration_list = {} 

    while idx < len(data):
        # idx points to the first of the events associated to the current recombination
        # first identify the corresponding recombination event
        # this may not exist, in the case of the initial tree
        end_idx = idx + 1
        while (end_idx < len(data) and 
               data[end_idx].type != 'R' and 
               data[end_idx - 1].type != 'C' and 
               data[end_idx].pos == data[end_idx - 1].pos):
            end_idx = end_idx + 1
            
        # if we have a recombination event, process that first.
        # if not, assume we are initializing a tree
        if end_idx < len(data) and data[end_idx].type == 'R':
            #print ("Recombination ",data[end_idx])
            if STORE_RECOMBINATION:
                ## enter population information into recombination node, using first event in block
                nodelist[ data[end_idx].nodeid ] = setPopulation( nodelist[ data[end_idx].nodeid ], nodelist[ data[idx].nodeid ].population )
            start_height = data[end_idx].height
            lineage = data[end_idx].descendants
            next_idx = end_idx + 1
        else:
            start_height = 0
            lineage = 2 ** (len(bin(data[end_idx - 1].descendants)) - 3)  ## rightmost lineage of those participating in coalescent
            #print("Novel lineage?  Desc=",data[end_idx - 1].descendants," lineage=",lineage)
            next_idx = end_idx

        global back_coals
        back_coals['rtime'] = start_height
            
        # process the events
        for event in data[idx : end_idx]:
            #print ("Processing ",event)
            check( start_height < event.height, "****** Recombination occurred above coalescent or migration")
            check( lineage & event.descendants == lineage, "****** Inconsistent tree structure" )
            tree = update(tree, event, start_height, lineage)
            start_height = event.height
            lineage = event.descendants

        # normalize tree
        tree = normalize( tree )

        # update edges
        if d:
            print(data[next_idx-1].pos)
            print(treeRepr(tree))
            
        
        current_edges, migration_list, segments = update_edges_and_migrations( current_edges, tree, edgelist, data[next_idx - 1].pos, migration_list, segments, hap = hap,d = d)
        
        if d:
            print(len(segments))

        # done
        idx = next_idx

    # enter last position into remaining edges
    
    #update_edges_and_migrations( current_edges, [], edgelist, data[next_idx - 1].pos, migration_list, segments)
    #print(segments)

    # build migration table.  tskit does not expected Nodes to be migration events,
    # and the documentation says that the migration time should be _strictly_ between
    # the parent and child node times.  We break that assumption, since migration events
    # coincide with nodes 
    migrationlist = list(segments.values())
    global r_count
    #print("Total: " + str(r_count))
    r_count = 0

 
    with open('back_coals.txt', 'a') as f:
        pd.DataFrame(back_coals, index = [0]).to_csv(f, header = False, index = False)
 
    back_coals = {'nFC': 0, 'nBC': 0,'BC': 0, 'FC': 1, 'current': False, 'time': 0, 'rtime':0, 'descendants':0}


    #nodeid2event = { datum.nodeid : datum for datum in data }
    #for edge in edgelist:
    #    child = edge.child
    #    event = nodeid2event[ child ]
    #    if event.type == 'M':
    #        migrationlist.append( Migration( edge.left, edge.right, child, event.frm, event.to, event.height ) )

    ## write the output
    # Disable all this for now:
    if write_all:
        fout = open(prefix + ".node","w")
        fout.write("is_sample\ttime\tpopulation\tindividual\tmetadata\n")
        for node in nodelist:
            fout.write("\t".join(map(str,list(node))) + "\n")
        fout.close

        fout = open(prefix + ".edge","w")
        fout.write("left\tright\tparent\tchild\n")
        for edge in edgelist:
            fout.write("\t".join(map(str,list(edge))) + "\n")
        fout.close()

        fout = open(prefix + ".migration","w")
        fout.write("left\tright\tnode\tsource\tdest\ttime\n")
        for migration in migrationlist:
            fout.write("\t".join(map(str,list(migration))) + "\n")
        fout.close()

        populations = max( event.frm for event in data ) + 1
        fout = open(prefix + ".population","w")
        fout.write("metadata\n")
        for pop in range(populations):
            fout.write("0\n")
        fout.close()

        if False:
            for i,n in enumerate(nodelist):
                print (i,n)
            for e in edgelist:
                print (e, e.right-e.left)
            for m in migrationlist:
                print (m)

    return Output(nodelist, edgelist, migrationlist)
