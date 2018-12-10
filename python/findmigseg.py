import gzip
import sys
from optparse import OptionParser


class Event:

    def process_l(self, l):
        o = []
        for i in range(0, len(l)):
            if(l[i] == '1'): o.append(i)
        return o


    def __init__(self, e,p,T,f,t,l):
        self.e = [e]
        self.p = [p]
        self.T = [T]
        self.l = []
        self.l.append(self.process_l( l))
        self.f = [f]
        self.t = [t]
        self.nm = 0

    def update(self, e,p,T,f,t,l):
        self.e.append(e)
        self.p.append(p)
        self.T.append(T)
        self.f.append(f)
        self.t.append(t)
        self.l.append(self.process_l(l)) 
        if(e == "M"): self.nm += 1

    def cout(self):
        print "Printing current event:"
        for i in range(0, len(self.e)): 
           print "\t", self.e[i], self.p[i], self.T[i], self.f[i], self.t[i], self.l
        print("\n")

def print_h (h) : 
    for i in range(len(h)) : 
        print h[i]
   
def make_bed (h, pos) : 
    for j in range(4) : 
        start =0 
        end = 0
        for i in range(1, len(h)) : 
            if( h[i][j] - h[i-1][j] == 1) :
                start = i
            elif ( h[i][j] - h[i-1][j] == -1) :
                end = i

            if ( start > 0 and end > 0): 
                print j, float(pos[start-1]) ,float(pos[end-1]), float(pos[start-1]) - float(pos[end-1])
                start = 0
                end = 0


parser = OptionParser()
parser.add_option("-f",
                  "--from",
                  dest="fp",
                  action="store",
                  help="From population X",
                  type="string")
parser.add_option("-t",
                  "--to",
                  dest="tp",
                  action="store",
                  help="To population Y",
                  type="string")
parser.add_option("-i",
                  "--input",
                  action="store",
                  dest="tree_path",
                  help="Path to gzipped tree file",
                  type="string")
parser.add_option("-s",
                  "--start",
                  dest="start_years",
                  action="store",
                  help="Upper end of era (i.e. 70 kya).",
                  type="float")
parser.add_option("-e",
                  "--end",
                  dest="end_years",
                  action="store",
                  help="Lower end of era (i.e. 50 kya).",
                  type="float")

(options, args) = parser.parse_args()

tree = gzip.open(options.tree_path, 'r')

# Start and end in generations
start = options.start_years / 29
end = options.end_years / 29

fp = options.fp
tp = options.tp

# 2D Array of haplotype states
h = []
h_init = [0,0,0,0]
h.append(h_init[:])

# List of events, positions (for convenience), and the times of the 
# last events for comparison with R events
events = []
pos = []
times = [0,0,0,0]
# So the rest of the loop works

# Read in the first line and initialize an object
line = tree.readline()
e,p,T,f,t,l = line.rstrip().split("\t") 
c = Event(e,p,T,f,t,l)

for line in tree:
    e,p,T,f,t,l = line.rstrip().split("\t")    

    if(e == "R"): 
        events.append(c)
        pos.append(c.p[0])
        c = Event(e,p,T,f,t,l)
    else:
        c.update(e,p,T,f,t,l)

for c in events : 
    upd = h[-1] [:]
    if( len(c.T) == 1 ): 
        h.append(upd[:])
        continue 
  
    # For each recombination event
    #   blank the record
    #   but only if it is BEFORE the last event (otherwise the event is still on the tree) 
    for j in c.l[0] :
        if c.T[0] < times [j] and c.T[1] > times [j]:
            upd [j] = 0
            times [ j ] = 0
        # Then check if there's a migration event that is
        #   i. confirming it (so keep it at 1) 
        for i in range(len(c.e)) : 
            if (c.e[i] == "M"):  
               # print "MIGRATION"
                if ( (float(c.T[i]) < start) and (float(c.T[i]) > end) and (c.t[i] == tp) ) :
                        upd [j] = 1
                        times [j] = c.T[i]
                elif ( ( float(c.T[i]) < start) and ( float(c.T[i]) > 350) and (c.t[i] == fp) ) : 
                        upd [j] = 0                
                        times [j] = 0

    h.append(upd[:])

make_bed (h, pos)
