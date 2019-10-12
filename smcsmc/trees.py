#import smcsmc
import gzip

#class NodeTable:
#    def __init__(self)

#class EdgeTable: 
#    def __init__(self)

class Event:
    def __init__(self, tokens):
        self.pos = tokens[1]
        self.time = tokens[2]
        self.lineages = self.parse_lineages(tokens[5])
        self.migration_events = []
        self.nmig = 0

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
        print("Event at " + self.pos + " in " + ''.join([str(int(s)) for s in self.lineages])+"\n" +
               "-\tTime: " + self.time + '\n' +
               '-\tCoal Time: ' + self.coal + '\n' +
               '-\tN Migration: ' + str(self.nmig))
        for e in self.migration_events:
            e.print()


class Migration(Event):
    def __init__(self, tokens):
        assert(tokens[0] == "M")
        super().__init__(tokens) 
        self.To = tokens[3]
        self.From = tokens[4]
        self.migration_events = None

    def print(self):
        print("*\t[" + ''.join([str(s) for s in self.lineages]) + "] \tMigration " + self.From + " -> " + self.To + " @ " + str(self.time) + '\n')

        
class Record:i
    def __init__(self, tree):
        self.events = {}

        with gzip.open(tree, 'rb') as f:
            current_event = Event([s.decode("utf-8") for s in f.readline().strip().split()])
            for line in f:
                tokens = [s.decode("utf-8") for s in line.strip().split(b'\t')]
                if tokens[0] == "R":
                    #import pdb; pdb.set_trace()
                    self.addEvent(current_event)
                    current_event = Event(tokens)
                else:
                    current_event.append(tokens)

    def addEvent(self, event):
        self.events[event.pos] = event


