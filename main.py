import csv
import argparse
import pdb

###DBG
import sys
sys.argv = ['main.py','data.csv']

_NAME   = 'Antibody'
_WEIGHT = 'Molecular weight (kD)'
_HOST   = 'HOST'
# currently hosts are treated as green={mouse,goat} red={rabbit}
hostMap = {'Mouse':0,'Goat':0,'Rabbit':1}

class ProteinInfo(object):
    _WEIGHT_MARGIN = 0
    def __init__(self, name, weight, host):
        self.name = name
        self.weight = weight
        self.host = host

    def hasWeightCollisionWith(self, anotherProtein):
        margin = ProteinInfo._WEIGHT_MARGIN
        myCeiling = self.maxWeight() + margin
        myFloor   = self.maxWeight() - margin
        otherCeiling = anotherProtein.maxWeight() + margin
        otherFloor   = anotherProtein.maxWeight() - margin
        
##        print self.name
##        print anotherProtein.name
##        print "me:  %f ---- %f" % (myCeiling, myFloor)
##        print "him: %f ---- %f" % (otherCeiling, otherFloor)
##        print (myCeiling > anotherProtein.getWeight() > myFloor)
##        print (myCeiling > otherCeiling > myFloor)
##        print (myCeiling > otherFloor > myFloor)
        
        return (myCeiling > anotherProtein.maxWeight() > myFloor) \
                       or (myCeiling > otherCeiling > myFloor) \
                       or (myCeiling > otherFloor > myFloor)
        
    def hasHostCollisionWith(self, anotherProtein):
        return hostMap[self.getHost()] == hostMap[anotherProtein.getHost()]
        
    def collidesWith(self, anotherProtein):
        return self.hasWeightCollisionWith(anotherProtein) and self.hasHostCollisionWith(anotherProtein)

    def getName(self):
        return self.name
                
    def getWeight(self):
        return self.weight

    def maxWeight(self):
        return max(self.weight) if self.weight else float('-inf')

    def getHost(self):
        return self.host

    def __str__(self):
        return "%s (%s kDa), %s" % (self.name,
                                    ' '.join([str(v) for v in self.weight]),
                                    self.host)

def initArgs():
    parser = argparse.ArgumentParser(description="Sort some proteins")
    parser.add_argument(
                'csv',
                metavar='csv', 
        type=argparse.FileType('r'),
                help='Input list of proteins as a .csv')
    parser.add_argument(
                '-m',
                metavar='m',
                type=float, help='Atomic weight (kDa) margin',
                default=7)
    return parser.parse_args()

def initProteins(reader, weightMargin):
    ProteinInfo._WEIGHT_MARGIN = weightMargin
    proteins = []
    for line in reader:
        proteins.append(
            ProteinInfo(
                line[_NAME],
                [float(v) for v in line[_WEIGHT].split(',') if v],
                line[_HOST]))
    return proteins
    
def binProteins(proteins):
    '''
    Given a list of ProtenInfo objects, return them sorted
    into bins so that there are no collisions (based on
    weightMargin and color)
    
    currently doing naiive O(n^2) just to get up and running
    '''
    # list that is populated with sets of proteins
    bins = []
    for protein in proteins:
        assigned = False
        for bin in bins:
            collisions = False
            for placedProtein in bin:
                if placedProtein.collidesWith(protein):
                    collisions = True
                    break
            if not collisions:
                bin.add(protein)
                assigned = True
                break
        if not assigned:
            bins.append(set([protein]))
    return bins
    

def run(filepath, weightMargin):
    proteins = initProteins(csv.DictReader(filepath), weightMargin)
    bins     = binProteins(proteins)
    for b in bins:
        for p in b:
            print p
        print '-------------'

if __name__ == "__main__":
    args = initArgs()
    run(args.csv, args.m)
