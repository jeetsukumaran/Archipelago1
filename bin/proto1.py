#! /usr/bin/env python

import random
import sys
from optparse import OptionGroup
from optparse import OptionParser

from dendropy.utility import probability


class Lineage(object):
    """
    Tracks a lineage in the simulation.
    """

    def __init__(self, parent=None):
        """
        Instantiates a new Lineage, with parent lineage given by `parent`. If
        `parent` is not given, assumed to be root lineage.
        """
        self.parent = parent
        self.children = []
        self.generations = 0

    def split(self):
        self.children = []
        self.children.append(Lineage(self))
        self.children.append(Lineage(self))
        return self.children

    def __repr__(self):
        return hex(id(self))

class Region(object):
    """
    An atomic biogeographical unit.
    """

    def __init__(self, label=None, carrying_capacity=30, rng=None):
        self.label = label
        if rng is None:
            self.rng = random.Random()
        else:
            self.rng = rng
        self.lineages = []
        self.connections = {}
        self.carrying_capacity = carrying_capacity

    def add_connection(self, dest, prob):
        self.connections[dest] = prob

    def get_migration_destination(self):
        dests = self.connections.keys()
        probs = [self.connections[d] for d in dests]
        k = 1.0 - sum(probs)
        if k >= 0:
            dests.append(self)
            probs.append(k)
        else:
            raise ValueError('Sum of probabilities of connections > 1')
        return probability.lengthed_choice(dests, probs)

class DiversificationModel(object):

    def __init__(self, rng=None):
        if rng is None:
            self.rng = random.Random()
        else:
            self.rng = rng

    def process(self, world):
        raise NotImplementedError()

class LocalBirthDeathDiversificationModel(DiversificationModel):

    def __init__(self,
            birth_rate,
            death_rate,
            rng=None):
        DiversificationModel.__init__(self, rng=rng)
        self.birth_rate = birth_rate
        self.death_rate = death_rate

    def process(self, world):
        for region in world.regions:
            for lineage in region.lineages:
                u = self.rng.uniform(0, 1)
                if u < self.birth_rate:
                    world.split_lineage_in_region(lineage, region=region)
                elif u > self.birth_rate and u < (self.birth_rate + self.death_rate):
                    world.remove_lineage_from_region(lineage, region)

class World(object):

    def __init__(self, title, div_model, output_prefix="archipelago-results"):
        self.title = title
        self.div_model = div_model
        self.rng = self.div_model.rng
        self.regions = []
        self.lineages = []
        self._lineage_regions = None
        self.started = False

    def new_region(self, label):
        self.regions.append(Region(label=label, rng=self.rng))
        return self.regions[-1]

    def get_lineage_regions(self, lineage):
        if not self._lineage_regions:
            self._lineage_regions = {}
            for region in self.regions:
                for lin in region.lineages:
                    if lin not in self._lineage_regions:
                        self._lineage_regions = [lin]
                    else:
                        self._lineage_regions.append(lin)
        if lineage not in self._lineage_regions:
            return []
        else:
            return self._lineage_regions[lineage]

    def split_lineage_in_region(self, lineage, region):
        assert lineage in self.lineages
        sys.stderr.write("Splitting ...\n")
        child1, child2 = lineage.split()
        occurring_regions = [r for r in self.regions if lineage in r.lineages]
        if len(occurring_regions) == 1:
            return
        self.lineages.remove(lineage)
        self.lineages.append(child1)
        self.lineages.append(child2)
        for r in occurring_regions:
            r.lineages.remove(lineage)
            if r is not region:
                r.lineages.append(child1)
            else:
                r.lineages.append(child2)

    def remove_lineage_from_region(self, lineage, region):
        assert lineage in self.lineages
        sys.stderr.write("Removing ...\n")
        if region is not None and lineage in region.lineages:
            region.lineages.remove(lineage)
            for r in self.regions:
                if lineage in r.lineages:
                    return
            self.lineages.remove(lineage)

    def migrate(self):
        for r in self.regions:
            for lineage in r.lineages:
                assert lineage in self.lineages
                dest_region = r.get_migration_destination()
                if dest_region is not r \
                        and lineage not in dest_region.lineages\
                        and len(dest_region.lineages) < dest_region.carrying_capacity:
                    dest_region.lineages.append(lineage)

    def start(self, seed_lineage=None, seed_region=None):
        assert self.regions
        sys.stderr.write("Bootstrapping simulation %s...\n" % self.title)
        if seed_lineage is None:
            seed_lineage = Lineage()
        if seed_region is None:
            seed_region = self.rng.choice(self.regions)
        self.lineages.append(seed_lineage)
        seed_region.lineages.append(seed_lineage)

    def cycle(self, ngen):
        if not self.started:
            self.start()
        try:
            for i in xrange(ngen):
                sys.stderr.write("Run %s, Generation %d/%d ...\n" % (self.title, i+1, ngen))
                self.migrate()
                self.div_model.process(self)
                if len(self.lineages) == 0:
                    sys.stderr.write("All lineages extinct: terminating\n")
                    break
        except KeyboardInterrupt:
            sys.stderr.write("Terminating on keyboard interrupt.\n")

    def write_report_region_diversity_header(output=sys.stdout):
        output.write("Region\tLineages\tEndemic\tTotal\n")
    write_report_region_diversity_header = staticmethod(write_report_region_diversity_header)

    def report_region_diversity(self, output=sys.stdout):
#        print self.lineages
        for r in self.regions:
            num_lineages = len(r.lineages)
            num_endemics = 0
            other_region_lineages = []
            for o in self.regions:
                if o is not r:
                    other_region_lineages.extend(o.lineages)
            for lineage in r.lineages:
                if lineage not in other_region_lineages:
                    num_endemics += 1
            output.write("%s\t%d\t%d\t%d\n" % (r.label, num_lineages, num_endemics, len(self.lineages)))
#            print r.label, r.lineages

def create_world(title, div_model, migration_rate):
    world = World(title, div_model)
    region_A = world.new_region('A')
    region_B = world.new_region('B')
    region_C = world.new_region('C')
    region_A.add_connection(region_B, migration_rate)
    region_B.add_connection(region_A, migration_rate)
    region_B.add_connection(region_C, migration_rate)
    region_C.add_connection(region_B, migration_rate)
    return world

def run1(**kwargs):
    rng = kwargs.get("rng", random.Random())
    div_model = LocalBirthDeathDiversificationModel(
            birth_rate=kwargs.get("birth_rate", 0.1),
            death_rate=kwargs.get("death_rate", 0.1),
            rng=rng)
    output = open(kwargs.get("output_path", "archipelago_results") + ".txt", "w")
    World.write_report_region_diversity_header(output)
    nreps = kwargs.get('nreps', 100)
    ngens_per_rep = kwargs.get('ngens_per_rep', 100)
    rep = 0
    while rep < nreps:
        world = create_world("R%03d" % (rep+1), div_model, kwargs.get("migration_rate", 0.1))
        world.cycle(ngens_per_rep)
        if len(world.lineages) > 0:
            world.report_region_diversity(output)
            rep += 1
        else:
            sys.stderr.write("[ALL LINEAGES EXTINCT: REPEATING SIMULATION]\n")

def main():
    """
    Main CLI handler.
    """

    parser = OptionParser(add_help_option=True)

    parser.add_option('-b', '--birth-rate',
        action='store',
        dest='birth_rate',
        type='float',
        default=0.1,
        metavar='LAMBDA',
        help="probability of speciation (default=%default)")

    parser.add_option('-d', '--death-rate',
        action='store',
        dest='death_rate',
        type='float',
        default=0.1,
        metavar='MU',
        help="probability of extinction (default=%default)")

    parser.add_option('-m', '--migration-rate',
        action='store',
        dest='migration_rate',
        type='float',
        default=0.1,
        metavar='RHO',
        help="probability of migration (default=%default)")

    parser.add_option('-z', '--random-seed',
        action='store',
        dest='random_seed',
        type='int',
        default=42,
        metavar='SEED',
        help="random number seed (default=%default)")

    parser.add_option('-o', '--output-prefix',
        action='store',
        dest='output_prefix',
        type='string', # also 'float', 'string' etc.
        default='archipelago_results',
        metavar='OUTPUT-FILE-PREFIX',
        help="prefix for output files (default='%default')")

    (opts, args) = parser.parse_args()

    rng = random.Random(opts.random_seed)

    run1(birth_rate=opts.birth_rate,
         death_rate=opts.death_rate,
         migration_rate=opts.migration_rate,
         rng=rng)

if __name__ == '__main__':
    main()
