#! /usr/bin/env python

import random
import sys
from optparse import OptionGroup
from optparse import OptionParser

import dendropy
from dendropy import treemanip
from dendropy.utility import probability

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

class Archipelago(object):

    SYMPATRIC_RANGE_INHERITANCE = 0
    VICARIANT_RANGE_INHERITANCE = 1
    PARTITION_RANGE_INHERITANCE = 2
    LOCAL_DIVERSIFICATION = 0
    GLOBAL_DIVERSIFICATION = 1

    def __init__(self, *args, **kwargs):
        self.title = kwargs.get("title", "Archipelago Run")
        self.rng = kwargs.get("rng", random.Random())
        self.birth_rate = kwargs.get("birth_rate", 0.2)
        self.death_rate = kwargs.get("death_rate", 0.0)
        self.migration_rate = kwargs.get("migration_rate", 0.1)
        self.diversification_mode = kwargs.get("diversification_mode", \
                Archipelago.LOCAL_DIVERSIFICATION)
        self.range_inheritance = kwargs.get("range_inheritance", \
                Archipelago.VICARIANT_RANGE_INHERITANCE)
        self.max_gens = kwargs.get("max_gens", 10000)
        self.max_lineages = kwargs.get("max_lineages", 30)
        self.output_prefix = kwargs.get("output_prefix", "archipelago_run")
        self.regions = []
        self.create_regions()
        self.tree = None
        self.bootstrapped = False
        self.run_log = sys.stderr
        self.diversify = None

    def create_regions(self):
        self.regions = []
        for label in ["A", "B", "C", "D", "E"]:
            self.regions.append(Region(label, rng=self.rng))
        # A <-> B <-> C <-> D <-> E
        self.regions[0].add_connection(self.regions[1], self.migration_rate)
        self.regions[1].add_connection(self.regions[0], self.migration_rate)
        self.regions[1].add_connection(self.regions[2], self.migration_rate)
        self.regions[2].add_connection(self.regions[1], self.migration_rate)
        self.regions[2].add_connection(self.regions[3], self.migration_rate)
        self.regions[3].add_connection(self.regions[2], self.migration_rate)
        self.regions[3].add_connection(self.regions[4], self.migration_rate)
        self.regions[4].add_connection(self.regions[3], self.migration_rate)

    def num_lineages_in_region(self, region):
        return sum([1 for nd in self.tree.leaf_iter() if region in nd.regions])

    def region_lineages(self):
        region_lineage_map = {}
        for leaf in self.tree.leaf_iter():
            for region in leaf.regions:
                if region not in region_lineage_map:
                    region_lineage_map[region] = [leaf]
                else:
                    region_lineage_map[region].append(leaf)
        return region_lineage_map

    def bootstrap(self, initial_region=None):
        if initial_region is None:
            initial_region = self.rng.choice(self.regions)
        self.tree = dendropy.Tree()
        self.tree.seed_node.regions = set([initial_region])
        self.tree.seed_node.edge.length = 1
        self.bootstrapped = True
        if self.diversification_mode == Archipelago.LOCAL_DIVERSIFICATION:
            self.diversify = self.local_diversification
        else:
            self.diversify = self.global_diversification

    def migrate(self):
        for leaf in self.tree.leaf_iter():
            to_be_added = []
            for region in leaf.regions:
                dest = region.get_migration_destination()
                if dest is not region:
                    if dest not in leaf.regions:
                        to_be_added.append(dest)
            for region in to_be_added:
                if region not in leaf.regions \
                        and (region.carrying_capacity > 0 \
                            and region.carrying_capacity > self.num_lineages_in_region(region)):
                    leaf.regions.add(region)

    def local_diversification(self):
        region_lineage_map = self.region_lineages()
        for region in region_lineage_map:
            for lineage in region_lineage_map[region]:
                u = self.rng.uniform(0, 1)
                if u < self.birth_rate:
                    child1 = lineage.new_child(edge_length=0)
                    child2 = lineage.new_child(edge_length=0)
                    if self.range_inheritance == Archipelago.SYMPATRIC_RANGE_INHERITANCE:
                        child1.regions = set(lineage.regions)
                        child2.regions = set(lineage.regions)
                    elif self.range_inheritance == Archipelago.VICARIANT_RANGE_INHERITANCE:
                        child1.regions = set([self.rng.choice(list(lineage.regions))])
                        child2.regions = set([r for r in lineage.regions if r not in child1.regions])
                    elif self.range_inheritance == Archipelago.PARTITION_RANGE_INHERITANCE:
                        child1.regions = set([self.rng.sample(list(lineage.regions), len(lineage.regions))])
                        child2.regions = set([r for r in lineage.regions if r not in child1.regions])
                    else:
                        raise ValueError("Invalid value for range inheritance mode: %s" % self.range_inheritance)
                elif u > self.birth_rate and u < (self.birth_rate + self.death_rate):
                    treemanip.prune_subtree(self.tree, lineage)

    def global_diversification(self):
        raise NotImplementedError()

    def run(self):
        if not self.bootstrapped:
            self.bootstrap()
        try:
            ngen = 0
            leaf_nodes = self.tree.leaf_nodes()
            while ngen < self.max_gens and len(leaf_nodes) < self.max_lineages:
                self.migrate()
                self.diversify()
                leaf_nodes = self.tree.leaf_nodes()
                if len(leaf_nodes) == 0:
                    self.run_log.write("All lineages extinct: terminating\n")
                    break
                for nd in leaf_nodes:
                    nd.edge.length += 1
                ngen += 1
        except KeyboardInterrupt:
            sys.stderr.write("Terminating on keyboard interrupt.\n")

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

    arch = Archipelago(
            birth_rate=opts.birth_rate,
            death_rate=opts.death_rate,
            migration_rate=opts.migration_rate,
            max_lineages=30,
            rng=rng)
    arch.run()

if __name__ == '__main__':
    main()
