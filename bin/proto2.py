#! /usr/bin/env python

import random
import sys
import re
from cStringIO import StringIO
from optparse import OptionGroup
from optparse import OptionParser
import logging

import dendropy
from dendropy import treemanip
from dendropy.utility import probability

_LOGGING_LEVEL_ENVAR = "ARCHIPELAGO_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR = "ARCHIPELAGO_LOGGING_FORMAT"

class RunLogger(object):

    _LOGGER_SET = False

    def __init__(self, **kwargs):
        if not RunLogger._LOGGER_SET:
            self.name = kwargs.get("name", "RunLog")
            self._log = logging.getLogger(self.name)
            self._log.setLevel(logging.INFO)
            if kwargs.get("log_to_stderr", True):
                ch1 = logging.StreamHandler()
                stderr_logging_level = self.get_logging_level(kwargs.get("stderr_logging_level", logging.INFO))
                ch1.setLevel(stderr_logging_level)
                ch1.setFormatter(self.get_default_formatter())
                self._log.addHandler(ch1)
            if kwargs.get("log_to_file", True):
                log_stream = kwargs.get("log_stream", \
                    open(kwargs.get("log_path", self.name + ".log"), "w"))
                ch2 = logging.StreamHandler(log_stream)
                file_logging_level = self.get_logging_level(kwargs.get("file_logging_level", logging.DEBUG))
                ch2.setLevel(file_logging_level)
                ch2.setFormatter(self.get_default_formatter())
                self._log.addHandler(ch2)
            RunLogger._LOGGER_SET = True

    def get_logging_level(self, level=None):
        if level in [logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING,
            logging.ERROR, logging.CRITICAL]:
            return level
        elif level is not None:
            level_name = str(level).upper()
        elif _LOGGING_LEVEL_ENVAR in os.environ:
            level_name = os.environ[_LOGGING_LEVEL_ENVAR].upper()
        else:
            level_name = "NOTSET"
        if level_name == "NOTSET":
            level = logging.NOTSET
        elif level_name == "DEBUG":
            level = logging.DEBUG
        elif level_name == "INFO":
            level = logging.INFO
        elif level_name == "WARNING":
            level = logging.WARNING
        elif level_name == "ERROR":
            level = logging.ERROR
        elif level_name == "CRITICAL":
            level = logging.CRITICAL
        else:
            level = logging.NOTSET
        return level

    def get_default_formatter(self):
        f = logging.Formatter("[%(asctime)s] %(message)s")
        f.datefmt='%Y-%m-%d %H:%M:%S'
        return f

    def get_rich_formatter(self):
        f = logging.Formatter("[%(asctime)s] %(filename)s (%(lineno)d): %(levelname) 8s: %(message)s")
        f.datefmt='%Y-%m-%d %H:%M:%S'
        return f

    def get_simple_formatter(self):
        return logging.Formatter("%(levelname) 8s: %(message)s")

    def get_raw_formatter(self):
        return logging.Formatter("%(message)s")

    def get_logging_formatter(self, format=None):
        if format is not None:
            format = format.upper()
        elif _LOGGING_FORMAT_ENVAR in os.environ:
            format = os.environ[_LOGGING_FORMAT_ENVAR].upper()
        if format == "RICH":
            logging_formatter = self.get_rich_formatter()
        elif format == "SIMPLE":
            logging_formatter = self.get_simple_formatter()
        elif format == "NONE":
            logging_formatter = self.get_raw_formatter()
        else:
            logging_formatter = self.get_default_formatter()
        if logging_formatter is not None:
            logging_formatter.datefmt='%H:%M:%S'

    def debug(self, msg, *args, **kwargs):
        self._log.debug(msg, *args, **kwargs)

    def info(self, msg, *args, **kwargs):
        self._log.info(msg, *args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        self._log.warning(msg, *args, **kwargs)

    def error(self, msg, *args, **kwargs):
        self._log.error(msg, *args, **kwargs)

    def critical(self, msg, *args, **kwargs):
        self._log.critical(msg, *args, **kwargs)

class TotalExtinctionException(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class Region(object):
    """
    An atomic biogeographical unit.
    """

    def matrix_from_file(f):
        if isinstance(f, str):
            stream = open(f, "rU")
        else:
            stream = f
        row_strings = [r.strip() for r in f.readlines() if r.strip()]
        region_names = [r.strip() for r in re.split("[\t ]+", row_strings[0]) if r.strip()]
        rows = [[c.strip() for c in re.split("[\t ]+", r) if c.strip()] for r in row_strings[1:]]
        regions = []
        region_connections = {}
        for i, row in enumerate(rows):
            if row[0] != region_names[i]:
                raise ValueError("Region names must match in column order: expecting '%s' in in first column of row %d, but found '%s'" \
                    % (region_names[i], i+1, row[0]))
            region = Region(label=region_names[i])
            regions.append(region)
            region_connections[region] = []
            for j, c in enumerate(row[1:]):
                c = c.strip()
                try:
                    region_connections[region].append(float(c))
                except ValueError:
                    if c == '-':
                        region_connections[region].append('-')
                    else:
                        raise ValueError("Invalid value for migration probability in row %d, column %d: '%s'" % (i+1, j+1, c))
#            region_connections[region] = [float(c) for c in row[1:]]
            if len(region_connections[region]) != len(region_names):
                raise ValueError("Expecting %d columns in row %d, but found %d" % (len(region_names), i+1, len(region_connections[region])))
            c = region_connections[region].count('-')
            if c == 1:
                region_connections[region][region_connections[region].index('-')] = 1.0 - sum([x for x in region_connections[region] if x != '-'])
            elif c > 1:
                raise ValueError("Multiple columns with '-' specified in row %d" % (i+1))
            if abs(1.0 - sum(region_connections[region])) > 0.001:
                raise ValueError("Probabilities for row %d do not sum to 1.0: %s" % (i+1, [str(p) for p in region_connections[region]]))
        for i, region1 in enumerate(regions):
            for j, region2 in enumerate(regions):
                region1.add_connection(region2, region_connections[region1][j])
#        for region1 in regions:
#            print
#            print region1.label
#            print "%s" % ("   ".join([("%s: %s" % (r.label, region1.connections[r])) for r in regions]))
        return regions
    matrix_from_file = staticmethod(matrix_from_file)

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
        self.run_title = kwargs.get("run_title", "ArchipelagoRun")
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
        self.incidence_log = kwargs.get("incidence_log", \
                open(self.output_prefix + "." + self.run_title + ".incidences.txt", "w"))
        self.diversity_stacked_log = kwargs.get("diversity_stacked_log", open(self.output_prefix + ".summary.stacked.txt", "w"))
        self.diversity_unstacked_log = kwargs.get("diversity_unstacked_log", open(self.output_prefix + ".summary.unstacked.txt", "w"))
        self.tree_log = kwargs.get("tree_log", open(self.output_prefix + ".summary.trees", "w"))
        self.logger = kwargs.get("run_logger", RunLogger(log_path=self.output_prefix + "." + self.run_title + ".log"))
        self.regions = Region.matrix_from_file(kwargs.get("regions_file", self.default_regions()))
        self.tree = None
        self.bootstrapped = False
        self.diversify = None

    def default_regions(self):
        r = """
            A           B           C           D           E

        A   -           %(m)s       0.0         0.0         0.0
        B   %(m)s       -           %(m)s       0.0         0.0
        C   0.0         %(m)s       -           %(m)s       0.0
        D   0.0         0.0         %(m)s       -           %(m)s
        E   0.0         0.0         0.0         %(m)s       -
        """ % {"m" : self.migration_rate}
        return StringIO(r)

    def num_lineages_in_region(self, region):
        return sum([1 for nd in self.tree.leaf_iter() if region in nd.regions])

    def _get_num_lineages(self):
        return len(self.tree.leaf_nodes())
    num_lineages = property(_get_num_lineages)

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
                    lineage.regions.discard(region)
                    if not lineage.regions:
                        if lineage is self.tree.seed_node:
                            raise TotalExtinctionException()
                        else:
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
                for nd in leaf_nodes:
                    nd.edge.length += 1
                ngen += 1
        except KeyboardInterrupt:
            self.logger.warning("%s: Keyboard interrupt: terminating" % self.run_title)
            return False
        except TotalExtinctionException:
            self.logger.warning("%s: All lineages extinct: terminating" % self.run_title)
            return False
        self.logger.info("%s: Completed run after %d generations, with %d lineages in system" % (self.run_title, ngen, len(leaf_nodes)))
        return True

    def report(self, write_summary_headers=True):
        pass
#        self.incidence_log
#        self.diversity_stacked_log
#        self.diversity_unstacked_log
#        self.tree_log

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

    parser.add_option('-X', '--max-lineages',
        action='store',
        dest='max_lineages',
        type='int',
        default=30,
        metavar='NUM-TAXA',
        help="end each simulation replicate when this number of lineages are reached (default=%default)")

    parser.add_option('-n', '--num-reps',
        action='store',
        dest='num_reps',
        type='int',
        default=10,
        metavar='NUM-REPS',
        help="number of simulation replicates (default=%default)")

    parser.add_option('-z', '--random-seed',
        action='store',
        dest='random_seed',
        type='int',
        default=None,
        metavar='SEED',
        help="random number seed (default=NONE)")

    parser.add_option('-o', '--output-prefix',
        action='store',
        dest='output_prefix',
        type='string', # also 'float', 'string' etc.
        default='archipelago_results',
        metavar='OUTPUT-FILE-PREFIX',
        help="prefix for output files (default='%default')")

    (opts, args) = parser.parse_args()

    logger = RunLogger(log_path=opts.output_prefix + ".log")
    rng = random.Random(opts.random_seed)
    rep = 0
    while rep < opts.num_reps:
        arch = Archipelago(
                birth_rate=opts.birth_rate,
                death_rate=opts.death_rate,
                migration_rate=opts.migration_rate,
                max_lineages=opts.max_lineages,
                run_title="R%03d" % (rep+1),
                diversity_stacked_log=open(opts.output_prefix + ".summary.stacked.txt", "w"),
                diversity_unstacked_log=open(opts.output_prefix + ".summary.unstacked.txt", "w"),
                tree_log=open(opts.output_prefix + ".summary.trees", "w"),
                run_logger=logger,
                rng=rng)
        success = arch.run()
        if arch.num_lineages < opts.max_lineages:
            logger.warning("Failed to reach target diversity: re-running replicate %d ('%s')" % (rep+1, arch.run_title))
        elif not success:
            logger.warning("Total extinction of all lineages: re-running replicate %d ('%s')" % (rep+1, arch.run_title))
        else:
            rep += 1

if __name__ == '__main__':
    main()
