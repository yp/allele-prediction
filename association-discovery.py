#!/usr/bin/env python

##########
#
#                               HapAllelePred
#   Haplotype-based prediction of gene alleles using SNP genotypes
#
#  Copyright (C) 2012,2013  Yuri Pirola <yuri.pirola(-at-)gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#
#  This file is part of HapAllelePred.
#
#  HapAllelePred is free software: you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  HapAllelePred is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with HapAllelePred.
#  If not, see <http://www.gnu.org/licenses/>.
#
##########

from __future__ import print_function

INVALID_IDS = ['.', '0']
INVALID_ALLELES = ['-']
HD_GENOTYPE_MAPPING = { '0': '1 1',
                        '1': '1 2',
                        '2': '2 2',
                        '5': '0 0'  }
MISSING = '0 0'
GENDER_MALE = 1
GENDER_FEMALE = 2

# SAT solver: clasp
REHCSTAR_CHECK_ERRORS = \
    './reHCstar -3' + \
    ' --pedigree {infile} --haplotypes {outfile}' + \
    ' --global-error --global-error-number {errors}' + \
    ' --sat-cmdline "nice -n19 clasp -t 1 --time-limit=3600 > %%OUTPUT%%"' + \
    ' --pipe'

REHCSTAR_MGR = \
    './reHCstar-mgr.py' + \
    ' --verbose' + \
    ' --block-length 10000' + \
    ' --pedigree {genpedfile}' + \
    ' --recomb-upper-limit {recombul}' + \
    ' --results {hapoutfile}' + \
    ' --cmd=\'./reHCstar -3' + \
    '  --global-error --global-error-number {errors}' + \
    '  --pedigree "{{pedigree}}" --haplotypes "{{haplotypes}}" --assumptions "{{assumptions}}"' + \
    '  --sat-cmdline "nice -n19 ./clasp -t 1 --time-limit=1200 > %%OUTPUT%%"' + \
    '  --pipe\''

# SAT solver: internal
#REHCSTAR_CHECK_ERRORS = \
#    './reHCstar -4' + \
#    ' --pedigree {infile} --haplotypes {outfile}' + \
#    ' --global-error --global-error-number {errors}'

#REHCSTAR_MGR = \
#    './reHCstar-mgr.py' + \
#    ' --verbose' + \
#    ' --block-length 10000' + \
#    ' --pedigree {genpedfile}' + \
#    ' --recomb-upper-limit {recombul}' + \
#    ' --results {hapoutfile}' + \
#    ' --cmd=\'./reHCstar -4' + \
#    '  --global-error --global-error-number {errors}' + \
#    '  --pedigree "{{pedigree}}" --haplotypes "{{haplotypes}}" --assumptions "{{assumptions}}"\''


import argparse
import collections
import itertools
import logging
import math
import re
import subprocess
import sys
import os
import glob

class Individual:
    def __init__(self, iid, progrid=-1):
        self.iid = iid
        self.progrid = progrid
        self.father = None
        self.mother = None
        self.gender = GENDER_MALE
        self.gender_given = False

    def __str__(self):
        return "<ID: {0} [ {1} ]>".format(
            self.iid,
            ", ".join([ "{0}: '{1}'".format(key, value)
                        for key,value in self.__dict__.items() ] ))

    def __repr__(self):
        return str(self)

    def is_founder(self):
        assert (self.father is None) == (self.mother is None)
        return not self.father

    def is_HD_genotyped(self):
        return 'HD_genotype' in self.__dict__

    def is_LD_genotyped(self):
        return 'LD_genotype' in self.__dict__ and self.__dict__['LD_genotype'] != MISSING


class Pedigree(dict):
    def __init__(self):
        dict.__init__(self)
        self.LD_mapping = collections.defaultdict(itertools.count(start=1).next)
        self.n_HD_loci = 0

    def __getitem__(self, key):
        if key not in self:
            ind = Individual(key, len(self)+1)
            self[key] = ind
        return dict.__getitem__(self, key)

    def __setitem__(self, key, val):
        if key in self:
            raise Exception('The pedigree is immutable. '
                            'Trying to replace individual \'{0}\.'.format(dict.__getitem__[key]))
        if not isinstance(val, Individual):
            raise TypeError('Trying to insert a non-individual into the pedigree.')
        return dict.__setitem__(self, key, val)

def non_negative_integer(string):
    value = None
    try:
        value = int(string)
    except:
        raise argparse.ArgumentTypeError("'{0}' is not a non-negative integer.".format(string))
    if value < 0:
        raise argparse.ArgumentTypeError("'{0}' is not a non-negative integer.".format(string))
    return value

def input_file_arg(string):
    if string == '-':
        return sys.stdin
    value = None
    try:
        value = file(string, mode='r')
    except:
        raise argparse.ArgumentTypeError("Impossible to open file '{0}' for reading.".format(string))
    return value


def interesting_lines(f, maxsplit=-1, skipfirstlines=0):
    i = 0
    for line in f:
        i = i+1
        if i<=skipfirstlines:
            continue
        if line.startswith("#"):
            continue   # Skip comments
        line = line.strip().split(None, maxsplit) # Split on white-spaces
        if not line:
            continue   # Skip empty lines
        yield line


def read_pedigree(pedf):
    pedigree = Pedigree()
    for line in interesting_lines(pedf):
        assert len(line)==3
        cind = pedigree[line[0]]
        if ( line[1] in INVALID_IDS ) != (line[2] in INVALID_IDS):
            logging.debug("Discarded line '%s'", " ".join(line))
            continue
        if line[1] not in INVALID_IDS:
            cind.father = pedigree[line[1]]
            cind.mother = pedigree[line[2]]
            if cind.father.gender_given and cind.father.gender == GENDER_FEMALE:
                sys.exit("Individual '{0}' who should be the father of '{1}' is female!".format(line[1], line[0]))
            if cind.mother.gender_given and cind.mother.gender == GENDER_MALE:
                sys.exit("Individual '{0}' who should be the mother of '{1}' is male!".format(line[2], line[0]))
            cind.father.gender = GENDER_MALE
            cind.father.gender_given = True
            cind.mother.gender = GENDER_FEMALE
            cind.mother.gender_given = True
    return pedigree


def add_HD_genotypes_to_pedigree(pedigree, genHDf):
    discarded_genotypes = []
    for line in interesting_lines(genHDf, 1):
        assert len(line) == 2
        if line[0] not in pedigree:
            discarded_genotypes.append(line[0])
            continue
        individual = pedigree[line[0]]
        individual.HD_genotype = [ HD_GENOTYPE_MAPPING[g] for g in line[1] ]
        assert not pedigree.n_HD_loci or pedigree.n_HD_loci == len(individual.HD_genotype)
        pedigree.n_HD_loci = len(individual.HD_genotype)

    if discarded_genotypes:
        logging.warn("Discarded the high-density genotype of %d individuals "
                     "which are not present in the pedigree!", len(discarded_genotypes))


def add_LD_genotypes_to_pedigree(pedigree, genLDf):
    discarded_genotypes = []
    # Only the first two columns are considered
    for line in interesting_lines(genLDf, 2):
        if line[1] in INVALID_ALLELES:
            continue
        if line[0] not in pedigree:
            discarded_genotypes.append(line[0])
            continue
        individual = pedigree[line[0]]
        individual.LD_genotype = sorted([ pedigree.LD_mapping[g.group()]
                                          for g in re.finditer('[A-Z][0-9]*', line[1]) ])
        assert len(individual.LD_genotype) == 2

    if discarded_genotypes:
        logging.warn("Discarded the low-density genotype of %d individuals "
                     "which are not present in the pedigree!", len(discarded_genotypes))

def write_pedigree(pedigree, outf, genotype_fn):
    for individual in pedigree.values():
        genotype = genotype_fn(individual)
        outf.write(
            '1\t{id}\t{fid}\t{mid}\t{gender}\tphenotype\t{genotype}\n'.format(
            id=individual.iid,
            fid=individual.father.iid if individual.father else 0,
            mid=individual.mother.iid if individual.mother else 0,
            gender=individual.gender,
            genotype=" ".join([ str(g) for g in genotype]) ) )


parser = argparse.ArgumentParser()
parser.add_argument('--pedigree',
                    type=str, metavar='<file>',
                    default='pedigree.ped')
parser.add_argument('--genotypes-hd',
                    type=str, metavar='<file>',
                    default='genotypes-hd.txt')
parser.add_argument('--map-hd',
                    type=str, metavar='<file>',
                    default='map-hd.txt')
parser.add_argument('--genotypes-ld',
                    type=str, metavar='<file>',
                    default='genotypes-ld.txt')
parser.add_argument('--ld-position',
                    help='the position of the low-density genotype w.r.t. '
                    'the high-density genotype',
                    metavar='<position>',
                    type=non_negative_integer, required=True)
parser.add_argument('--set-ld-always-missing',
                    help='Does not consider the LD genotype while phasing (USE WITH CARE!!)',
                    dest='discard_ld',
                    action='store_true')
parser.add_argument('--max-errors',
                    help='the maximum number of errors allowed in the '
                    'phasing process (default=%(default)s)',
                    metavar='<integer>',
                    dest='max_errors',
                    type=non_negative_integer, default=0)
parser.add_argument('--max-recombinations',
                    help='the maximum number of recombinations allowed in the '
                    'phasing process (default=%(default)s)',
                    metavar='<integer>',
                    dest='max_recombinations',
                    type=non_negative_integer, default=255)
parser.add_argument('output-suffix',
                    nargs='?',
                    type=str,
                    default="output")
parser.add_argument('-v', '--verbose',
                    help='log message verbosity',
                    action='store_true')

args = parser.parse_args()
args = vars(args)

for f in ['pedigree', 'genotypes_hd', 'map_hd', 'genotypes_ld']:
    args[f] = input_file_arg(args[f])

if not args['verbose']:
    log_level = logging.INFO
else:
    log_level = logging.DEBUG

if 'output' not in args:
    args['output'] = 'results-'


logging.basicConfig(level=log_level,
                    format='%(levelname)-8s [%(asctime)s]  %(message)s',
                    datefmt="%y%m%d %H%M%S")




logging.info('Allele association discovery via genotype phasing')

logging.info('Reading pedigree from file \'%s\'...', args['pedigree'].name)
pedigree = read_pedigree(args['pedigree'])
logging.info('Read %d individuals (%d founders).',
             len(pedigree), sum(( 1
                                  for ind in pedigree.values()
                                  if ind.is_founder() ) ))

logging.info('Reading high-density genotypes from file \'%s\'...',
             args['genotypes_hd'].name)
add_HD_genotypes_to_pedigree(pedigree, args['genotypes_hd'])
logging.info('Read %d high-density genotypes on %d loci.',
             sum(( 1 for ind in pedigree.values() if ind.is_HD_genotyped() ) ),
             pedigree.n_HD_loci)
assert args['ld_position'] <= pedigree.n_HD_loci

logging.info('Reading high-density genotype map from file \'%s\'...',
             args['map_hd'].name)
maphd = []
maphdpos = {}
for line in interesting_lines(args['map_hd']):
    maphd.append([line[2], line[1], line[3]])
    maphdpos[line[1]] = len(maphdpos)
assert pedigree.n_HD_loci == len(maphd)
assert pedigree.n_HD_loci == len(maphdpos)

logging.info('Reading low-density genotypes from file \'%s\'...',
             args['genotypes_ld'].name)
add_LD_genotypes_to_pedigree(pedigree, args['genotypes_ld'])
logging.info('Read %d low-density genotypes (%d of them have also high-density genotypes).',
             sum(( 1 for ind in pedigree.values() if ind.is_LD_genotyped() ) ),
             sum(( 1 for ind in pedigree.values() if ind.is_HD_genotyped() and ind.is_LD_genotyped() ) ))
logging.info('Low-density genotypes have %d alleles [%s].',
             len(pedigree.LD_mapping.keys()),
             ", ".join(pedigree.LD_mapping.keys()) )

#logging.info("Preparing files for Mendelian error check with plink (only bi-allelic)...")
#logging.info("...low-density genotypes to file 'tmp-plink-ld.{ped,map}'...")
#
#with file('tmp-plink-ld.ped', mode="w") as outf:
#    write_pedigree(pedigree, outf, lambda individual: ( individual.LD_genotype
#                                                        if ( individual.is_LD_genotyped() and
#                                                             all([int(ldg)<=2 for ldg in individual.LD_genotype]) )
#                                                        else MISSING.split(' ') ) )
#with file('tmp-plink-ld.map', mode="w") as outf:
#    outf.write('1\tLDgen\t1\t1\n')
#
#
#
#logging.info("Executing Mendelian errors checks via 'plink'...")
#subprocess.call("p-link --file tmp-plink-ld --mendel --out tmp-plink-ld-out", shell=True)
#
#
#logging.info("Setting as MISSING genotypes reported with Mendelian errors...")
wrong_ld_genotypes = []
n_errors = 0
#with file('tmp-plink-ld-out.mendel', mode='r') as inf:
#    for line in interesting_lines(inf, skipfirstlines=1):
#        assert line[1] in pedigree
#        logging.warn("Setting low-density genotype of individual %s as missing "
#                     "since it induces a Mendelian error.", line[1])
#
#        pedigree[line[1]].LD_genotype = MISSING
#        wrong_ld_genotypes.append('{id}\t{original}\t{new}\n'.format(id=line[1],
#                                                                     original=line[9].replace("/",""),
#                                                                     new=MISSING))
#        n_errors += 1
#if n_errors>0:
#    logging.warn("Found %d Mendelian errors.", n_errors)



logging.info("Preparing files for Mendelian error check with plink (only bi-allelic)...")
logging.info("...high-density genotypes to file 'tmp-plink-hd.{ped,map}'...")

with file('tmp-plink-hd.ped', mode="w") as outf:
    write_pedigree(pedigree, outf, lambda individual: ( individual.HD_genotype
                                                        if individual.is_HD_genotyped()
                                                        else [ MISSING ] * pedigree.n_HD_loci ) )
with file('tmp-plink-hd.map', mode='w') as outf:
    for snp in maphd:
        outf.write('\t'.join(snp)+'\n')

logging.info("Executing Mendelian errors checks via 'plink'...")
subprocess.call("p-link --file tmp-plink-hd --mendel --out tmp-plink-hd-out --map3", shell=True)

logging.info("Setting as MISSING genotypes reported with Mendelian errors...")
n_errors = 0
with file('tmp-plink-hd-out.mendel', mode='r') as inf:
    for line in interesting_lines(inf, maxsplit=4, skipfirstlines=1):
        assert line[1] in pedigree
        assert line[3] in maphdpos
        logging.warn("Setting high-density genotype of individual %s at position %d "
                     "as missing since it induces a Mendelian error.",
                     line[1], maphdpos[line[3]])
        pedigree[line[1]].HD_genotype[maphdpos[line[3]]] = MISSING
        n_errors += 1
if n_errors>0:
    logging.warn("Found %d Mendelian errors.", n_errors)

logging.info("Preparing files for Mendelian error check with reHCstar...")
logging.info("...low-density genotypes to file 'tmp-rehcstar-ld.ped'...")

with file('tmp-rehcstar-ld.ped', mode="w") as outf:
    write_pedigree(pedigree, outf, lambda individual: ( individual.LD_genotype
                                                        if individual.is_LD_genotyped()
                                                        else MISSING.split(' ') ) )


# Invert LD_mapping
inv_LD_mapping = dict([ (pedigree.LD_mapping[x], x) for x in pedigree.LD_mapping ])
def to_LD_genotype(g):
    return "".join([ inv_LD_mapping[x] for x in g ])

min_errors = -1
max_errors = 0
optimum_found = False
best_errors = None
while not optimum_found:
    logging.info("Checking if a solution with at most %d errors exists..", max_errors)
    retcode= subprocess.call(REHCSTAR_CHECK_ERRORS.format(infile='tmp-rehcstar-ld.ped',
                                                          outfile='tmp-rehcstar-ld.out',
                                                          errors=max_errors),
                             shell=True)
    if retcode == 0: # A solution exists
        logging.info("A solution with %d errors has been found!", max_errors)
        best_errors = []
        with file('tmp-rehcstar-ld.out', mode='r') as inf:
            for line in interesting_lines(inf):
                assert line[1] in pedigree
                if not pedigree[line[1]].is_LD_genotyped():
                    continue
                solution_genotype = sorted([int(g) for g in line[6].split('|')])
                if solution_genotype != pedigree[line[1]].LD_genotype:
                    best_errors.append( (line[1],
                                         to_LD_genotype(pedigree[line[1]].LD_genotype),
                                         to_LD_genotype(solution_genotype) ) )
        max_errors = (len(best_errors)+min_errors) // 2
    else:
        logging.info("A solution with %d errors has NOT been found!", max_errors)
        if best_errors is None:
            min_errors = max_errors
            max_errors = max(1, 2 * max_errors)
        else:
            min_errors = max_errors
            max_errors = (len(best_errors)+min_errors) // 2

    if best_errors is not None and (min_errors + 1 >= len(best_errors)):
        logging.info("An optimum solution with %d errors has been found!", len(best_errors))
        optimum_found = True

logging.info("Saving LD errors (if present) to file '%s'...", outf.name)
if len(best_errors) > 0:
    logging.info("Reading the solution and setting as missing the wrong genotypes...")
    for error in best_errors:
        logging.warn("Low-density genotype of individual %s is wrong.",
                     error[0])
        wrong_ld_genotypes.append(
            '{id}\t{original}\t{new}\n'.format(id=error[0],
                                               original=error[1],
                                               new=error[2]))
        pedigree[error[0]].LD_genotype = MISSING.split(" ")

with file('{0}wrong-ld-genotypes.txt'.format(args['output']),
          mode='w') as outf:
    for x in wrong_ld_genotypes:
        outf.write(x)



logging.info('Removing temporary files...')
for fl in glob.glob("./tmp-*"):
    logging.info('...removing {0}'.format(fl))
    os.remove(fl)



rehcstarf_name = '{0}gen-ped-hd_plus_ld-ld_pos{1}.txt'.format(args['output'],
                                                              args['ld_position'])
with file(rehcstarf_name, mode='w') as rehcstarf:
    logging.info("Preparing the reHCstar-mgr instance in file '%s'...",
                 rehcstarf.name)
    for individual in pedigree.values():
        HD_genotype = ( individual.HD_genotype if individual.is_HD_genotyped()
                        else [ MISSING ] * pedigree.n_HD_loci )
        LD_genotype = ( individual.LD_genotype if ( individual.is_LD_genotyped() and
                                                    not args['discard_ld'] )
                        else MISSING.split(' ') )
        rehcstarf.write(
            '1\t{id}\t{fid}\t{mid}\t{gender}\tphenotype\t{genotype}\n'.format(
                id=individual.iid,
                fid=individual.father.iid if individual.father else 0,
                mid=individual.mother.iid if individual.mother else 0,
                gender=individual.gender,
                genotype="\t".join( (HD_genotype[:args['ld_position']] +
                                     [" ".join([str(ldg) for ldg in LD_genotype])] +
                                     HD_genotype[args['ld_position']:] ) )
                )
            )
    logging.info('Instance prepared!')


logging.info("Executing reHCstar-mgr on file '%s'...", rehcstarf_name)
hapoutf_name = '{0}hap-gen-ped-hd_plus_ld-ld_pos{1}.txt'.format(args['output'],
                                                                args['ld_position'])

n_haplotyping_errors = 0
recomb_ul = args['max_recombinations']
solution_found = False
retcode = 0
while ( not solution_found and
        n_haplotyping_errors <= args['max_errors'] ):
    cmd_line = REHCSTAR_MGR.format(genpedfile=rehcstarf_name,
                                   hapoutfile=hapoutf_name,
                                   recombul=recomb_ul,
                                   errors=n_haplotyping_errors)
    logging.info("Command-line: >>>%s<<<", cmd_line)
    retcode = subprocess.call(cmd_line, shell=True)
    if retcode != 0:
        logging.warning("A solution with %d errors has not been found.%s",
                        n_haplotyping_errors,
                        " Retrying allowing more errors..." if ( n_haplotyping_errors < args['max_errors']) else "")
        n_haplotyping_errors += 1
        recomb_ul = int(math.ceil((recomb_ul + 1.0) / 2.0))-1
    else:
        solution_found = True

if not solution_found:
    logging.error("reHCstar-mgr failed on input '%s'! Please check and adjust the parameters...",
                  rehcstarf_name)
    sys.exit(128+retcode)

logging.info('reHCstar-mgr terminated!')

logging.info('Removing temporary files...')
for fl in glob.glob("./tmp-*"):
    logging.info('...removing {0}'.format(fl))
    os.remove(fl)

logging.info('Cleaning the output haplotypes...')
hapouthdf_name = '{0}hap-gen-ped-hd-ld_pos{1}.txt'.format(args['output'],
                                                          args['ld_position'])
with file(hapoutf_name, mode='r') as inf:
    with file(hapouthdf_name, mode='w') as outf:
        LDpos = args['ld_position']
        for line in inf:
            if re.search("^## GENOTYPE   IND\[", line):
                outf.write(line[0:29+(LDpos*3)])
                outf.write(line[32+(LDpos*3):])
            elif ( re.search("^## ([PM]AT_(HAPLOT|SOUR_V)|MASK_ERR  ) IND\[", line) and
                   len(line)>31+LDpos ):
                outf.write(line[0:30+LDpos])
                outf.write(line[31+LDpos:])
            elif line.startswith("1\t"):
                vline = line.strip().split("\t", 6)
                if vline[1] not in pedigree:
                    logging.debug("Read an extraneous individual '%s'. Copying...", vline[1])
                    outf.write(line)
                else:
                    outf.write("\t".join(vline[0:6]))
                    outf.write("\t")
                    vline = vline[6].split("\t")
                    outf.write("\t".join(vline[0:LDpos] + vline[LDpos+1:]))
                    outf.write("\n")
            else:
                outf.write(line)


logging.info('Completed!')
