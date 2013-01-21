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

import collections
import gzip
import itertools
import json
import logging
import optparse
import re
import sys

import hapallelepred_utils

import cplex

genotype_mapping = { '00': '0',
                     '11': '1',
                     '12': '3',
                     '21': '3',
                     '22': '2' }

class VarType:
    (E, A, P)= (0, 1, 2)

class Individual:
    def __init__(self, iid):
        self.iid= iid
        self.high_genotyped= False
        self.low_genotyped= False

    def __str__(self):
        return "<ID: {0}\n{1}>\n".format(self.iid,
                                         "\n".join([ "\t\t{0}: '{1}'".format(key,self.__dict__[key])
                                                     for key in self.__dict__ ] ))

    def __repr__(self):
        return str(self)

class LogFile:
    def __init__(self, prefix=""):
        self._buff= []
        self._prefix= prefix

    def _do_print(self):
        logging.debug("%s%s",
                      self._prefix,
                      "".join(self._buff))
        self._buff= []

    def __enter__(self):
        if self._buff:
            self._do_print()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self._buff:
            self._do_print()
        if exc_type is not None:
            logging.warn("%sERROR: %s", self._prefix, exc_value)
        return False

    def write(self, string):
        outstr= "{0}".format(string)
        if "\n" not in outstr:
            self._buff.append(outstr)
        else:
            lines= outstr.split("\n", 1)
            assert len(lines)==2
            self._buff.append(lines[0])
            self._do_print()
            self.write(lines[1])
        return self

    def flush(self):
        return self

def hamming_distance(s1, s2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def print_haplotype_distances(haplotype_distances,
                              haplotype_length,
                              filename="similarity-matrix.tsv"):

    sortperm= [ (i, sum(haplotype_distances[i])) for i in range(len(haplotype_distances)) ]
    sortperm.sort(key=lambda x: x[1])
    new_haplotype_distances= [ [ haplotype_distances[sortperm[i][0]][sortperm[j][0]]/float(haplotype_length)
                                 for j in range(len(haplotype_distances)) ]
                               for i in range(len(haplotype_distances)) ]
    with open(filename, "w") as sdf:
        print("\n".join([ "\t".join(["{0}".format(e) for e in row])
                          for rown, row in zip(itertools.count(),
                                               new_haplotype_distances) ]), file=sdf)

    return None


def prepare_ILP_variables(pedigree, haplotypes, low_density_alleles):
    variable_list= []
    variable_names= []
    variable_obj= []
# Variables E: 2 for each individual
    for iid in pedigree:
        variable_list.extend([ (VarType.E, iid, 0), (VarType.E, iid, 1) ])
        variable_names.extend([ "E_{0}_0".format(iid), "E_{0}_1".format(iid) ])
        variable_obj.extend( [ 1.0, 1.0 ] )

# Variables A: 1 for each low-density allele and each haplotype
    for hid in range(len(haplotypes)):
        for ldaid in range(len(low_density_alleles)):
            variable_list.append((VarType.A, hid, ldaid))
            variable_names.append("A_{0}_{1}".format(hid, ldaid))
            variable_obj.append( 0.0 )

# Variables P: 2 for heterozygous genotypes and zero for homozygous genotypes
    for ind in pedigree.values():
        if ind.low_genotype[0] != ind.low_genotype[1]: # heterozygous
            variable_list.extend([ (VarType.P, ind.iid, 0), (VarType.P, ind.iid, 1) ])
            variable_names.extend([ "P_{0}_0".format(ind.iid), "P_{0}_1".format(ind.iid) ])
            variable_obj.extend( [ 0.0, 0.0 ] )

    variable_idx= { v:i for i,v in enumerate(variable_list) }

    return (variable_names, variable_obj, variable_idx)




def prepare_ILP_constraints(pedigree, haplotypes, low_density_alleles,
                            variable_idx):

    def get_constr_coeff(variable_idx, *variables):
        coeff= cplex.SparsePair(ind= [ variable_idx[var] for c,var in variables ],
                                val= [ c for c,var in variables ])
        return coeff

    constr_mat= []
    constr_rhs= []
    constr_sense= []
    for iid,ind in pedigree.items():
        # Phase constraints
        if ind.low_genotype[0] != ind.low_genotype[1]: # heterozygous
            constr_mat.append( get_constr_coeff(variable_idx,
                                                (1, (VarType.P, iid, 0)),
                                                (1, (VarType.P, iid, 1)) ) )
            constr_sense.append('E')
            constr_rhs.append(1)

    # Error constraints
        if ind.low_genotype[0] != ind.low_genotype[1]: # heterozygous
            constr_mat.append( get_constr_coeff(variable_idx,
                                                (1, (VarType.E, iid, 0)),
                                                (1, (VarType.A,
                                                     haplotypes[ind.pat_haplotype],
                                                     low_density_alleles[ind.low_genotype[0]])),
                                                (1, (VarType.P, iid, 0)) ) )
            constr_sense.append('G')
            constr_rhs.append(1)
            constr_mat.append( get_constr_coeff(variable_idx,
                                                (1, (VarType.E, iid, 0)),
                                                (1, (VarType.A,
                                                     haplotypes[ind.pat_haplotype],
                                                     low_density_alleles[ind.low_genotype[1]])),
                                                (1, (VarType.P, iid, 1)) ) )
            constr_sense.append('G')
            constr_rhs.append(1)

            constr_mat.append( get_constr_coeff(variable_idx,
                                                (1, (VarType.E, iid, 1)),
                                                (1, (VarType.A,
                                                     haplotypes[ind.mat_haplotype],
                                                     low_density_alleles[ind.low_genotype[0]])),
                                                (1, (VarType.P, iid, 1)) ) )
            constr_sense.append('G')
            constr_rhs.append(1)
            constr_mat.append( get_constr_coeff(variable_idx,
                                                (1, (VarType.E, iid, 1)),
                                                (1, (VarType.A,
                                                     haplotypes[ind.mat_haplotype],
                                                     low_density_alleles[ind.low_genotype[1]])),
                                                (1, (VarType.P, iid, 0)) ) )
            constr_sense.append('G')
            constr_rhs.append(1)

        else:  # homozygous
            constr_mat.append( get_constr_coeff(variable_idx,
                                                (1, (VarType.E, iid, 0)),
                                                (1, (VarType.A,
                                                     haplotypes[ind.pat_haplotype],
                                                     low_density_alleles[ind.low_genotype[0]])) ) )
            constr_sense.append('G')
            constr_rhs.append(1)
            constr_mat.append( get_constr_coeff(variable_idx,
                                                (1, (VarType.E, iid, 1)),
                                                (1, (VarType.A,
                                                     haplotypes[ind.mat_haplotype],
                                                     low_density_alleles[ind.low_genotype[0]])) ) )
            constr_sense.append('G')
            constr_rhs.append(1)

    for hid in range(len(haplotypes)):
        a_constr= []
        for ldaid in range(len(low_density_alleles)):
            a_constr.append((1, (VarType.A, hid, ldaid)))
        constr_mat.append(get_constr_coeff(variable_idx, *a_constr))
        constr_sense.append('E')
        constr_rhs.append(1)

    return (constr_mat, constr_sense, constr_rhs)



parser= optparse.OptionParser(usage="usage: "
                              "%prog -i <HAPLOTYPE FILE> -g <LOW-DENSITY GENOTYPE FILE> [-v]")
parser.add_option("-i", "--haplotypes",
                  action="store", dest="hapfile",
                  type="string", default=None,
                  help="The file containing the input phased genotypes (as produced by reHCstar).",
                  metavar="FILE")
parser.add_option("-g", "--low-density-genotypes",
                  action="store", dest="lowgenfile",
                  type="string", default=None,
                  help="The file containing the input low-density genotypes.",
                  metavar="FILE")
parser.add_option("--keep-not-genotyped",
                  action="store_true", dest="keep_ungen",
                  default=False,
                  help="Do not discard not-genotyped individuals.")
parser.add_option("--trim-borders",
                  action="store", dest="border",
                  type="int", default=-1,
                  help="The number of loci to discard at the haplotype borders.")
parser.add_option("--trim-left-border",
                  action="store", dest="lborder",
                  type="int", default=-1,
                  help="The number of loci to discard at the haplotype left border.")
parser.add_option("--trim-right-border",
                  action="store", dest="rborder",
                  type="int", default=-1,
                  help="The number of loci to discard at the haplotype right border.")
parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose",
                  default=False,
                  help="Print additional log messages.")

(options, args)= parser.parse_args()

log_level= logging.DEBUG if options.verbose else logging.INFO

logging.basicConfig(level=log_level,
                    format='%(levelname)-8s [%(asctime)s]  %(message)s',
                    datefmt="%y%m%d %H%M%S")


logging.info("Haplotype-based Low-density Genotype Phasing")

hapallelepred_utils.check_file(options.hapfile)
hapallelepred_utils.check_file(options.lowgenfile)

if options.border<0:
    options.border= 0
if options.lborder<0:
    options.lborder= options.border
if options.rborder<0:
    options.rborder= options.border

logging.info("Reading haplotypes from file '%s'...", options.hapfile)
pedigree= {}
haplotype_length= -1
with open(options.hapfile, "r") as hf:
    for line in hf:
        if ( not line.startswith("## ")
             or " IND[" not in line ):
            continue
        line= line.strip()

        iid = line[18:28].strip()
        if iid not in pedigree:
            pedigree[iid] = Individual(iid)
        ind = pedigree[iid]
        data = line[29:].strip()

        if line.startswith("## INDIVIDUAL IND["):
            ind.header= data
            a= data.split()
            ind.father= a[2]
            ind.mother= a[3]

        if line.startswith("## GENOTYPE   IND["):
            data= "".join([genotype_mapping[g] for g in data.split(" ")])
            assert haplotype_length==-1 or len(data)==haplotype_length
            haplotype_length= len(data)
            ind.high_genotype= data
            ind.high_genotyped= "1" in data or "2" in data or "3" in data

        if line.startswith("## PAT_HAPLOT IND["):
            assert haplotype_length==-1 or len(data)==haplotype_length
            haplotype_length= len(data)
            ind.pat_haplotype= data
        if line.startswith("## MAT_HAPLOT IND["):
            assert haplotype_length==-1 or len(data)==haplotype_length
            haplotype_length= len(data)
            ind.mat_haplotype= data

        if line.startswith("## PAT_SOUR_V IND["):
            ind.pat_source_vector= data
        if line.startswith("## MAT_SOUR_V IND["):
            ind.mat_source_vector= data

logging.info("Reading low-density genotypes from file '%s'...", options.lowgenfile)
with open(options.lowgenfile, "r") as lgf:
    for line in lgf:
        line= line.strip()
        line= line.split(None, 1)

        if line[0] not in pedigree:
#            logging.debug("Low-density genotyped individual '%s' has no haplotypes. Discarding...",
#                          line[0])
            continue

        pedigree[line[0]].low_genotype= line[1].strip()
        pedigree[line[0]].low_genotyped= True

logging.info("No. of individuals: %6d", len(pedigree))

# Trim haplotype borders
if options.lborder + options.rborder > 0:
    logging.debug("Trimming haplotype borders...")
    lborder=options.lborder
    rborder=options.rborder
    assert lborder + rborder < haplotype_length, "Borders must be less than haplotypes' length."
    if rborder>0:
        for ind in pedigree.values():
            ind.high_genotype= ind.high_genotype[lborder:-rborder]
            ind.pat_haplotype= ind.pat_haplotype[lborder:-rborder]
            ind.mat_haplotype= ind.mat_haplotype[lborder:-rborder]
            ind.pat_source_vector= ind.pat_source_vector[lborder:-rborder]
            ind.mat_source_vector= ind.mat_source_vector[lborder:-rborder]
    else:
        for ind in pedigree.values():
            ind.high_genotype= ind.high_genotype[lborder:]
            ind.pat_haplotype= ind.pat_haplotype[lborder:]
            ind.mat_haplotype= ind.mat_haplotype[lborder:]
            ind.pat_source_vector= ind.pat_source_vector[lborder:]
            ind.mat_source_vector= ind.mat_source_vector[lborder:]
    for ind in pedigree.values():
        ind.high_genotyped= "1" in ind.high_genotype or "2" in ind.high_genotype or "3" in ind.high_genotype
    haplotype_length= haplotype_length - lborder - rborder

# Remove un-genotyped individuals
if not options.keep_ungen:
    logging.debug("Removing not-genotyped individuals...")
    ungen= [ iid for iid in pedigree if not pedigree[iid].high_genotyped]
    for iid in ungen:
        del pedigree[iid]
    logging.info("Removed %d not-genotyped individuals.", len(ungen))
    logging.info("No. of remaining individuals: %6d", len(pedigree))


logging.debug("Removing individuals without low-density genotypes...")
ungen= [ iid for iid in pedigree if not pedigree[iid].low_genotyped ]
for iid in ungen:
    del pedigree[iid]
logging.info("Removed %d individuals without low-density genotypes.", len(ungen))
logging.info("No. of remaining individuals: %6d", len(pedigree))
del ungen

logging.debug("Parsing low-density genotypes...")
lowgen_pattern= re.compile("^(?P<allele1>[A-Za-z]+[0-9]*)\s*(?P<allele2>[A-Za-z]+[0-9]*)$")
remove_iid= []
low_density_alleles= {}
for ind in pedigree.values():
    assert ind.low_genotyped
    match_info= lowgen_pattern.match(ind.low_genotype)
    if match_info is None:
        logging.warn("Low-density genotype '%s' of individual '%s' NOT RECOGNIZED. Discarding individual...",
                     ind.low_genotype, ind.iid)
        remove_iid.append(ind.iid)
    else:
        ind.low_genotype= match_info.group('allele1', 'allele2')
        low_density_alleles[ind.low_genotype[0]]= 0
        low_density_alleles[ind.low_genotype[1]]= 0

del lowgen_pattern

if remove_iid:  ## Some low-density genotypes have not been parsed correctly
    for iid in remove_iid:
        logging.debug("Removing individual '%s'...", iid)
        del pedigree[iid]
logging.info("Removed %d individuals with invalid low-density genotypes.", len(remove_iid))
logging.info("No. of remaining individuals: %6d", len(pedigree))
del remove_iid

logging.info("The low-density genotypes have %d alleles: '%s'.",
             len(low_density_alleles), "', '".join(low_density_alleles.keys()))

for i,allele in enumerate(low_density_alleles.keys()):
    low_density_alleles[allele]= i


haplotype_list= ( [ ind.pat_haplotype for ind in pedigree.values() ] +
                  [ ind.mat_haplotype for ind in pedigree.values() ] )
haplotypes= dict.fromkeys(haplotype_list)

for i,hap in enumerate(haplotypes.keys()):
    haplotypes[hap]= i

haplotype_frequency= collections.Counter(haplotype_list)
del haplotype_list

logging.info("No. of distinct haplotypes:   %6d", len(haplotypes))
logging.debug("The most frequent haplotypes are:")
for el in haplotype_frequency.most_common(10):
    logging.debug("(%6d indiv.) %s", el[1], el[0])
private_haplotypes= [x for x, cnt in haplotype_frequency.items() if cnt==1]
logging.info("No. of 'private' haplotypes:  %6d", len(private_haplotypes))
logging.debug("The haplotypes that are not shared by two individuals are:")
for i, el in enumerate(sorted(private_haplotypes)):
    logging.debug("(%13d) %s", i+1, el)
del private_haplotypes


logging.debug("Computing the (dis)similarity matrix among the haplotypes...")
haplotype_distances= [ [0]*len(haplotypes) for i in range(len(haplotypes)) ]
for i1,h1 in zip(itertools.count(), haplotypes):
    for i2,h2 in zip(itertools.count(), haplotypes):
        if i1==i2:
            break
        haplotype_distances[i1][i2]= haplotype_distances[i2][i1]= hamming_distance(h1, h2)

logging.debug("Saving the (dis)similarity matrix among the haplotypes...")
print_haplotype_distances(haplotype_distances, haplotype_length,
                          "similarity-matrix.tsv")

logging.info("Discovering associations...")


logging.debug("Preparing the integer-linear program...")
## Prepare the integer program

logging.debug("Computing the set of variables...")

(variable_names, variable_obj, variable_idx) = prepare_ILP_variables(pedigree, haplotypes, low_density_alleles)
N_VARS= len(variable_names)

logging.debug("The program has %d variables.", N_VARS)


logging.debug("Computing the set of constraints...")
(constr_mat, constr_sense, constr_rhs) = prepare_ILP_constraints(pedigree, haplotypes, low_density_alleles,
                                                                 variable_idx)

logging.debug("The program has %d constraints.", len(constr_mat))

c = cplex.Cplex()
with LogFile("CPLEX_LOG == ") as lf:

    c.set_results_stream(lf)
    c.set_log_stream(lf)
    c.set_warning_stream(lf, lambda x: "!! WARNING !! == " + x)
    c.set_error_stream(lf,   lambda x: "!!  ERROR  !! == " + x)

    c.objective.set_sense(c.objective.sense.minimize)

    c.variables.add(names = variable_names,
                    obj = variable_obj,
                    lb = [ 0.0 ] * N_VARS,
                    ub = [ 1.0 ] * N_VARS,
                    types = [ c.variables.type.binary ] * N_VARS)



    c.linear_constraints.add(names = [ "c{0}".format(i) for i in range(len(constr_mat)) ],
                             lin_expr = constr_mat,
                             senses = constr_sense,
                             rhs = constr_rhs)

    if options.verbose:
        logging.debug("Saving problem to file 'ilp-problem.mps'...")
        c.write('ilp-problem.mps', 'mps')

    logging.info("Solving the problem...")
    c.solve()

logging.info("Search process terminated!")

STATUS = cplex.Cplex.solution.status
OPTIMAL = [ STATUS.optimal, STATUS.MIP_optimal, STATUS.MIP_optimal_relaxed_sum ]
INFEASIBLE = [ STATUS.infeasible, STATUS.MIP_infeasible, STATUS.MIP_infeasible_or_unbounded ]
logging.debug("The solver status is: %s.", STATUS[c.solution.get_status()])
logging.debug("Solving method: %s.", cplex.Cplex.solution.method[c.solution.get_method()])

def str_src_vect(src):
    if not src:
        return ""
    prev_c = src[0]
    result = []
    for c in src:
        result.append(" " if prev_c == c else "/")
        prev_c = c
    return "".join(result)

if c.solution.get_status() in OPTIMAL:
    logging.info("Optimum found! -- There are %d wrong associations over %d.",
                 c.solution.get_objective_value(), len(pedigree)*2)
    haplotype_assignment = {}
    for h in haplotypes:
        for lda in low_density_alleles:
            if c.solution.get_values(variable_idx[(VarType.A, haplotypes[h],
                                                   low_density_alleles[lda])]) > 0 :
                assert h not in haplotype_assignment
                haplotype_assignment[h]= lda

    logging.info("Analyzing haplotype assignments...")
    with open("ilp-problem-solution.txt", mode="w") as outf:
        for var in variable_names:
            print("{:12}= {:5}".format(var, "True" if c.solution.get_values(var)>0.0 else "False"), file=outf)

    assert set(haplotype_assignment.keys()) == set(haplotype_frequency.keys())
    with gzip.open("haplotype-assignment.json.gz", mode="wb") as outf:
        haplotype_data = collections.defaultdict(dict)
        for haplotype, allele in haplotype_assignment.items():
            haplotype_data[haplotype]["assigned_allele"]= allele
        for haplotype, frequency in haplotype_frequency.items():
            haplotype_data[haplotype]["frequency"]= frequency
        for haplotype, frequency in haplotype_frequency.items():
            tmp_ind = []
            for ind in pedigree.values():
                if ind.pat_haplotype==haplotype:
                    tmp_ind.append({"individual": ind.iid,
                                    "source": "paternal",
                                    "LD genotype": "".join(ind.low_genotype),
                                    "parent": ind.father,
                                    "parent LD genotype": "".join(pedigree[ind.father].low_genotype) if ind.father in pedigree else "--",
                                    "wrong association": c.solution.get_values(variable_idx[(VarType.E, ind.iid, 0)]) > 0} )
                    tmp_src= str_src_vect(ind.pat_source_vector)
                    if "/" in tmp_src:
                        tmp_ind[-1]["source vector"] = tmp_src
            for ind in pedigree.values():
                if ind.mat_haplotype==haplotype:
                    tmp_ind.append({"individual": ind.iid,
                                    "source": "maternal",
                                    "LD genotype": "".join(ind.low_genotype),
                                    "parent": ind.mother,
                                    "parent LD genotype": "".join(pedigree[ind.mother].low_genotype) if ind.mother in pedigree else "--",
                                    "wrong association": c.solution.get_values(variable_idx[(VarType.E, ind.iid, 1)]) > 0} )
                    tmp_src= str_src_vect(ind.mat_source_vector)
                    if "/" in tmp_src:
                        tmp_ind[-1]["source vector"] = tmp_src
            haplotype_data[haplotype]["supporting info"] = tmp_ind
        json_out = { "border_left": options.lborder,
                     "border_right": options.rborder,
                     "haplotypes": haplotype_data,
                     }
        json.dump(json_out, fp=outf, sort_keys=True, indent=2)

    with open("haplotype-assignment-analysis.txt", mode="w") as outf:
        MAX_DISTANCE= max(10, int(0.05 * haplotype_length))
        print("## Considering neighbors at distance at most {0}.".format(MAX_DISTANCE),
              file=outf)
        print("=" * (haplotype_length + 28), file=outf)
#        for h1,a in haplotype_assignment.items():
        for h1, h1count in sorted(haplotype_frequency.items(), key=lambda x: -x[1]):
            a= haplotype_assignment[h1]
            h1id= haplotypes[h1]
            print("Haplotype id:{0:<4}    ({1:<2})   {2}".format(h1id, a, h1), file=outf)
            print("-" * (haplotype_length + 28), file=outf)
            print("Neighbors:", file=outf)
            for h2,h2id in sorted(haplotypes.items(), key=lambda x: haplotype_distances[h1id][x[1]]):
                if h2id == h1id:
                    continue
                if haplotype_distances[h1id][h2id] <= MAX_DISTANCE:
                    print("    id:{0:<4} dist:{1:<3} ({2:<2}) {3} {4}".format(
                            h2id,
                            haplotype_distances[h1id][h2id],
                            haplotype_assignment[h2],
                            " " if haplotype_assignment[h1]==haplotype_assignment[h2] else "#",
                            "".join( [ " " if a1==a2 else a2 for a1,a2 in zip([ch for ch in h1], [ch for ch in h2]) ])),
                          file=outf)
            print("-" * (haplotype_length + 28), file=outf)
            print("Individuals:", file=outf)
            tmp_ind = collections.defaultdict(list)
            for ind in pedigree.values():
                if ind.pat_haplotype==h1:
                    tmp_ind[ind.father].append("  Individual {:10} (paternal)"
                                               "  Low-density genotype: {:4}"
                                               "  Parent (father): {:10} ({:4})".format(ind.iid,
                                                                                        "".join(ind.low_genotype),
                                                                                        ind.father,
                                                                                        "".join(pedigree[ind.father].low_genotype) if ind.father in pedigree else "--"))
                    if c.solution.get_values(variable_idx[(VarType.E, ind.iid, 0)]) > 0:
                        tmp_ind[ind.father][-1] = "@@" + tmp_ind[ind.father][-1][2:]
                    tmp_src= str_src_vect(ind.pat_source_vector)
                    if "/" in tmp_src:
                        tmp_ind[ind.father][-1] = tmp_ind[ind.father][-1] + "\n    Paternal source vector: {0}".format(tmp_src)
            for par in sorted(tmp_ind.keys(), key=int):
                print("\n".join(tmp_ind[par]), file=outf)
            tmp_ind = collections.defaultdict(list)
            for ind in pedigree.values():
                if ind.mat_haplotype==h1:
                    tmp_ind[ind.mother].append("  Individual {:10} (maternal)"
                                               "  Low-density genotype: {:4}"
                                               "  Parent (mother): {:10} ({:4})".format(ind.iid,
                                                                                        "".join(ind.low_genotype),
                                                                                        ind.mother,
                                                                                        "".join(pedigree[ind.mother].low_genotype) if ind.mother in pedigree else "--"))
                    if c.solution.get_values(variable_idx[(VarType.E, ind.iid, 1)]) > 0:
                        tmp_ind[ind.mother][-1] = "@@" + tmp_ind[ind.mother][-1][2:]
                    tmp_src= str_src_vect(ind.mat_source_vector)
                    if "/" in tmp_src:
                        tmp_ind[ind.mother][-1] = tmp_ind[ind.mother][-1] + "\n    Maternal source vector: {0}".format(tmp_src)
            for par in sorted(tmp_ind.keys(), key=int):
                print("\n".join(tmp_ind[par]), file=outf)
            print("", file=outf)
            print("=" * (haplotype_length + 28), file=outf)
            print("", file=outf)





elif c.solution.get_status() in INFEASIBLE:
    logging.warn("No feasible solution exists!! Aborting...")
    sys.exit("No feasible solution exists.")
else:
    logging.error("Solver status '%s' unknown! Aborting...", STATUS[c.solution.get_status()])
    sys.exit("Unknown solver status.")





# get the solution
#lpstat = CPX.getstat(env,lp)
#if lpstat in OPTIMAL:
#    # can call getx, getobjval, too
#    x, objval = CPX.getmipx(env,lp), CPX.getmipobjval(env,lp)
#    print(x, objval)
#    for h in haplotypes:
#        ass= []
#        for lda in low_density_alleles:
#            if x[variables[(VarType.A, haplotypes[h], low_density_alleles[lda])]] > 0:
#                ass.append(lda)
#        assert len(ass)==1
#        print("{0} ---> {1}".format(h, ass[0]))


logging.info("Haplotype-based Low-density Genotype Phasing -- Completed")
