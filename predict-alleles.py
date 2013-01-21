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
import json
import logging
import math
import optparse
import sys

import hapallelepred_utils

## Try to import some time-consuming routines written with 'pyx'.
## If it is not available, use the regular versions
using_fast = False
try:
    import pyximport
    pyximport.install()
    import predict_alleles_utils_fast
    predict_alleles_utils = predict_alleles_utils_fast
    using_fast = True
except:
    import predict_alleles_utils


genotype_mapping = { '00': '5',
                     '11': '0',
                     '12': '1',
                     '21': '1',
                     '22': '2' }



parser= optparse.OptionParser(usage="usage: "
                              "%prog -i <HAPLOTYPE FILE> -g <LOW-DENSITY GENOTYPE FILE> [-v]")
parser.add_option("-i", "--haplotype-assignment",
                  action="store", dest="hapfile",
                  type="string", default=None,
                  help="The file containing the haplotype assignments (JSON format).",
                  metavar="FILE")
parser.add_option("-g", "--high-density-genotypes",
                  action="store", dest="highgenfile",
                  type="string", default=None,
                  help="The file containing the input high-density genotypes.",
                  metavar="FILE")
parser.add_option("--left-border",
                  action="store", dest="lborder",
                  type="int", default=-1,
                  help="The number of loci to discard at the haplotype left border.")
parser.add_option("--right-border",
                  action="store", dest="rborder",
                  type="int", default=-1,
                  help="The number of loci to discard at the haplotype right border.")
parser.add_option("--alpha",
                  action="store", dest="alpha",
                  type="float", default=0.05,
                  help="The dampening parameter alpha (default: %default).")
parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose",
                  default=False,
                  help="Print additional log messages.")

(options, args)= parser.parse_args()

log_level= logging.DEBUG if options.verbose else logging.INFO

logging.basicConfig(level=log_level,
                    format='%(levelname)-8s [%(asctime)s]  %(message)s',
                    datefmt="%y%m%d %H%M%S")


logging.info("Haplotype-based Low-density Genotype Imputation")
if using_fast:
    logging.info("Using fast computation routines.")
logging.debug("Using dampening parameter alpha=%f", options.alpha)

hapallelepred_utils.check_file(options.hapfile)
hapallelepred_utils.check_file(options.highgenfile)

logging.info("Reading haplotype assignments from file '%s'...", options.hapfile)
haplotype_assignment = {}
with gzip.open(options.hapfile, "rb") as hf:
    haplotype_assignment = json.load(hf)

logging.info("Processing the haplotype assignments...")
haps = haplotype_assignment["haplotypes"]
no_orig_ind = sum(h["frequency"] for h in haps.values())

def combine_haplotypes(h1, h2):
    return "".join(genotype_mapping[h1i+h2i] for h1i,h2i in zip(h1,h2))

pgens = {}
for h1,h2 in hapallelepred_utils.lower_triangle(haps.keys()):
#    logging.debug("Considering haplotype pair:")
#    logging.debug("H1:  %s   allele=%s freq=%d",
#                  h1, haps[h1]["assigned_allele"], haps[h1]["frequency"])
#    logging.debug("H2:  %s   allele=%s freq=%d",
#                  h2, haps[h2]["assigned_allele"], haps[h2]["frequency"])
    gen = combine_haplotypes(h1, h2)
    newgen = { "haplotype1": h1,
               "haplotype2": h2,
               "frequency":
                   float(haps[h1]["frequency"]*haps[h2]["frequency"])/(no_orig_ind*no_orig_ind),
               "ld_gen": "".join(sorted(haps[h1]["assigned_allele"]+haps[h2]["assigned_allele"])),
              }
    if gen not in pgens:
        pgens[gen] = {}
    if newgen["ld_gen"] not in pgens[gen]:
        pgens[gen][newgen["ld_gen"]] = newgen
    else:
        pgens[gen][newgen["ld_gen"]]["frequency"] += newgen["frequency"]

npgens = {}
for gen in pgens.keys():
    if len(pgens[gen])==1:
        npgens[gen] = pgens[gen].values()[0]
    else:
        logging.debug("Conflicting associations for genotype %s.", gen)
        # DISCARD
        logging.debug("Discarding...")
        # KEEP
        #associations = sorted(pgens[gen].values(), key= lambda x: -x["frequency"])
        #logging.debug(" given associations: %s", ", ".join([ "{0} [{1:.7f}]".format(x["ld_gen"], x["frequency"])
        #                                                     for x in associations]))
        #if associations[0]["frequency"]>2*associations[1]["frequency"]:
        #    npgens[gen] = associations[0]
        #    logging.debug("Keeping that (%s) with largest frequency %.7f.",
        #                  npgens[gen]["ld_gen"], npgens[gen]["frequency"])
        #else:
        #    logging.debug("The two most frequent associations have the same frequency. Discard the genotype...")
        #del associations

#for gen in sorted(pgens.keys(), key=lambda x: -pgens[x]["frequency"]):
#    logging.debug("Genotype: %s  frequency=%.7f ld_gen=%s", gen, pgens[gen]["frequency"], pgens[gen]["ld_gen"])

del pgens
pgens = { tuple([int(x) for x in gen ]):npgens[gen] for gen in npgens }


logging.info("Reading high-density genotypes from file '%s'...",
             options.highgenfile)
hdgen = {}
with open(options.highgenfile, "r") as gf:
    for line in gf:
        if line.startswith("#"):
            continue
        line = line.strip()
        (iid, gen) = line.split()
        hdgen[iid] = gen

logging.info("Read %d high-density genotypes.", len(hdgen))


logging.info("Processing the high-density genotypes...")


for iid,gen in sorted(hdgen.items(), key=lambda (iid,gen): iid):
#    logging.debug("Processing individual '%s'...", iid)
#    logging.debug("Obs. G:    %s", gen)
    genv = tuple([int(x) for x in gen])
    possible_ld_gens = predict_alleles_utils.possible_associations(genv, pgens, options.alpha)
    preds = sorted(possible_ld_gens, key=lambda x: -possible_ld_gens[x]["frequency"])
    if preds:
        prediction = { "ld_gen": preds[0],
                       "frequency": possible_ld_gens[preds[0]]["frequency"],
                       "log_odd_ratio": ( ( math.log10(possible_ld_gens[preds[0]]["frequency"]) -
                                            math.log10(possible_ld_gens[preds[1]]["frequency"]) )
                                          if len(possible_ld_gens)>1 else float("nan") ) }
        logging.debug("%10s %s   frequency=%.7f  log_odd_ratio=%.3f",
                      iid,
                      prediction["ld_gen"],
                      prediction["frequency"],
                      prediction["log_odd_ratio"])
        print("{0} {1} {2:.9f} {3:.2f}".format(iid,
                                               prediction["ld_gen"],
                                               prediction["frequency"],
                                               prediction["log_odd_ratio"])
              )
    else:
        logging.warning("No prediction for individual '%s'.", iid)

