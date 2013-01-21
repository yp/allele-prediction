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

import collections
import logging
import os
import os.path
import sys

def check_file(filename):
    if ( filename is None
         or not os.path.isfile(filename)
         or not os.access(filename, os.R_OK) ):
        logging.fatal("Input file not specified or invalid. Given: '%s'.", filename)
        sys.exit("Input file '{0}' not valid.".format(filename))


def lower_triangle(n):
    v = None
    if not isinstance(n, collections.Iterable):
        v = list(range(n))
    else:
        v = list(n)

    for i2 in xrange(1, len(v)):
        for i1 in xrange(i2):
            yield (v[i1], v[i2])


