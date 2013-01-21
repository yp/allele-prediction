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

genotype_difference = [
    # 0
    #   0 1 2 3 4 5
      [ 0,1,2,0,0,0 ],
    # 1
    #   0 1 2 3 4 5
      [ 1,0,1,0,0,0 ],
    # 2
    #   0 1 2 3 4 5
      [ 2,1,0,0,0,0 ],
    # 3
    [],
    # 4
    [],
    # 5
    #   0 1 2 3 4 5
      [ 0,0,0,0,0,0 ]
]    


cdef int generalized_hamming_distance(tuple s1, tuple s2, gd=genotype_difference):
     cdef int i
     cdef int s = 0
     cdef int N = len(s1)
     s = 0
     for i in range(N):
          s += gd[s1[i]][s2[i]]
     return s


def possible_associations(genv, pgens, alpha=0.05):
    possible_ld_gens = collections.defaultdict(lambda: {"frequency": 0.0})
    for pgen in pgens:
        hd = generalized_hamming_distance(genv, pgen)
        if (hd > 20):
            continue
        possible_ld_gens[pgens[pgen]["ld_gen"]]["frequency"] += pgens[pgen]["frequency"] * (alpha ** hd)
    return possible_ld_gens
