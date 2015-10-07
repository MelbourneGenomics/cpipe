#!/usr/bin/env python
###########################################################################
#
# This file is part of Cpipe.
# 
# Cpipe is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, under version 3 of the License, subject
# to additional terms compatible with the GNU General Public License version 3,
# specified in the LICENSE file that is part of the Cpipe distribution.
#
# Cpipe is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Cpipe.  If not, see <http:#www.gnu.org/licenses/>.
#
###########################################################################
# This script takes a sample metadata file as input (stdin) and
# prints out the field values in a human readable form
###########################################################################


import sys

headers = sys.stdin.readline().strip().split('\t')
for idx, line in enumerate(sys.stdin):
  print "===== Sample %i =====" % idx
  fields = line.strip('\n').split('\t')
  for jdx, field in enumerate(fields):
    print "%24s: %s" % ( headers[jdx], field )
