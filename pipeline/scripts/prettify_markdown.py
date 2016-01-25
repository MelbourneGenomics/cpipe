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
###################################################################################
#
# Purpose:
#   Turn output from markdown2.py into something that looks better
# Usage:
#   prettify < html > html
####################################################################################

import re
import sys

sys.stdout.write( '<html>\n<head>\n<link rel="stylesheet" href="http://yui.yahooapis.com/pure/0.6.0/pure-min.css">\n</head>\n<body>\n<div class="pure-g">\n<div class="pure-u-1-24"></div><div class="pure-u-22-24">\n' )

in_thead = None
is_empty = True
for line in sys.stdin:
  line = line.replace( '<table>', '<table class="pure-table pure-table-bordered">' )
  # remove empty table headers
  if in_thead is not None: # already in thead
    if '</thead>' in line:
      if not is_empty:
        sys.stdout.write(in_thead)
        sys.stdout.write(line)
      else:
        pass # don't print empty
      in_thead = None
    else:
      # determine if line contains content
      if re.search('<th>[^<]+</th>', line) is not None:
        is_empty = False
      in_thead += line
  else:
    if '<thead>' in line:
      in_thead = line
      is_empty = True
    else:
      sys.stdout.write( line )

sys.stdout.write( '</div>\n<div class="pure-u-1-24"></div>\n</div>\n</body>\n</html>\n' )

