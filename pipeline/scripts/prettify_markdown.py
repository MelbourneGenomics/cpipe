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

import sys

sys.stdout.write( '<html>\n<head>\n<link rel="stylesheet" href="http://yui.yahooapis.com/pure/0.6.0/pure-min.css">\n</head>\n<body>\n<div class="pure-g">\n<div class="pure-u-1">\n' )

for line in sys.stdin:
  line = line.replace( '<table>', '<table class="pure-table pure-table-bordered">' )
  sys.stdout.write( line )

sys.stdout.write( '</div>\n</div>\n</body>\n</html>\n' )

