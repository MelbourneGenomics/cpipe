# vim: ts=4:expandtab:sw=4:cindent
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
import os
import unittest
import csv

class QCExcelReportTest(unittest.TestCase):
    
    def __init__(self,n):
        unittest.TestCase.__init__(self,n)

    def run_script( self, s ):
        os.system( 'JAVA_OPTS="-Xmx4g -Djava.awt.headless=true" %s -cp %s/groovy-ngs-utils.jar:%s/excel.jar ./%s' % ( os.environ['GROOVY'], os.environ['GROOVY_NGS'], os.environ['EXCEL'], s ) )

    def testSimple(self):
        self.run_script( '../scripts/qc_excel_report.groovy -s 012345678 -o qctest -p qcprefix -resultsdir . ./data/012345678.*' )
 
if __name__ == '__main__':
    unittest.main()             
