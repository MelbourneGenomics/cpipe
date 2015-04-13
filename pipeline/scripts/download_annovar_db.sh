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
# 
# This script downloads the necessary Annovar databases
#
###########################################################################

if [ -z "$1" ] || [ -z "$2" ];
then
        echo
        echo "Usage: download_annovar_db.sh <annovar dir> <db dir>"
        echo
        exit 1
fi

ANNOVAR="$1"
DBDIR="$2"

#$ANNOVAR/annotate_variation.pl -downdb ALL.sites.2010_11 -buildver hg19 "$DBDIR"
$ANNOVAR/annotate_variation.pl -downdb 1000g2010nov -buildver hg19  "$DBDIR"
$ANNOVAR/annotate_variation.pl -downdb snp138 -buildver hg19  "$DBDIR"
$ANNOVAR/annotate_variation.pl -downdb avsift -webfrom annovar -buildver hg19  "$DBDIR"

$ANNOVAR/annotate_variation.pl -downdb refGene -buildver hg19  "$DBDIR"
$ANNOVAR/annotate_variation.pl -downdb knownGene -buildver hg19  "$DBDIR"
$ANNOVAR/annotate_variation.pl -downdb genomicSuperDups -buildver hg19  "$DBDIR"
$ANNOVAR/annotate_variation.pl -downdb ljb_pp2 -webfrom annovar  -buildver hg19  "$DBDIR"
$ANNOVAR/annotate_variation.pl -downdb ljb_all -webfrom annovar  -buildver hg19  "$DBDIR"
$ANNOVAR/annotate_variation.pl -downdb esp5400_all   -buildver hg19 -webfrom annovar   "$DBDIR"
$ANNOVAR/annotate_variation.pl -downdb phastConsElements46way   -buildver hg19   "$DBDIR"
$ANNOVAR/annotate_variation.pl -buildver hg19  -downdb -webfrom annovar exac03  "$DBDIR"
