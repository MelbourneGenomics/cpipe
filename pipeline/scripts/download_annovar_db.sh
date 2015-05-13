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

function all_exist() {
    for i in $*;
    do
        if [ ! -e $DBDIR/$i ];
        then	
            return 1
        fi
        return 0
    done
}

ANNOVAR="$1"
DBDIR="$2"

ONEKG_FILES="hg19_AFR.sites.2014_10.txt hg19_AFR.sites.2014_10.txt.idx hg19_ALL.sites.2014_10.txt hg19_ALL.sites.2014_10.txt.idx hg19_AMR.sites.2014_10.txt hg19_AMR.sites.2014_10.txt.idx hg19_EAS.sites.2014_10.txt hg19_EAS.sites.2014_10.txt.idx hg19_EUR.sites.2014_10.txt hg19_EUR.sites.2014_10.txt.idx hg19_SAS.sites.2014_10.txt hg19_SAS.sites.2014_10.txt.idx"

#$ANNOVAR/annotate_variation.pl -downdb ALL.sites.2010_11 -buildver hg19 "$DBDIR"
#$ANNOVAR/annotate_variation.pl -downdb 1000g2010nov -buildver hg19  "$DBDIR"

all_exist $ONEKG_FILES || {
    $ANNOVAR/annotate_variation.pl -downdb 1000g2014oct -buildver hg19  "$DBDIR"
}

all_exist hg19_snp138.txt || {
    $ANNOVAR/annotate_variation.pl -downdb snp138 -buildver hg19  "$DBDIR"
}

#$ANNOVAR/annotate_variation.pl -downdb avsift -webfrom annovar -buildver hg19  "$DBDIR"
all_exist hg19_refGeneMrna.fa  hg19_refGene.txt  hg19_refGene.txt.gz || {
    $ANNOVAR/annotate_variation.pl -downdb refGene -buildver hg19  "$DBDIR"
}

#$ANNOVAR/annotate_variation.pl -downdb knownGene -buildver hg19  "$DBDIR"
$ANNOVAR/annotate_variation.pl -downdb genomicSuperDups -buildver hg19  "$DBDIR"
$ANNOVAR/annotate_variation.pl -downdb ljb_pp2 -webfrom annovar  -buildver hg19  "$DBDIR"

$ANNOVAR/annotate_variation.pl -downdb ljb26_all -webfrom annovar  -buildver hg19  "$DBDIR"

all_exist hg19_esp6500siv2_all.txt || {
    $ANNOVAR/annotate_variation.pl -downdb esp6500siv2_all   -buildver hg19 -webfrom annovar   "$DBDIR"
}

all_exist hg19_phastConsElements46way.txt || {
    $ANNOVAR/annotate_variation.pl -downdb phastConsElements46way   -buildver hg19   "$DBDIR"
}

all_exist hg19_exac03.txt  hg19_exac03.txt.idx || {
    $ANNOVAR/annotate_variation.pl -buildver hg19  -downdb -webfrom annovar exac03  "$DBDIR"
}
