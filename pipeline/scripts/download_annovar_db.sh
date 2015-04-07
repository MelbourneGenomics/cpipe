#######################################################
# 
# This script downloads the necessary Annovar databases
#
#######################################################

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
