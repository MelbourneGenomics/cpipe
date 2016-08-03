if [ -z "$BASE" ];
then
    BASE="."
fi
CONFIG=`sed 's/\/\/.*$//' $BASE/pipeline/config.groovy` 
eval "$CONFIG"
