if [ $# -eq 0 ]; then
    echo "path to file is required. To run this script: sh killjobs.sh PATHTOFILE"
    exit 1
fi

awk 'BEGIN { FS = ":" } ; {print $2}' $1 > killSelectedjobs.txt

while read line
do
qdel $line

done<killSelectedjobs.txt

echo "jobs killed!"
rm killSelectedjobs.txt
