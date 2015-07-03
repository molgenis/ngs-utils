function usage () {
echo "
Arguments
        Required:
        -s|--stop              At which number to stop (double digits, e.g. 09)

        Optional:
        -b|--begin    		What is the startnumber (double digits) (default: 01)
"
}

PARSED_OPTIONS=$(getopt -n "$0"  -o b:s: --long "stop:,begin:"  -- "$@")

#
# Bad arguments, something has gone wrong with the getopt command.
#
if [ $? -ne 0 ]; then
        usage
        echo "FATAL: Wrong arguments."
        exit 1
fi

eval set -- "$PARSED_OPTIONS"

#
# Now goes through all the options with a case and using shift to analyse 1 argument at a time.
# $1 identifies the first argument, and when we use shift we discard the first argument, so $2 becomes $1 and goes again through the case.
#

while true; do
  case "$1" in
        -b|--begin)
                case "$2" in
                "") shift 2 ;;
                *) BEGIN=$2 ; shift 2 ;;
            esac ;;
        -s|--stop)
                case "$2" in
                *) STOP=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
  esac
done

THISDIR=`pwd`
echo "THIS DIRECTORY: $THISDIR"
#
# Check required options were provided.
if [[ -z "${STOP-}" ]]; then
        usage
        echo "FATAL: missing required parameter."
        exit 1
fi

if [[ -z "${BEGIN-}" ]]; then
        BEGIN="01"
fi

for i in $(printf "%02d " $(seq $BEGIN $STOP))

do
        ls -1 ${THISDIR}/s${i}*.sh >> ${THISDIR}/filenaam.txt

done

while read line

do
        touch $line.finished
	ENV=${line%???}
        touch $ENV.env

done<${THISDIR}/filenaam.txt

rm ${THISDIR}/filenaam.txt
