#!/bin/bash
#
# Script for syncing:
#  * Installed tools
#  * Deployed resources.
# from primary install/deploy location (${SOURCE_ROOT_DIR}) 
# to tmp file systems (${DESTINATION_ROOT_DIRS[@]}) 
# in the GCC cluster environment.
#


#
##
### Functions.
##
#
function showHelp() {
  #
  # Display commandline help on STDOUT.
  #
  cat <<EOH

Usage:

   $(basename $0) [-l] -a
   $(basename $0) [-l] -r relative/path/to/resource/
   $(basename $0) [-l] -t toolname/toolversion

Details:

 -l   List: Do not perform actual sync, but only list changes instead (dry-run).

 -a   All: syncs all tools and resources from ${SOURCE_ROOT_DIR} to ${DESTINATION_ROOT_DIRS[@]}

 -r   Resource: syncs only the specified resource to ${DESTINATION_ROOT_DIRS[@]}.
      The specified path must be relative to ${SOURCE_ROOT_DIR}resources/

 -t   Tool: syncs only the specified tool to ${DESTINATION_ROOT_DIRS[@]}.
      The tool must be deployed as "module" and specified using name/version as per "module" command syntax.

EOH
  #
  # Clean up.
  #
  rm -Rf ${TMP_DIR}
  #
  # Reset trap and exit.
  #
  trap - EXIT
  exit 0
}

function reportError() {
  local SCRIPT_NAME=$(basename $0 .bash)
  local PROBLEMATIC_LINE=$1
  local exit_status=${2:-$?}
  local ERROR_MESSAGE=$(cat ${TMP_LOG} 2> /dev/null) || true
  local ERROR_MESSAGE=${ERROR_MESSAGE:-Unknown error.}
  local errorMessage=${3:-"${ERROR_MESSAGE}"}
  local LOG2STDERR=${LOG2STDERR:-1}
  if [ "${LOG2STDERR}" -eq 0 ]; then
    LOG2STDERR=' '
  else
    LOG2STDERR='-s'
  fi
  local SOURCE="${RSYNC_SOURCES[@]:-"${SOURCE_ROOT_DIR}"}"
  local DESTINATION="${RSYNC_DESTINATION:-"${DESTINATION_ROOT_DIRS[@]}"}"
  local DETAILED_LOGS="$(cat ${RSYNC_LOG} 2> /dev/null)" || true
  #
  # Notify syslog.
  #
  logger ${LOG2STDERR} "$(hostname) - ${SCRIPT_NAME}:${PROBLEMATIC_LINE}: FATAL: rsync of source(s) (${SOURCE}) to destination(s) (${DESTINATION}) FAILED!"
  logger ${LOG2STDERR} "$(hostname) - ${SCRIPT_NAME}:${PROBLEMATIC_LINE}: Exit code = $exit_status"
  logger ${LOG2STDERR} "$(hostname) - ${SCRIPT_NAME}:${PROBLEMATIC_LINE}: Error message = ${errorMessage}"
  logger ${LOG2STDERR} "$(hostname) - ${SCRIPT_NAME}:${PROBLEMATIC_LINE}: Details = ${DETAILED_LOGS:-none.}"
  #
  # Notify admins by e-mail.
  #
  echo "
$(hostname) - ${SCRIPT_NAME}:${PROBLEMATIC_LINE}: FATAL: rsync of ${SOURCE} to ${DESTINATION} FAILED!
$(hostname) - ${SCRIPT_NAME}:${PROBLEMATIC_LINE}: Exit code = $exit_status
$(hostname) - ${SCRIPT_NAME}:${PROBLEMATIC_LINE}: Error message = ${errorMessage}
===============================================================================
${DETAILED_LOGS:-}
" | \
  mail -s "rsync of ${SOURCE} to ${DESTINATION} FAILED!" \
       -r "${EMAIL_FROM}" \
       "${EMAIL_TO}" \
       || logger ${LOG2STDERR} "$(hostname) - ${SCRIPT_NAME}:${LINENO}: FATAL: Could not send email."
  #
  # Clean up.
  #
  rm -Rf ${TMP_DIR}
  #
  # Reset trap and exit.
  #
  trap - EXIT
  exit $exit_status
}

function createConfigTemplate () {
  (cat > "${SCRIPT_CONFIG}.template"  <<EOCT 

##########################################################
# Configuration file for the ${SCRIPT_NAME} script.
#
#   * Listing variables in bash syntax
#   * To activate this config:
#     * Edit this file and 
#     * Remove the .template suffix from the filename  
#
##########################################################

#
# System account, group and perms for tools + resources.
#
#  * Group on SOURCE will be recursively changed to this one before sync.
#  * These permissions will be applied recursively on SOURCE before sync.
#
#SYS_USER='envsync'
#SYS_GROUP='depad'
#SYS_FILE_PERMS='g+rwX,o+rX,o-w'
#SYS_FOLDER_PERMS='g+rwXs,o+rX,o-w'

#
# Original location where we deployed our tools + resources.
#
#SOURCE_ROOT_DIR='/gcc/'

#
# Locations on tmp* file system where we want a copy of our deployed tools + resources.
#
#declare -a DESTINATION_ROOT_DIRS=('/gcc/tmp01/' '/gcc/tmp02/' '/gcc/tmp03/')

#
# Should the script delete old stuff in DESTINATION when it is no longer present in SOURCE?
#
DELETE_OLD=0

#
# Errors are always logged to syslog.
# Errors are in addition logged to STDERR (default).
# The latter can be disabled by setting LOG2STDERR to 0.
#
LOG2STDERR=1

#
# Email reporting of failures.
#
#EMAIL_FROM='sysop.gcc.groningen@gmail.com'
#EMAIL_TO='gcc-groningen@googlegroups.com'

EOCT
) || {
    echo "FATAL: Cannot find/access ${SCRIPT_CONFIG} and could not create a template config file with disabled options either."
    trap - EXIT
    exit 1
  }
}

#
##
### Bash sanity and error trapping.
##
#

#
# Bash sanity.
#
set -u
set -e
umask 0027

#
# Trap all exit signals: HUP(1), INT(2), QUIT(3), TERM(15), ERR
#
trap 'reportError $LINENO' HUP INT QUIT TERM EXIT ERR

#
##
### Configure rsync job.
##
#

#
# Get path to directory where this script is located.
#
SCRIPT_DIR=$(cd -P "$( dirname "$0" )" && pwd)
SCRIPT_NAME=$(basename "$0" .bash)

#
# Script's config must be in the same location.
#
SCRIPT_CONFIG="${SCRIPT_DIR}/${SCRIPT_NAME}.cfg"
if [[ -r "${SCRIPT_CONFIG}" && -f "${SCRIPT_CONFIG}" ]]; then
  source "${SCRIPT_CONFIG}" || reportError ${LINENO} $? "Cannot source ${SCRIPT_CONFIG}."
else
  createConfigTemplate
  logger -s "$(hostname) - ${SCRIPT_NAME}:${LINENO}: FATAL: Cannot find/access ${SCRIPT_CONFIG}!"
  logger -s "$(hostname) - ${SCRIPT_NAME}:${LINENO}: INFO:  Created a template config file with disabled options: ${SCRIPT_CONFIG}.template."
  logger -s "$(hostname) - ${SCRIPT_NAME}:${LINENO}: INFO:  Edit + rename template and try again."
  trap - EXIT
  exit 1
fi

START_TS=$(date "+%Y-%m-%d-T%H%M")
TMP_DIR="${TMPDIR:-/tmp}/${SCRIPT_NAME}/"
#echo "DEBUG: Using TMP_DIR: ${TMP_DIR}."
RSYNC_LOG="${TMP_DIR}/${SCRIPT_NAME}-${START_TS}.log"
TMP_LOG="${TMP_DIR}/tmp-${START_TS}.log"

#
# Create tmp dir.
#
mkdir -p "${TMP_DIR}/"       || reportError ${LINENO} $? "Cannot create ${TMP_DIR}."
test -d "${TMP_DIR}"         || reportError ${LINENO} $? "Cannot access ${TMP_DIR}."
touch ${TMP_LOG}             || reportError ${LINENO} $? "Cannot create ${TMP_LOG}."

#
# Initialise empty rsync log file, so emailing the logs won't fail, because the log does not yet exist.
#
touch ${RSYNC_LOG} 2> ${TMP_LOG}

#
##
### Process commandline arguments.
##
#

#
# Get commandline arguments.
#
ALL=0
RESOURCE=0
TOOL=0
LIST=0
SOURCE=''
while getopts ":halr:t:" opt; do
  case $opt in
    h)
      showHelp
      ;;
    a)
      ALL=1
      ;;
    l)
      LIST=1
      ;;
    r)
      RESOURCE=1
      SOURCE="${OPTARG}"
      ;;
    t)
      TOOL=1
      SOURCE="${OPTARG}"
      ;;
    \?)
      reportError ${LINENO} '1' "Invalid option -${OPTARG}. Try \"$(basename $0) -h\" for help."
      ;;
    :)
      reportError ${LINENO} '1' "Option -${OPTARG} requires an argument. Try \"$(basename $0) -h\" for help."
      ;;
  esac
done

#
# Make sure there are no extra arguments we did not expect nor need.
#
shift $(($OPTIND - 1))
if [ ! -z ${1:-} ]; then
  reportError ${LINENO} '1' "Invalid argument \"$1\". Try \"$(basename $0) -h\" for help."
fi

#
# Check commandline arguments.
#
ARG_SUM=$((${ALL}+${RESOURCE}+${TOOL}))

if [ "${ARG_SUM}" -eq 0 ]; then
  #
  # No commandline arguments specified.
  #
  showHelp
elif [ "${ARG_SUM}" -gt 1 ]; then
  reportError ${LINENO} '1' "Too many mutually exclusive arguments specified. Try \"$(basename $0) -h\" for help."
fi

#
##
### Create ToDo list.
##
#
if [ ${ALL} -eq 1 ]; then
  #
  # Add all tools, modules and resources to list of data to rsync.
  #
  RSYNC_SOURCES[0]='tools'
  RSYNC_SOURCES[1]='modules'
  RSYNC_SOURCES[2]='resources'
elif [ ${RESOURCE} -eq 1 ]; then
  #
  # Find and add only specified resource to list of data to rsync.
  #
  cd "${SOURCE_ROOT_DIR}/resources/" 2> ${TMP_LOG} || reportError ${LINENO} $?
  if [ -e ${SOURCE} ]; then
    echo "INFO: Found resource ${SOURCE}."
  else
    reportError ${LINENO} $? "Cannot find resource ${SOURCE} in ${SOURCE_ROOT_DIR}/resources/."
  fi
  # Create list of RSYNC SOURCES
  RSYNC_SOURCES[0]="resources/${SOURCE}"
elif [ ${TOOL} -eq 1 ]; then
  #
  # Find and add only specified tool to list of data to rsync.
  #
  cd "${SOURCE_ROOT_DIR}/tools/" 2> ${TMP_LOG}     || reportError ${LINENO} $?
  MODULE_TOOL_SPEC=(${SOURCE//\// })               || reportError ${LINENO} $?
  if [ ${#MODULE_TOOL_SPEC[@]} -ne 2 ]; then
    reportError ${LINENO} $? "Illegal tool specification ${SOURCE}. Tool spec must be in format TOOL_NAME/TOOL_VERSION."
  fi
  TOOL_NAME=${MODULE_TOOL_SPEC[0]}
  #echo "BEDUG: TOOL_NAME    = ${TOOL_NAME}"
  TOOL_VERSION=${MODULE_TOOL_SPEC[1]}
  #echo "DEBUG: TOOL_VERSION = ${TOOL_VERSION}"
  VERSIONED_TOOL=$(ls -1 | grep ${TOOL_NAME}*${TOOL_VERSION})               || reportError ${LINENO} $? "Cannot find tool ${TOOL_NAME}*${TOOL_VERSION} in ${SOURCE_ROOT_DIR}/tools/."
  VERSIONED_TOOL_COUNT=$(ls -1 | grep ${TOOL_NAME}*${TOOL_VERSION} | wc -l) || reportError ${LINENO} $? "Cannot count instances of tool ${TOOL_NAME}*${TOOL_VERSION} in ${SOURCE_ROOT_DIR}/tools/."
  if [ ${VERSIONED_TOOL_COUNT} -ne 1 ]; then
    reportError ${LINENO} '1' "Found multiple directories in ${SOURCE_ROOT_DIR}/tools/ for ${TOOL_NAME}/${TOOL_VERSION}: ${VERSIONED_TOOL}"
  else
    echo "INFO: Found tool ${SOURCE}."
  fi
  # Create list of RSYNC SOURCES
  RSYNC_SOURCES[0]="tools/${VERSIONED_TOOL}"
  RSYNC_SOURCES[1]="modules/${TOOL_NAME}/${TOOL_VERSION}"
fi 

echo "INFO: RSYNC_SOURCES contains ${RSYNC_SOURCES[@]}"

#
# Define rsync options.
#
#RSYNC_OPTIONS='-avRK' Archive mode (-a) = -rlptgoD.
# We don't sync ownership of the files.
# Instead all secondary copies on the destinations are owned by ${SYS_USER}.
RSYNC_OPTIONS='-rlptgDvRK'
if [ "${DELETE_OLD}" -eq 1 ]; then
  echo "WARN: Cleanup of outdated ${SOURCE_ROOT_DIR} data is enabled for ${DESTINATION_ROOT_DIRS[@]}."
  RSYNC_OPTIONS="${RSYNC_OPTIONS} --delete-after"
fi
if [ "${LIST}" -eq 1 ]; then
  echo 'WARN: List mode enabled: will only list what is out of sync and needs to be updated, but will not perform actual sync.'
  RSYNC_OPTIONS="${RSYNC_OPTIONS} -nu"
else
  RSYNC_OPTIONS="${RSYNC_OPTIONS} -q"
fi

echo "INFO: RSYNC_OPTIONS contains ${RSYNC_OPTIONS}"

#
##
### Check environment.
##
#

#
# Check if we are running with the correct account + permissions.
#
CURRENT_USER=$(whoami)
if [ ${CURRENT_USER} != ${SYS_USER} ]; then
  reportError ${LINENO} '1' "This script must be executed by user ${SYS_USER}, but you are ${CURRENT_USER}."
fi
CURRENT_GROUP=$(id -gn)
if [ ${CURRENT_GROUP} != ${SYS_GROUP} ]; then
  reportError ${LINENO} '1' "This script must be executed by user ${SYS_USER} with primary group ${SYS_GROUP}, but your current primary group is ${CURRENT_GROUP}."
fi

#
# Recursively fix group + permissions on SOURCE (should not be necessary, but just in case :))
#
cd ${SOURCE_ROOT_DIR}
for (( i = 0 ; i < ${#RSYNC_SOURCES[@]:-0} ; i++ ))
do
  echo "INFO: Trying to fix group and permissions on ${SOURCE_ROOT_DIR}${RSYNC_SOURCES[${i}]} recursively before sync."
  echo '      Should not be necessary, but just in case...'
  echo "      This may fail (depending on current group and permissions) if user '${SYS_USER}' does not own the files/folders."
  #
  # We use find to try to fix group + perms only when they are not correct.
  # This prevents permission denied errors when there is no need to change group or perms and we do not own the files/folders. 
  #
  find "${RSYNC_SOURCES[${i}]}" ! -group "${SYS_GROUP}" -exec chgrp "${SYS_GROUP}" '{}' \;                     2> ${TMP_LOG} || reportError ${LINENO} $?
  find "${RSYNC_SOURCES[${i}]}" ! -type d ! -perm -${SYS_FILE_PERMS} -exec chmod "${SYS_FILE_PERMS}" '{}' \;   2> ${TMP_LOG} || reportError ${LINENO} $?
  find "${RSYNC_SOURCES[${i}]}" -type d ! -perm -${SYS_FOLDER_PERMS} -exec chmod "${SYS_FOLDER_PERMS}" '{}' \; 2> ${TMP_LOG} || reportError ${LINENO} $?
done

#
# Check if all destinations are available and remove destinations, which are offline!
# 
# This is critically essential as syncing to a mount point with missing mount would add the data to the disk containing the mount point, 
# which is usually a relatively small disk containing the OS. Running out of space on the local system disk, will crash a server!
#
declare -a AVAILABLE_DESTINATION_ROOT_DIRS
for (( i = 0 ; i < ${#DESTINATION_ROOT_DIRS[@]:-0} ; i++ ))
do 
  if [ -e ${DESTINATION_ROOT_DIRS[${i}]}/modules ] && [ -e ${DESTINATION_ROOT_DIRS[${i}]}/resources ] && [ -e ${DESTINATION_ROOT_DIRS[${i}]}/tools ] \
  && [ -r ${DESTINATION_ROOT_DIRS[${i}]}/modules ] && [ -r ${DESTINATION_ROOT_DIRS[${i}]}/resources ] && [ -r ${DESTINATION_ROOT_DIRS[${i}]}/tools ] \
  && [ -w ${DESTINATION_ROOT_DIRS[${i}]}/modules ] && [ -w ${DESTINATION_ROOT_DIRS[${i}]}/resources ] && [ -w ${DESTINATION_ROOT_DIRS[${i}]}/tools ]; then
    if [ "${#AVAILABLE_DESTINATION_ROOT_DIRS[@]:-0}" -eq 0 ]; then
      AVAILABLE_DESTINATION_ROOT_DIRS=("${DESTINATION_ROOT_DIRS[${i}]}")
    else
      AVAILABLE_DESTINATION_ROOT_DIRS=("${AVAILABLE_DESTINATION_ROOT_DIRS[@]:-}" "${DESTINATION_ROOT_DIRS[${i}]}")
    fi
  else
    echo "WARN: ${DESTINATION_ROOT_DIRS[${i}]} not available (symlink dead or mount missing). Skipping rsync to ${DESTINATION_ROOT_DIRS[${i}]}."
  fi
done

if [ "${#AVAILABLE_DESTINATION_ROOT_DIRS[@]:-0}" -gt 0 ]; then
  echo "INFO: AVAILABLE_DESTINATION_ROOT_DIRS contains ${AVAILABLE_DESTINATION_ROOT_DIRS[@]}"
else
  reportError ${LINENO} '1' "None of the destinations is available: Giving up!"
fi

#
##
### Rsync.
##
#

#
# Perform the rsync for all sources that need to be synced to all destinations.
#
cd ${SOURCE_ROOT_DIR}
for (( i = 0 ; i < ${#RSYNC_SOURCES[@]:-0} ; i++ ))
do
  for (( j = 0 ; j < ${#AVAILABLE_DESTINATION_ROOT_DIRS[@]:-0} ; j++ ))
  do
    echo "INFO: Rsyncing ${RSYNC_SOURCES[${i}]} to ${AVAILABLE_DESTINATION_ROOT_DIRS[${j}]}..."
    if [ "${LIST}" -eq 1 ]; then
      echo '================================================================================================================' >> "${RSYNC_LOG}"
      echo "  Dry run stats for syncing ${RSYNC_SOURCES[${i}]} to ${AVAILABLE_DESTINATION_ROOT_DIRS[${j}]}:" >> "${RSYNC_LOG}"
      echo '================================================================================================================' >> "${RSYNC_LOG}"
    fi
    set +e
    rsync ${RSYNC_OPTIONS} \
    "${RSYNC_SOURCES[${i}]}" \
    "${AVAILABLE_DESTINATION_ROOT_DIRS[${j}]}" \
    >> "${RSYNC_LOG}" 2>&1
    XVAL=$?
    set -e
    if [[ ${XVAL} -ne 0 && ${XVAL} -ne 24 ]]; then
      reportError ${LINENO} ${XVAL} "Rsync of sources (${RSYNC_SOURCES[@]}) to destinations (${AVAILABLE_DESTINATION_ROOT_DIRS[${j}]}) started on ${START_TS} failed."
    fi
  done
done

#
##
### Sanity check.
##
#

#
# Parse log: rsync log should exist and should be empty.
#
if [ "${LIST}" -eq 1 ]; then
  cat "${RSYNC_LOG}" || reportError ${LINENO} $? "Listing differences between sources (${RSYNC_SOURCES[@]}) and destinations (${AVAILABLE_DESTINATION_ROOT_DIRS[@]}) started on ${START_TS} failed: cannot display ${RSYNC_LOG} contents!"
elif [[ ! -f "${RSYNC_LOG}" || -s "${RSYNC_LOG}" ]]; then
  reportError ${LINENO} $? "Rsync of sources (${RSYNC_SOURCES[@]}) to destinations (${AVAILABLE_DESTINATION_ROOT_DIRS[@]}) started on ${START_TS} failed: error log not empty!"
fi

#
# Cleanup.
#
(rm ${TMP_LOG} ; rm ${RSYNC_LOG} ; rmdir ${TMP_DIR}) || reportError ${LINENO} $? "Cannot cleanup tmp dir ${TMP_DIR}."

#
# Signal succes.
#
echo "INFO: Rsync finished successfully for sources (${RSYNC_SOURCES[@]}) to destinations (${AVAILABLE_DESTINATION_ROOT_DIRS[@]})."

#
# Reset trap and exit.
#
trap - EXIT
exit 0
