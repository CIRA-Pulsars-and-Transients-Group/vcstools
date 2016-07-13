#!/bin/bash -l

# The first argument ($1) must specify the root directory of the MWACUtils source code
# The second argument ($2) must specify a valid remote branch name
# The third optional argument ($3) has the following meanings, with regard to
# version numbers:
#    (absent) --> patch++
#       -v    --> minor++
#       -V    --> major++

deploy_install_dir=/group/mwaops/PULSAR

function usage {
  echo "usage: $0 <MWACUtils source code directory> <remote branch> [-v|-V]"
  echo "  -V = increment major version number,"
  echo "  -v = increment minor version number,"
  echo "  otherwise, increment patch version number"
  exit 1
}

# Parse arguments
if [ -z $1 ]
then
  usage
fi

if [ ! -d $1 ]
then
  echo "Cannot find directory \"$1\""
  usage
fi

src_dir=`readlink -f "$1"`
branch=$2
install_dir=${src_dir}/install/${branch}

# Setup log file
logdir=`readlink -f "../build_logs"`
mkdir -p "$logdir"
logfile="${logdir}/import_mwacutils.${branch}.log"

{
  echo "---------------------------------------------"
  echo "Running \"$0 $@\" on `date`:"
  echo
} >> ${logfile}

# Get it!
echo "Fetching remote repository \"${branch}\": " | tee -a ${logfile}
echo "  ./fetch_and_checkout.sh \"${src_dir}\" \"${branch}\"" | tee -a ${logfile}
if ! ./fetch_mwacutils.sh "${src_dir}" "${branch}" >> "${logfile}" 2>&1
then
  echo "  fetch_and_checkout.sh failed. See logs for details."
  exit 1
fi

# Get version of currently installed make_beam
#current_version=`make_beam -V`
current_version=`/group/mwaops/smcsweeney/bin/make_beam -V` # Temporary hack before -V system is built into the deployed version
major_version=`echo ${current_version} | cut -d. -f1`
minor_version=`echo ${current_version} | cut -d. -f2`
patch_version=`echo ${current_version} | cut -d. -f3`

# Increment version number and commit to git repo
if [ "$3" = "-V" ]
then
  major_version=`echo "${major_version} + 1" | bc`
elif [ "$3" = "-v" ]
then
  minor_version=`echo "${minor_version} + 1" | bc`
else
  patch_version=`echo "${patch_version} + 1" | bc`
fi
new_version="${major_version}.${minor_version}.${patch_version}"

# Build it!
{
  echo "Building \"${branch}\"..."
  echo "  (Current version: ${current_version}; Proposed version: ${new_version})"
  echo "  ./build_mwacutils.sh \"${src_dir}\" \"${install_dir}\""
} | tee -a ${logfile}

if ! ./build_mwacutils.sh "${src_dir}" "${install_dir}" >> "${logfile}" 2>&1
then
  echo "  build_mwacutils.sh failed. See logs for details."
  exit 1
fi

# Test it!
{
  echo "Testing \"${branch}\"..."
  echo "  ./test_mwacutils.sh \"${install_dir}\""
} | tee -a ${logfile}
if ! ./test_mwacutils.sh \"${install_dir}\" >> "${logfile}" 2>&1
then
  echo "  test_mwacutils.sh failed. See logs for details."
  exit 1
fi

#
# THIS SCRIPT SHOULD HAVE EXITED BY NOW. THE FOLLOWING CODE WILL NOT BE RUN (UNTIL TESTING IS IMPLEMENTED).
#

# Input new version into make_beam source code

#old_dir=`pwd`
#cd "${src_dir}"
#version_file="${src_dir}/mwac_utils/beamer_version.h"
#sed -i "s;\"[0-9\.]*\";\"${new_version}\";" "${version_file}"
#git add "${version_file}"
#git commit -m "Update to version ${new_version}"
#cd "${old_dir}"

# Deploy it!
#mkdir -p "${deploy_install_dir}"
#./build.sh \"${src_dir}\" "${deploy_install_dir}"

# Checkout the master branch again (we were never here...)
#git checkout master
