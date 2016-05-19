#!/bin/bash -l

# The first argument ($1) must specify the root directory of the MWACUtils source code
# The second argument ($2) must specify a valid remote branch name

function usage {
  echo "usage: $0 [MWACUtils source code directory] [remote branch]"
  exit 1
}

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
install_dir=${src_dir}/test/${branch}

# Setup log file
logdir="${src_dir}/install/log"
mkdir -p "$logdir"
logfile="${logdir}/import.${branch}.log"

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

# Build it!
echo "Building \"${branch}\"..." | tee -a ${logfile}
echo "  ./build.sh \"${src_dir}\" \"${install_dir}\"" | tee -a ${logfile}
if ! ./build_mwacutils.sh "${src_dir}" "${install_dir}" >> "${logfile}" 2>&1
then
  echo "  build.sh failed. See logs for details."
  exit 1
fi

#version=`${install_dir}/bin/make_beam -V`

# Test it!
echo "Testing \"${branch}\"..." | tee -a ${logfile}
echo "  (Automatic testing not implemented. Quitting...)" | tee -a ${logfile}
exit 0

#
# THIS SCRIPT SHOULD HAVE EXITED BY NOW. THE FOLLOWING CODE WILL NOT BE RUN (UNTIL TESTING IS IMPLEMENTED).
#

# Deploy it!
mkdir -p "${install_dir}"
./build.sh "${install_dir}"

# Checkout the master branch again (we were never here...)
git checkout master
