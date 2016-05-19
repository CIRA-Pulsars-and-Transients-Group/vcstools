#!/bin/bash -l

# Make sure argument $1 is an existing directory
src_dir="$1"
if [ ! -d ${src_dir} ]
then
  echo "Cannot find directory \"${src_dir}\""
  exit 1
fi

cd "${src_dir}"

# Make sure argument $2 is a valid remote branch
branch="$2"
if [ ! `git branch -r | grep "^ *${branch}/master$"` ]
then
  echo "Cannot find remote branch \"${branch}\" in git repository"
  echo "git branch -r:"
  git branch -r
  echo "Exiting..."
  exit 1
fi

# Fetch the branch and check it out
git fetch "${branch}"
git checkout "${branch}/master"
