#!/usr/bin/env bash

# This is a wrapper script designed to make your like exceptionally
# convenient when checking your personal or group disk usage and quotas.
#
# It is important to note though, that the underlying commands are
# easy to handle themselves, thus we recommend users become familiar
# with them, as they are handy to know for all Lustre Filesystems.

# If the user wants to provide a specific partition, then allow it,
# otherwise default to checking the /astro mount.
partition=${1:-/astro}

# Disk space used on /group and /astro by $USER, in human-readable formats
echo "==== Disk usage for user: ${USER} on /astro ===="
lfs quota -h -u "$USER" "$partition" | head -n3 | tail -n+2
echo

echo "==== Disk usage for group: mwavcs on /astro ===="
lfs quota -h -g mwavcs "$partition" | head -n3 | tail -n+2
echo
