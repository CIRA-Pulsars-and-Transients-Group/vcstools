#!/bin/bash

groupUsage=`lfs quota -u $USER /group | tail -1 | awk '{print $2/1e9}'`
astroUsage=`lfs quota -u $USER /astro | tail -1 | awk '{print $2/1e9}'`
echo "=== Stats for user ${USER} ==="
echo "You're using ${groupUsage} TB on /group"
echo "You're using ${astroUsage} TB on /astro"

echo ""

groupUsage=`lfs quota -g mwaops /group | tail -1 | awk '{print $2/1e9}'`
astroUsage=`lfs quota -g mwaops /astro | tail -1 | awk '{print $2/1e9}'`
echo "=== Stats for the group mwaops ==="
echo "${groupUsage} TB on /group"
echo "${astroUsage} TB on /astro"
