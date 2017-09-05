#!/bin/bash
groupUsage=`lfs quota -u $USER /group | tail -1 | awk '{print $2/1e9}'`
scratch2Usage=`lfs quota -u $USER /scratch2 | tail -1 | awk '{print $2/1e9}'`
astroUsage=`lfs quota -u $USER /astro | tail -1 | awk '{print $2/1e9}'`
echo "Stats for user ${USER}"
echo "You're using ${groupUsage} TB on /group"
echo "You're using ${scratch2Usage} TB on /scratch2"
echo "You're using ${astroUsage} TB on /astro"
groupUsage=`lfs quota -g mwaops /group | tail -1 | awk '{print $2/1e9}'`
scratch2Usage=`lfs quota -g mwaops /scratch2 | tail -1 | awk '{print $2/1e9}'`
astroUsage=`lfs quota -g mwaops /astro | tail -1 | awk '{print $2/1e9}'`
echo ""
echo "Stats for the group mwaops"
echo "${groupUsage} TB on /group"
echo "${scratch2Usage} TB on /scratch2"
echo "${astroUsage} TB on /astro"
