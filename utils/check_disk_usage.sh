#!/usr/bin/env bash

# NOTE: lfs quota returns the quota in kilobytes...
#       To convert this to terabytes: multiply by 1024, divide by 1000^4


# Disk space used on /group and /astro by $USER, in terabytes
groupUsage=$(lfs quota -q -u $USER /group | awk '{print $2*1024/1000000000000}')
astroUsage=$(lfs quota -q -u $USER /astro | awk '{print $2*1024/1000000000000}')
echo "==== Disk usage for user: ${USER} ===="
echo "  Used on /group: $(printf "%7.3f" "$groupUsage") TB"
echo "  Used on /astro: $(printf "%7.3f" "$astroUsage") TB"
echo


# Quotas for /group and /astro, in terabytes
groupQuota=$(lfs quota -q -g mwavcs /group | awk '{print $3*1024/1000000000000}')
astroQuota=$(lfs quota -q -g mwavcs /astro | awk '{print $3*1024/1000000000000}')

# Disk space used on /group and /astro for the mwavcs group, in terabytes
groupUsage=$(lfs quota -q -g mwavcs /group | awk '{print $2*1024/1000000000000}')
astroUsage=$(lfs quota -q -g mwavcs /astro | awk '{print $2*1024/1000000000000}')

# Percentage of quota filled
groupPercent=$(printf "%.2f" $(bc -l <<< "100*${groupUsage}/${groupQuota}"))
astroPercent=$(printf "%.2f" $(bc -l <<< "100*${astroUsage}/${astroQuota}"))
echo "==== Total disk usage for group: mwavcs ===="
echo "  Used on /group: $(printf "%7.3f / %7.3f TB (%.2f%%)" ${groupUsage} ${groupQuota} ${groupPercent})"
echo "  Used on /astro: $(printf "%7.3f / %7.3f TB (%.2f%%)" ${astroUsage} ${astroQuota} ${astroPercent})"
echo
