#!/usr/bin/env bash

# NOTE: lfs quota returns the quota in kilobytes...
#       To convert this to terabytes: multiply by 1024, divide by 1000^4


# Disk space used on /group and /astro by $USER, in terabytes
echo "==== Disk usage for user: ${USER} ===="
if [[ $HOSTNAME != garrawarla* ]]; then 
    groupUsage=$(lfs quota -q -u $USER /group | awk '{print $2*1024/1000000000000}')
    echo "  Used on /group: $(printf "%7.3f" "$groupUsage") TB"
fi
astroUsage=$(lfs quota -q -u $USER /astro | awk '{print $2*1024/1000000000000}')
echo "  Used on /astro: $(printf "%7.3f" "$astroUsage") TB"
echo

echo "==== Total disk usage for group: mwavcs ===="

if [[ $HOSTNAME != garrawarla* ]]; then 
    # Quotas for /group in terabytes
    groupQuota=$(lfs quota -q -g mwavcs /group | awk '{print $3*1024/1000000000000}')
    # Disk space used on /group for the mwavcs group, in terabytes
    groupUsage=$(lfs quota -q -g mwavcs /group | awk '{print $2*1024/1000000000000}')
    # Percentage of quota filled
    groupPercent=$(printf "%.2f" $(bc -l <<< "100*${groupUsage}/${groupQuota}"))
    echo "  Used on /group: $(printf "%7.3f / %7.3f TB (%.2f%%)" ${groupUsage} ${groupQuota} ${groupPercent})"
fi

# Quotas for /astro
astroQuota=$(lfs quota -q -g mwavcs /astro | awk '{print $3*1024/1000000000000}')
# Disk space used on /astro for the mwavcs group, in terabytes
astroUsage=$(lfs quota -q -g mwavcs /astro | awk '{print $2*1024/1000000000000}')
# Percentage of quota filled
astroPercent=$(printf "%.2f" $(bc -l <<< "100*${astroUsage}/${astroQuota}"))
echo "  Used on /astro: $(printf "%7.3f / %7.3f TB (%.2f%%)" ${astroUsage} ${astroQuota} ${astroPercent})"
echo
