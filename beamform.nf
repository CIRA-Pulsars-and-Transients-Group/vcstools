params.obsid = null
params.pointings = null
params.calid = null

params.begin = null
params.end = null

params.basedir = '/group/mwaops/vcs'
params.didir = "${params.basedir}/${params.obsid}/cal/${params.calid}/rts"
params.channels = null

process gps_to_utc {
    output:
    stdout utctime

    """
    #!/usr/bin/env python

    from process_vcs import gps_to_utc
    print(gps_to_utc($params.begin))
    """
}

process beamform{
    input:
    val channel from channels
    val utc from utctime

    //TODO add other beamform options and flags -F
    """
    make_beam -o $params.obsid -b $params.begin -e $params.end -a 128 -n 128 \
    -f $channel -J ${params.didir}/DI_JonesMatrices_node{1:0>3}.dat \
    -d ${params.basedir}/${params.obsid}/combined -P $params.pointings \
    -r 10000 -m ${params.basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits \
    -z $utc -p
    """
}