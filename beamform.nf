params.obsid = null
params.pointings = null
params.calid = null

params.begin = null
params.end = null

params.basedir = '/group/mwaops/vcs'
params.didir = "${params.basedir}/${params.obsid}/cal/${params.calid}/rts"
params.channels = null


pointings = Channel
    .from(params.pointings.split(","))
    .collect()
    .flatten()
    .collate( 15 )
    .view()


process ensure_metafits {

    """
    #!/usr/bin/env python3

    from process_vcs import ensure_metafits
    
    ensure_metafits("${params.basedir}/${params.obsid}", "${params.obsid}",\
                    "${params.basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits")
    """
}
process gps_to_utc {
    output:
    stdout utctime

    """
    #!/usr/bin/env python3

    from process_vcs import gps_to_utc

    print(gps_to_utc($params.begin)),
    """
}

process get_channels {
    output:
    file "${params.obsid}_channels.txt" into channelsfile

    """
    #!/usr/bin/env python3

    from mwa_metadb_utils import get_channels
    import csv

    channels = get_channels($params.obsid)
    with open("${params.obsid}_channels.txt", "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        spamwriter.writerow(channels)
    """
}

range = Channel.from( ['001', '002', '003', '004', '005', '006',\
                       '007', '008', '009', '010', '011', '012',\
                       '013', '014', '015', '016', '017', '018',\
                       '019', '020', '021', '022', '023', '024'] )

channelsfile
    .splitCsv()
    .into{ chanlist; chantemp }

chantemp
    .flatten()
    .merge(range)
    //.view()
    .set{channels}


process beamform {
    executor 'slurm'
    cpus 3
    queue 'gpuq'
    time '1h'
    memory '12 GB'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(channel), val(ch) from channels
    val utc from utctime
    each point from pointings

    output:
    file "*/*fits" into unspliced

    //TODO add other beamform options and flags -F
    """
    make_beam -o $params.obsid -b $params.begin -e $params.end -a 128 -n 128 \
-f $channel -J ${params.didir}/DI_JonesMatrices_node${ch}.dat \
-d ${params.basedir}/${params.obsid}/combined -P ${point.join(",")} \
-r 10000 -m ${params.basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits \
-p -z $utc
    """
}

unspliced
    .flatten()
    .map { it -> [it.baseName.split("ch")[0], it ] }
    .groupTuple(size: 24)
    .map { it -> it[1] }
    //.view()
    .set{ unspliced_files }

process splice {
    publishDir "${params.basedir}/${params.obsid}/pointings/${unspliced[0].baseName.split("_")[2]}_${unspliced[0].baseName.split("_")[3]}", mode: 'move'

    executor 'slurm'
    cpus 3
    queue 'gpuq'
    time '1h'
    memory '4 GB'

    input:
    val chan from chanlist
    each file(unspliced) from unspliced_files

    output:
    file "${params.obsid}*fits" into output_fits
    """
    splice_wrapper.py -o ${params.obsid} -c ${chan.join(" ")}
    """
}