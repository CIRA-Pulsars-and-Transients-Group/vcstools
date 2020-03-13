params.obsid = null
params.pointings = null
params.calid = null

params.begin = null
params.end = null
params.all = false

params.summed = false
params.vcstools_version = 'master'

params.basedir = '/group/mwaops/vcs'
params.didir = "${params.basedir}/${params.obsid}/cal/${params.calid}/rts"
params.channels = null


//----------------

if ( params.summed ) {
    bf_out = " -p -s "
}
else {
    bf_out = " -p "
}

pointings = Channel
    .from(params.pointings.split(","))
    .collect()
    .flatten()
    .collate( 15 )
    //.view()

// Handling begin and end times
process get_all_beg_end {
    when:
    params.all == true

    output:
    stdout obs_all

    """
    #!/usr/bin/env python3

    from mwa_metadb_utils import obs_max_min

    beg, end = obs_max_min(${params.obsid})
    print("{},{}".format(beg, end)),
    """
}

if ( params.all ) {
    obs_all
        .trim()
        .map { it.split(",") }
        .flatten()
        .collect()
        //.view()
        .into { obs_beg_end; obs_beg }
}
else {
    obs_temp = Channel
        .from( params.begin, params.end )
        .collect()
        //.view()
        .into { obs_beg_end; obs_beg }
}

//Working out how many output fits files there will be
total_time = params.end - params.begin + 1 
int nfiles = (int) total_time / 200
if ( total_time%200 != 0) {
    nfiles += 1
}
n_expected_splice_files = nfiles * 24

println n_expected_splice_files
exit 0

process ensure_metafits {

    """
    #!/usr/bin/env python3

    from process_vcs import ensure_metafits
    
    ensure_metafits("${params.basedir}/${params.obsid}", "${params.obsid}",\
                    "${params.basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits")
    """
}

//obs_beg.view()

process gps_to_utc {
    input:
    tuple val(begin), val(end) from obs_beg

    output:
    stdout utctime

    """
    #!/usr/bin/env python3

    from process_vcs import gps_to_utc

    print(gps_to_utc($begin)),
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
    label 'gpu'
    time '10h'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(channel), val(ch) from channels
    val utc from utctime
    each point from pointings
    tuple val(begin), val(end) from obs_beg_end

    output:
    file "*/*fits" into unspliced

    //TODO add other beamform options and flags -F
    """
    module use /group/mwa/software/modulefiles
    module load vcstools/$params.vcstools_version
    make_beam -o $params.obsid -b $begin -e $end -a 128 -n 128 \
-f $channel -J ${params.didir}/DI_JonesMatrices_node${ch}.dat \
-d ${params.basedir}/${params.obsid}/combined -P ${point.join(",")} \
-r 10000 -m ${params.basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits \
${bf_out} -z $utc
    """
}

unspliced
    .flatten()
    .map { it -> [it.baseName.split("ch")[0], it ] }
    //.view()
    .groupTuple(size:  n_expected_splice_files)
    .map { it -> it[1] }
    //.view()
    .set{ unspliced_files }

process splice {
    publishDir "${params.basedir}/${params.obsid}/pointings/${unspliced[0].baseName.split("_")[2]}_${unspliced[0].baseName.split("_")[3]}", mode: 'move'

    label 'cpu'
    time '1h'

    input:
    val chan from chanlist
    each file(unspliced) from unspliced_files

    output:
    file "${params.obsid}*fits" into output_fits
    """
    module use /group/mwa/software/modulefiles
    module load mwa_search
    module load vcstools
    splice_wrapper.py -o ${params.obsid} -c ${chan.join(" ")}
    """
}
