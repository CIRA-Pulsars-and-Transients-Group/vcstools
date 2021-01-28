"""
Functions to handle parsing the config file for multiple super computers
"""

#config data

GALAXY_CONFIG = {'base_data_dir' : '/astro/mwavcs/vcs/',
                 'base_product_dir' : '/group/mwavcs/vcs/',
                 'group_account' : {'cpuq':  '#SBATCH --account=pawsey0348',
                                    'gpuq':  '#SBATCH --account=mwavcs',
                                    'copyq': '#SBATCH --account=mwavcs',
                                    'zcpuq': '#SBATCH --account=mwavcs'},
                 'module_dir' : '/group/mwa/software/modulefiles',
                 'presto_module' : 'presto/master',
                 'psrcat_module' : 'psrcat/1.59',
                 'cpuq_cluster' : 'magnus',
                 'cpuq_partition' : 'workq',
                 'gpuq_cluster' : 'galaxy',
                 'gpuq_partition' : 'gpuq',
                 'gpu_beamform_mem' : '1024',
                 'zcpuq_cluster' : 'zeus',
                 'zcpuq_partition' : 'workq',
                 'copyq_cluster' : 'zeus',
                 'copyq_partition' : 'copyq',
                 'container_module' : '',
                 'container_command' : '',
                 'prschive_container' : '/pawsey/mwa/singularity/dspsr/dspsr.sif',
                 'ssd_dir' : None,
                 'gid' : 34858} # mwavcs

GARRAWARLA_CONFIG = {'base_data_dir' : '/astro/mwavcs/vcs/',
                 'base_product_dir' : '/group/mwavcs/vcs/',
                 'group_account' : {'cpuq':  '#SBATCH --account=mwavcs',
                                    'gpuq':  '#SBATCH --account=mwavcs',
                                    'copyq': '#SBATCH --account=mwavcs',
                                    'zcpuq': '#SBATCH --account=mwavcs'},
                 'module_dir' : '/pawsey/mwa/software/python3/modulefiles/',
                 'presto_module' : 'presto/master',
                 'psrcat_module' : 'psrcat/1.59',
                 'cpuq_cluster' : 'garrawarla',
                 'cpuq_partition' : 'workq',
                 'gpuq_cluster' : 'garrawarla',
                 'gpuq_partition' : 'gpuq',
                 'gpu_beamform_mem' : '25600',
                 'zcpuq_cluster' : 'zeus',
                 'zcpuq_partition' : 'workq',
                 'copyq_cluster' : 'zeus',
                 'copyq_partition' : 'copyq',
                 'container_module' : '', #'singularity',
                 'container_command' : '', #'singularity exec --nv /pawsey/mwa/singularity/vcstools_master.sif /bin/bash -c ',
                 'prschive_container' : '/pawsey/mwa/singularity/dspsr/dspsr.sif',
                 'ssd_dir' : '/nvmetmp',
                 'gid' : 34858} # mwavcs

OZSTAR_CONFIG = {'base_data_dir' : '/fred/oz125/vcs/',
                 'base_product_dir' : '/fred/oz125/vcs/',
                 'group_account' : {'cpuq':  '#SBATCH --account=oz125',
                                    'gpuq':  '#SBATCH --account=oz125',
                                    'copyq': '#SBATCH --account=oz125',
                                    'zcpuq': '#SBATCH --account=oz125'},
                 #'module_dir' : '/fred/oz125/software/modulefiles\nmodule use /apps/users/pulsar/skylake/modulefiles',
                 'module_dir' : '/fred/oz125/software/modulefiles',
                 'presto_module' : 'module use /apps/users/pulsar/skylake/modulefiles\nmodule load presto/no-python',
                 #'presto_module' : 'presto/no-python',
                 'psrcat_module' : 'psrcat/1.49',
                 'cpuq_cluster' : 'farnarkle',
                 'cpuq_partition' : 'skylake',
                 'gpuq_cluster' : 'farnarkle',
                 'gpuq_partition' : 'skylake-gpu',
                 'gpu_beamform_mem' : '25600',
                 'copyq_cluster' : 'farnarkle',
                 'copyq_partition' : 'skylake', #TODO check if there's a better one
                 'zcpuq_cluster' : 'farnarkle',
                 'zcpuq_partition' : 'skylake',
                 'container_module' : 'singularity/latest',
                 #removed since I've now installed it on Ozstar
                 #'container_command' : 'singularity exec -H /fred/oz125/vcs/1221832280/ --nv /fred/oz125/container_images/vcstools_multi-pixel.simg'
                 'container_command' : '',
                 'prschive_container' : '',
                 'ssd_dir' : '$JOBFS',
                 'gid' : 10169} # oz125

ARM_CONFIG =   {'base_data_dir' : '/ibo9000/Pulsar/vcs/',
                'base_product_dir' : '/ibo9000/Pulsar/vcs/',
                'group_account' : {'cpuq':  '',
                                   'gpuq':  '',
                                   'copyq': '',
                                   'zcpuq': ''},
                'module_dir' : '/home/app/modulefiles/',
                'presto_module' : 'presto/cpu-master',
                #'psrcat_module' : 'psrcat/1.49',
                'cpuq_cluster' : 'chess',
                'cpuq_partition' : 'arm',
                'gpuq_cluster' : 'chess',
                'gpuq_partition' : 'all-gpu',
                'gpu_beamform_mem' : '30720',
                'copyq_cluster' : 'chess',
                'copyq_partition' : 'arm', #TODO check if there's a better one
                'zcpuq_cluster' : 'chess',
                'zcpuq_partition' : 'arm',
                #None currently as we haven't worked out container software
                'container_module' : '',
                #'container_command' : 'docker run 192.168.6.123:5000/vcstools'}
                'container_command' : '',
                'prschive_container' : '',
                'ssd_dir' : None,
                'gid' : 10002} # pulsar




import logging
import socket
import argparse

logger = logging.getLogger(__name__)

def load_config_file():
    """
    Work out which supercomputer you are using and load the appropriate config file
    """
    #Work out which supercomputer you're using
    hostname = socket.gethostname()
    # galaxy head node, galaxy and magnus job nodes, zeus job nodes, garrawarla job nodes
    if hostname.startswith('galaxy') or hostname.startswith('nid') or hostname.startswith('z'):
        comp_config = GALAXY_CONFIG
    elif hostname.startswith('mwa') or hostname.startswith('garrawarla'):
        comp_config = GARRAWARLA_CONFIG
    elif hostname.startswith('john') or hostname.startswith('farnarkle'):
        comp_config = OZSTAR_CONFIG
    elif hostname.startswith('x86')  or hostname.startswith('arm'):
        comp_config = ARM_CONFIG
    else:
        logger.error('Unknown computer {}. Exiting'.format(hostname))
        quit()

    return comp_config


if __name__ == '__main__':

    # Dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING)

    # Option parsing
    parser = argparse.ArgumentParser("Creates a config file (only required to be run on install or when a new supercomputer is added) and has functions for reading them.")

    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO",
                        choices=loglevels.keys(), default="INFO")

    parser.add_argument("-V", "--version", action='store_true', help="Print version and quit")
    args = parser.parse_args()

    if args.version:
        try:
            import version
            print(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            print("Couldn't import version.py - have you installed vcstools?")
            print("ImportError: {0}".format(ie))
            sys.exit(0)

    # set up the logger for stand-alone execution
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.propagate = False


    #print config file
    config = load_config_file()
    for i in config.keys():
        logger.info("{0}\t{1}".format(i,config[i]))

