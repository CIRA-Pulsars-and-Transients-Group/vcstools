#!/usr/bin/env python

import os
import sys

def mdir(path, description, gid=30832):
    # the default groupID is mwaops which is 30832 in numerical
    # we try and make sure all directories created by process_vcs
    # end up belonging to the user and the group mwaops
    # with rwx permissions and the sticky bit set for both user and group
    try:
        os.makedirs(path)
        # we leave the uid unchanged but change gid to mwaops
        os.chown(path,-1,gid)
        os.chmod(path,0771)
        os.system("chmod -R g+s {0}".format(path))
    except:
        if (os.path.exists(path)):
            print "{0} Directory Already Exists\n".format(description)
        else:
            sys.exit(0)


