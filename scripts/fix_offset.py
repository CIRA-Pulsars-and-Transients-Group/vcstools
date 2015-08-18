#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4


def usage():
    print "fix_offset -o <output file> -f <input fits file> -u <update in place> -s [time per subint (s)] <update subints too>\n"
    print "Read the STT_SMJD and the STT_OFFS in the fitsfile - then set STT_OFFS to 0.0 and fixes the SMJD accordingly\n"

if __name__ == '__main__':

    import sys
    import getopt
    import astropy.io.fits as pyfits
   

    options = {'offset': 0.0, 'file' : "null", 'clobber': 0, 'outfile' : "test.fits", 'do_subints' : 0}
    try:
        opts, args = getopt.getopt(sys.argv[1:],"ho:f:us")
    except getopt.GetoptError:
        usage()
        sys.exit()
    finally:
        if len(sys.argv) < 2:
            usage()
            sys.exit()

    for opt,arg in opts:
        
        if (opt == "-h"):
            usage()
            sys.exit()
        elif (opt == "-o"):
            options['outfile'] = arg
            options['clobber'] = 0
        elif (opt == "-f"):
            options['file'] = arg
        elif (opt == "-u"):
            options['clobber'] = 1
        elif (opt == "-s"):
            options['do_subints'] = 1

    try:
        if (options['clobber'] == 1):
            hdulist = pyfits.open(options['file'],mode='update')
        else:
            hdulist = pyfits.open(options['file'],mode='readonly')
    except:
        print "Error opening the fitsfile"
        sys.exit()

    try:    
        prihdr = hdulist[0].header

        current_offset = float(prihdr['STT_OFFS'])
        current_seconds = int(prihdr['STT_SMJD'])
        current_days = int(prihdr['STT_IMJD'])

        print "before::DAYS = %.d\n" % current_days
        print "before::SECONDS = %.d\n" % current_seconds
        print "before::OFFSET = %.10f\n" % current_offset
        
        if (current_offset >= 0.5):
            prihdr['comment'] = "Changed SMJD from %d to %d" % (current_seconds,current_seconds+1)
            current_seconds = current_seconds+1
            if (current_seconds == 86400):
                current_seconds = 0
                current_days = current_days+1

        prihdr['STT_IMJD'] = current_days
        prihdr['STT_SMJD'] = current_seconds
        prihdr['STT_OFFS'] = 0.0

        prihdr['comment'] = "Changed STT_OFFS to 0.0 to reflect true start precision" 
        
        current_days = int(prihdr['STT_IMJD'])
        current_offset = float(prihdr['STT_OFFS'])
        current_seconds = int(prihdr['STT_SMJD'])
        
        print "after::DAYS = %.d\n" % current_days
        print "after::SECONDS = %.d\n" % current_seconds
        print "after::OFFSET = %.10f\n" % current_offset

        print repr(prihdr)
    except:
        print "Error getting/setting FITSKEYS in PRIMARY header\n"
        sys.exit()       

    if (options['do_subints'] == 1):
        try:
            subint_tbl = hdulist[1].data
            start_lst = float(prihdr['STT_LST'])
            for subnum,entry in enumerate(subint_tbl):
                tmp = round(float(entry[0]))
            
                entry[0]=tmp
                subint_offset = subnum*tmp + 0.5*tmp;
                entry[1]=subint_offset
                entry[2] = start_lst+subint_offset

        except:
            print "Error updating the subints"

    try:    
        if (options['clobber'] == 0):
            print "Will attempt to write %s\n " % options['outfile']
            hdulist.writeto(options['outfile'])
        else:
            print "Will attempt to clobber input file\n"
            hdulist.flush()


        hdulist.close()
        print "Done"
    except:
        print "Error writing out result\n"
        sys.exit()


