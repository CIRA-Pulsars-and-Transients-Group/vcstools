#!/usr/bin/env python
import sys, os, time, socket, json, struct
from optparse import OptionParser
import threading
import urllib2, urllib
import base64
import time

username = 'ngas'
password = base64.decodestring('bmdhcw==')
    
class PrintStatus():
    
   def __init__(self, numfiles):
      self.status = {}
      self.lock = threading.RLock()
      self.currentbytes = 0
      self.totalbytes = 0
      self.runtime = 0
      self.errors = []
      self.files = 0
      self.filesComplete = 0
      self.totalfiles = numfiles;
   
   
   def fileError(self, err):
       with self.lock:
           self.errors.append(err)
           print err
       
   def fileStarting(self, filename):
       with self.lock:
           print "%s [INFO] Downloading %s" % (time.strftime("%c"), filename)
           
   def fileComplete(self, filename):
       with self.lock:
           self.filesComplete = self.filesComplete + 1
           print "%s [INFO] %s complete [%d of %d]" % (time.strftime("%c"), filename, self.filesComplete, self.totalfiles)


def splitRawRecombined(filename):
   
   try:
      file = os.path.basename(filename)
      if '.dat' not in file:
         raise Exception('dat extension not found')
   
      part = file.split('_')
      if 'ch' not in part[2]:
         raise Exception('ch not found in 3rd part')
      
      obsid = int(part[0])
      
      try:
         # 1070978272_1401856338_ch164.dat variant
         tm = int(part[1])
         chan = part[2].split('.')[0]
         
      except:
         #1070978272_c_ch05_1386943943.dat variant
         chan = part[2]
         tm = int(part[3].split('.')[0])
      
      return obsid, tm, chan
            
   except Exception as e:
      raise Exception('invalid voltage recombined product filename %s' % file)


def splitRawVoltage(filename):
   try:
      file = os.path.basename(filename)
      if '.dat' not in file:
         raise Exception('dat extension not found')
      
      part = file.split('_')
      if 'vcs' not in part[2]:
         raise Exception('vcs not found in 3rd part')
      
      return (int(part[0]), int(part[1]), part[2], int(part[3].split('.')[0]))
   
   except Exception as e:
      raise Exception('invalid voltage data filename %s' % file)


def queryObs(obs, host, type, timefrom, duration):
   
   processRange = False
   
   if timefrom != None and duration != None:
      processRange = True
   
   url = 'http://%s/metadata/obs/?obs_id=%s&filetype=%s' % (host, str(obs), str(type))
   
   response = None
   
   try:
      request = urllib2.Request(url)
      response = urllib2.urlopen(request)
      
      resultbuffer = []
      while True:
        result = response.read(32768)
        if not result:
          break

        resultbuffer.append(result)

      keymap = {}
      files = json.loads(''.join(resultbuffer))['files']
      if processRange:
         time = None
         for f, v in files.iteritems():
            if type == 11:
               obsid, time, vcs, part = splitRawVoltage(f)
            elif type == 12:
               obsid, time, chan = splitRawRecombined(f)
               
            if time >= timefrom and time <= timefrom+duration:
               
               keymap[f] = v['size']
      else:
         for f, v in files.iteritems():
            keymap[f] = v['size']
      
      return keymap
   
   finally:
      if response:
         response.close()

 

def checkFile(filename, size, dir):
    path = dir + filename
        
    # check the file exists
    if os.path.isfile(path) is True:
        #check the filesize matches
        filesize = os.stat(path).st_size

        if filesize == int(size):
            return True
        
    return False


def worker(url, size, filename, s, out, stat, bufsize, prestage):
    u = None
    f = None
    
    try:
        stat.fileStarting(filename)

        # open file URL
        request = urllib2.Request(url)
        base64string = base64.encodestring('%s:%s' % (username, password)).replace('\n', '')
        request.add_header("Authorization", "Basic %s" % base64string)   
        request.add_header('prestagefilelist', prestage)
        
        u = urllib2.urlopen(request)
        u.fp.bufsize = bufsize
        
        # get file size
        meta = u.info()
        file_size = int(meta.getheaders("Content-Length")[0])
        
        # open file for writing
        f = open(out + filename, 'wb')
        
        file_size_dl = 0
        while True:
            buff = u.read(bufsize)
            if not buff:
              break

            f.write(buff)
            file_size_dl += len(buff)

        if file_size_dl != file_size:
          raise Exception("size mismatch %s %s" % str(file_size), str(file_size_dl))

        stat.fileComplete(filename)
        
    except urllib2.HTTPError as e:
        stat.fileError("%s [ERROR] %s %s" % (time.strftime("%c"), filename, str(e.read()) ))
    
    except urllib2.URLError as urlerror:
        if hasattr(urlerror, 'reason'):
            stat.fileError("%s [ERROR] %s %s" % (time.strftime("%c"), filename, str(urlerror.reason) ))
        else:
            stat.fileError("%s [ERROR] %s %s" % (time.strftime("%c"), filename, str(urlerror) ))
    
    except Exception as exp:
        stat.fileError("%s [ERROR] %s %s" % (time.strftime("%c"), filename, str(exp) ))
        
    finally:
        if u:
            u.close()
            
        if f:
            f.flush()
            f.close()
            
        s.release()


def main():
    stat = None
    
    try:
        parser = OptionParser(usage="usage: %prog [options]", version="%prog 1.0")
        parser.add_option("--obs", action="store", dest="obs", help="Observation ID")
        parser.add_option("--type", action='store', type = 'int', dest='filetype', help='Voltage data type (Raw = 11, Recombined Raw = 12)')
        parser.add_option("--from", action='store', type = 'int', dest='timefrom', help='Time from (taken from filename)')
        parser.add_option("--duration", default=0, type = 'int', dest='duration', help='Duration')
        parser.add_option("--ngas",  default='fe4.pawsey.ivec.org:7790', action="store", dest="ngashost", help="NGAS server (default: fe4.pawsey.ivec.org:7790)")
        parser.add_option("--dir", default= './', action="store", dest="out", help="Output directory (default: ./<Observation ID>")
        parser.add_option("--parallel", default='6', action="store", dest="td", help="Number of simultaneous downloads (default: 6)")
        
        bufsize = 4096
        
        (options, args) = parser.parse_args()
        
        if options.ngashost == None:
            print 'NGAS host not defined'
            sys.exit(-1)
            
        if options.obs == None:
            print 'Observation ID is empty'
            sys.exit(-1)
            
        if options.filetype == None:
            print 'File type not specified'
            sys.exit(-1)
        
        if options.timefrom == None:
            print 'Time from not specified'
            sys.exit(-1)
        
        if options.timefrom != None and options.duration != None:
            if options.duration < 0:
               print 'Duration must not be negative'
               sys.exit(-1)
            
        numdownload = int(options.td)
        
        if numdownload <= 0 or numdownload > 12:
            print 'Number of simultaneous downloads must be > 0 and <= 12'
            sys.exit(-1)
        
        print '%s [INFO] Finding observation %s' % (time.strftime("%c"), options.obs)

        fileresult = queryObs(options.obs, 'mwa-metadata01.pawsey.org.au', options.filetype, options.timefrom, options.duration)
        if len(fileresult) <= 0:
            print '%s [INFO] No files found for observation %s and file type %s' % (time.strftime("%c"), options.obs, int(options.filetype))
            sys.exit(1)
        
        print '%s [INFO] Found %s files' % (time.strftime("%c"), str(len(fileresult)))
            
        if len(fileresult) > 12000:
            print '%s [INFO] File limit exceeded 12000, please stagger your download' % (time.strftime("%c"))
            sys.exit(1)
     
        # advise that we want to prestage all the files
        filelist = []
        for key, value in fileresult.iteritems():
           filelist.append(key)
         
        prestageStr = json.dumps(filelist)
     
        if options.out == None or len(options.out) == 0:
            options.out = './' + options.out + '/'

        # check we have a forward slash before file
        if options.out[len(options.out)-1] != '/':
             options.out += '/'

        dir = options.out # + options.obs + '/'
        if not os.path.exists(dir):
            os.makedirs(dir)
        
        stat = PrintStatus(len(fileresult))
        urls = []
        
        for key, value in sorted(fileresult.iteritems()):
            url = "http://" + options.ngashost + "/RETRIEVE?file_id=" + key
            
            if checkFile(key, int(value), dir) is False:
                urls.append((url, value, key))  
            else:
                stat.fileComplete(key)
             
        s = threading.BoundedSemaphore(value=numdownload)        
        for u in urls:
            while True:
                if s.acquire(blocking=False):
                    t = threading.Thread(target=worker, args=(u[0],u[1],u[2],s,dir,stat,int(bufsize), prestageStr))
                    t.setDaemon(True)
                    t.start()
                    break
                else:
                    time.sleep(1)

        while True:
            main_thread = threading.currentThread()
            #if the main thread and print thread are left then we must be done; else wait join
            if len(threading.enumerate()) == 1:
                break;
            
            for t in threading.enumerate():
                #if t is main_thread or t is status_thread:
                if t is main_thread:
                    continue
                
                t.join(1)
                if t.isAlive():
                    continue

        print "%s [INFO] File Transfer Complete." % (time.strftime("%c"))
        
        # check if we have errors
        if len(stat.errors) > 0:
            print "%s [INFO] File Transfer Error Summary:" % (time.strftime("%c"))
            for i in stat.errors:
                print i
                
            raise Exception()
        else:
            print "%s [INFO] File Transfer Success." % (time.strftime("%c"))
            
    except KeyboardInterrupt as k:
        raise k
    
    except Exception as e:
        raise e

if __name__ == "__main__":
    try:
        main()
        sys.exit(0)
        
    except KeyboardInterrupt as k:
        print 'Interrupted, shutting down'
        sys.exit(2)
    
    except Exception as e:
        print e
        sys.exit(3)
