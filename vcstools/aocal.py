import numpy as np
import struct, os, logging
from collections import namedtuple

# AOCal variables
HEADER_FORMAT = b"8s6I2d"
HEADER_SIZE = struct.calcsize(HEADER_FORMAT)
HEADER_INTRO = b"MWAOCAL\0"

Header = namedtuple("header", "intro fileType structureType intervalCount antennaCount channelCount polarizationCount timeStart timeEnd")
Header.__new__.__defaults__ = (HEADER_INTRO, 0, 0, 0, 0, 0, 0, 0.0, 0.0)

"""
Collection of database related utilities that are used throughout the VCS processing pipeline
"""
class AOCal(np.ndarray):

    """
    AOCAl stored as a numpy array (with start and stop time stored as floats)

    Array is of dtype complex128 with the following dimensions:

    - calibration interval
    - antenna
    - channel
    - polarisation (order XX, XY, YX, YY)

    The following attributes are made available for convenience, however they
    are not stored explicitly, just read from the array shape.

    aocal.n_int
    aocal.n_ant
    aocal.n_chan
    aocal.n_pol
    """

    def __new__(cls, input_array, time_start=0.0, time_end=0.0):
        """
        See http://docs.scipy.org/doc/numpy-1.10.1/user/basics.subclassing.html
        """
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.time_start = float(time_start)
        obj.time_end = float(time_end)
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.time_start = getattr(obj, 'time_start', None)
        self.time_end = getattr(obj, 'time_end', None)

    def __getattr__(self, name):
        if name == 'n_int':
            return self.shape[0]
        elif name == 'n_ant':
            return self.shape[1]
        elif name == 'n_chan':
            return self.shape[2]
        elif name == 'n_pol':
            return self.shape[3]
        elif name == 'time_start':
            # required to avoid infinite recursion
            return object.__getattribute__(self, time_start)
        elif name == 'time_end':
            # required to avoid infinite recursion
            return object.__getattribute__(self, time_end)
        else:
            raise AttributeError("AOCal has no Attribute %s. Dimensions can be accessed via n_int, n_ant, n_chan, n_pol" % name)

    def strip_edge(self, n_chan):
        """
        return a copy of the array with edge channels removed

        useful for printing without nans but don't write out as calibration solution!
        """
        return self[:, :, n_chan:-n_chan, :]

    def tofile(self, cal_filename):
        if not (np.iscomplexobj(self) and self.itemsize == 16 and len(self.shape) == 4):
            raise TypeError("array must have 4 dimensions and be of type complex128")
        header = Header(intervalCount=self.shape[0], antennaCount = self.shape[1], channelCount = self.shape[2], polarizationCount = self.shape[3], timeStart = self.time_start, timeEnd = self.time_end)
        with open(cal_filename, "wb") as cal_file:
            header_string = struct.pack(HEADER_FORMAT, *header)
            cal_file.write(header_string)
            logging.debug("header written")
            cal_file.seek(HEADER_SIZE, os.SEEK_SET) # skip header. os.SEEK_SET means seek relative to start of file
            np.ndarray.tofile(self, cal_file)
            logging.debug("binary file written")

    def fit(self, pols=(0, 3), mode='model', amp_order=5):
        if not (np.iscomplexobj(self) and self.itemsize == 16 and len(self.shape) == 4):
            raise TypeError("array must have 4 dimensions and be of type complex128")
        for interval in range(self.shape[0]):
            for antenna in range(self.shape[1]):
                logging.debug("fitting antenna %d" % antenna)
                for pol in pols:
                    if sum(~np.isnan(self[interval, antenna, :, pol])) > 0:
                        self[interval, antenna, :, pol] = fit_complex_gains(self[interval, antenna, :, pol])


