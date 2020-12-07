import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

class Stickel:
    def __init__(self, data):
        """Create an instance of a Stickel class.

        Required arguments:
        data -- a (2D) NumPy array with either 2 or 3 columns. If there are only
                2 columns, they are interpreted as (x,y). If 3 columns, they are
                interpreted as (x,y,yerr), where yerr is the absolute size of the
                error (on y).
        """
        # Treat data differently if data has two columns (x, y) vs three columns (x, y, yerr)
        try:
            ncols = data.shape[1]
        except IndexError:
            print("data must be a 2 dimensional array")

        if ncols == 3:
            self.yerr = data[:,2]
        elif ncols == 2:
            self.yerr = None
        else:
            raise ValueError("data must have either 2 or 3 columns")

        self.xdata = data[:,0]
        self.ydata = data[:,1]
        self.yhat  = None
        self.lambda_param = None

        self.xlabel = None
        self.ylabel = None
        self.dylabel = None

        self.N = len(self.xdata)

        self.dyhat = None
        self.dyhatdx = None

        self.diff() # Sets dx, dy, xmid, and dydx

    def plot(self, smooth=False, lambda_param=None, **kwargs):
        plt.plot(self.xdata, self.ydata, '.', label="Orig. data", **kwargs)
        if smooth:
            if lambda_param is None:
                if self.lambda_param is None:
                    raise ValueError("no lambda_param set, required for plotting smoothed y-data")
            else:
                self.smooth_y(lambda_param)
            plt.plot(self.xdata, self.yhat, '.', label="Smoothed data", **kwargs)

        if self.xlabel is not None:
            plt.xlabel(self.xlabel)
        if self.ylabel is not None:
            plt.ylabel(self.ylabel)

    def diff(self):
        self.dx = np.diff(self.xdata)
        self.dy = np.diff(self.ydata)

        self.xmid = (self.xdata[:-1] + self.xdata[1:])/2
        self.dydx = self.dy / self.dx

    def diff_smooth(self, lambda_param=None):
        if lambda_param is None:
            if self.lambda_param is None:
                raise ValueError("no lambda_param set, required for calculating smoothed derivative")
        else:
            self.smooth_y(lambda_param)

        self.dyhat = np.diff(self.yhat)
        self.dyhatdx = self.dyhat / self.dx

    def plot_deriv(self, smooth=False, lambda_param=None, **kwargs):
        if self.xmid is None or self.dydx is None:
            self.diff()

        plt.plot(self.xmid, self.dydx, '.', label="Orig numerical deriv.", **kwargs)

        if smooth:
            self.diff_smooth(lambda_param=lambda_param)
            plt.plot(self.xmid, self.dyhatdx, '.', label="Smoothed deriv.", **kwargs)

        if self.xlabel is not None:
            plt.xlabel(self.xlabel)
        if self.dylabel is not None:
            plt.ylabel(self.dylabel)

    def interp_smooth_deriv(self, **kwargs):
        """kwargs applied to scipy.interpolate.interp1d.
        """
        if self.dyhatdx is None:
            raise ValueError("Must create smoothed deriv (e.g., using diff_smooth()) before calling interp_smooth_deriv()")
        self.dyhatdx_interp = interp1d(self.xmid, self.dyhatdx, kind='cubic', assume_sorted=True, **kwargs)

    def unwrap(self, discont=np.pi):
        self.ydata = np.unwrap(self.ydata, discont)

    def smooth_y(self, lambda_param):
        if lambda_param != self.lambda_param:
            d = 3 # Using two orders higher than the derivative I want (1st), as recommended in Stickel (2010).
            D = self._Dmatrix(d)
            U = self._midpointtruncmatrix(d)
            M = np.identity(self.N)
            delta = np.trace(np.matmul(np.transpose(D), D)) / (self.N**(d+2))

            if self.yerr is not None:
                W = np.diag(1 / self.yerr**2)
            else:
                W = np.identity(self.N)

            MW  = np.matmul(np.transpose(M), W)
            MWM = np.matmul(MW, M)
            MWy = np.matmul(MW, self.ydata)
            DUD = np.matmul(np.matmul(np.transpose(D), U), D)
            inv = np.linalg.pinv(MWM + lambda_param/delta*DUD)

            self.yhat = np.matmul(inv, MWy)
            self.lambda_param = lambda_param

    def _midpointmatrix(self, d):
        B1 = np.roll(self.xdata, 1)
        B2 = np.roll(self.xdata, -1)
        B1[0]  = self.xdata[0]
        B2[-1] = self.xdata[-1]

        return np.diag((-B1+B2)/2)

    def _midpointtruncmatrix(self, d):
        B = self._midpointmatrix(d)
        start = int(np.ceil(d/2))
        finish = self.N - int(np.floor(d/2))

        return B[start:finish,start:finish]

    def _Dmatrix(self, d):
        if d == 0:
            return np.identity(self.N)

        V     = self._Vmatrix(d)
        Dhat1 = Dhat1matrix(self.N-d+1)
        Dprev = self._Dmatrix(d-1)
        VD    = np.matmul(V, Dhat1)
        D     = d*np.matmul(VD, Dprev)

        return D

    def _Vmatrix(self, d):
        return np.diag(1/(self.xdata[d:] - self.xdata[:-d]))

def Dhat1matrix(N):
    ones_p = np.ones(N-1)
    ones_n = -np.ones(N)
    Dhat1 = np.diag(ones_n) + np.diag(ones_p,1)
    Dhat1 = Dhat1[:-1,:]

    return Dhat1

