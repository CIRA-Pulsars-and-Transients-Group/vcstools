import numpy as np
import logging
from scipy.interpolate import UnivariateSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from vcstools.prof_utils import est_sn_from_prof, sigmaClip, check_clip
# Error imports
from vcstools.prof_utils import LittleClipError, LargeClipError, NoComponentsError, ProfileLengthError, NoFitError, BadFitError

logger = logging.getLogger(__name__)


# Main class
class gfit:
    """
    This class is used to fit multiple Gaussians to a pulse profile

    User inputs:
    ------------
    raw_profile: list
        The pulse profile to be fit
    max_N: int
        The maximum number of gaussians to attempt to fit
        Default - 6
    plot_name: str
        The name of the output plot. Can be set with gfit.plot_name. If unsupplied, will use a generic name
        Default - None
    period: float
        The period of the pulsar. Used for SN calculations
        Default - None
    clip_type: str:
        The verbosity of clipping. Choose between 'regular', 'noisy' or 'verbose'.
        Default: 'regular'


    Functions for users:
    --------------------
    auto_gfit:
        Runs the gaussian fit evaluation for a range of values of alpha. This is necessary as there is no way to know
        a priori which alpha to use. Alpha is the input for sigmaClip() and can be interpreted as the level
        of verbosity in the noice-clipping. This function fills out the fit_dict dictionary.

    plot_fit:
        Plots the best chosen gaussian fit to a .png whose name is gfit.plot_name. This can only be run after the fit_dict
        dictionary has been filled (presumably by auto_gfit).


    Products:
    --------
    fit_dict: dictionary
        contains the following keys:
        W10: float
            The W10 width of the profile measured in phase
        W10_e: float
            The uncertainty in the W10
        W50: float
            The W50 width of the profile measured in phase
        W50_e: float
            The uncertainty in the W50
        Weq: float
            The equivalent width of the profile measured in phase
        Weq_e: float
            The uncertainty in the equivalent width
        Wscat: float
            The scattering width of the profile measured in phase
        Wscat_e: float
            The uncertainty in the scattering width
        maxima: list
            A lost of floats corresponding to the bin location of each maximum point
        maxima_e: list
            A list of floats, each correspinding to the error of the maxima of the same index. Measured in bins
        redchisq: float
            The reduced chi sqared of the fit
        num_gauss: int
            The number of gaussian components used in the best fit
        bic: float
            The Bayesian Information Criterion for the best fit
        gaussian_params: list
            A list of length 3*N there N is num_gauss. Each set of 3 parameters corresponds to the amp, centre and width of a guassian component
        cov_mat: np.matrix
            The covariance matrix from the fit
        comp_dict: dictionary
            dict["component_x"] contains an array of the component x
        comp_idx: dictionary
            dict["component_x"] contains an array of indexes of the original profile corresponding to component x
        alpha: float
            The alpha value used in sigmaClip()
        profile: list
            The input profile
        fit: list
            The best fit made into a list form
        sn: float
            The estimated signal to noise ratio, obtained from the profile. Will be None is period unsupplied
        sn_e: float
            The uncertainty in sn. Will be None is period unsupplied
        scattered: boolean
            True is the profile is scattered. Will be None is period unsupplied
    """

    def __init__(self, raw_profile, max_N=6, plot_name=None, min_comp_len=None, period=None, clip_type="regular"):
        self._raw_profile = raw_profile
        self._max_N = max_N
        self._plot_name = plot_name
        self._min_comp_len = min_comp_len
        self._period = period
        self._clip_type = clip_type

        self._fit_dict = {}
        self._best_chi = None
        self._best_alpha = None

        # Private stuff
        self._std_profile = []
        self._alpha = None
        self._noise_std = None
        self._clipped_prof = []

        # Initialize minimum component length if not done by user
        if self._min_comp_len is None:
            self._min_comp_len = int(len(self._raw_profile)/100 + 0.5) + 2
            if self._min_comp_len > 100:
                self._min_comp_len = 100
        if self._min_comp_len < 3:
            self._min_comp_len = 3


    # Setters and Getters
    @property
    def plot_name(self):
        return self._plot_name
    @plot_name.setter
    def plot_name(self, val):
        self._plot_name = val

    @property
    def max_N(self):
        return self._max_N
    @max_N.setter
    def max_N(self, val):
        self._max_N = val

    @property
    def period(self):
        return self._period
    @period.setter
    def period(self, val):
        self._period = val

    @property
    def fit_dict(self):
        return self._fit_dict

    @property
    def best_chi(self):
        return self._best_chi

    @property
    def best_alpha(self):
        return self._best_alpha

    @property
    def std_profile(self):
        return self._std_profile
    @std_profile.setter
    def std_profile(self, val):
        self._std_profile = val

    @property
    def alpha(self):
        return self._alpha
    @alpha.setter
    def alpha(self, val):
        self._alpha = val

    @property
    def noise_std(self):
        return self._noise_std
    @noise_std.setter
    def noise_std(self, val):
        self._noise_std = val

    @property
    def clipped_prof(self):
        return self._clipped_prof
    @clipped_prof.setter
    def clipped_prof(self, val):
        self._clipped_prof = val


    # Main function - intended for use by the user
    def auto_gfit(self):
        """Fits multiple gaussian profiles and finds the best combination of N_Gaussians and alpha"""
        if len(self._raw_profile)<100:
            raise ProfileLengthError("Profile must have length > 100")

        if self._clip_type == "regular":
            alphas = np.linspace(1, 5, 9)
        elif self._clip_type == "noisy":
            alphas = np.linspace(1, 3, 17)
        elif self._clip_type == "verbose":
            alphas = np.linspace(1, 5, 33)
        else:
            raise ValueError("cliptype not recognised. Options are: 'regular', 'noisy' or 'verbose'.")

        # Loop over the gaussian evaluation fucntion, excepting in-built errors
        attempts_dict = {}
        for alpha in alphas:
            try:
                self._alpha = alpha
                prof_dict = self._prof_eval_gfit()
                attempts_dict[alpha] = prof_dict
            except(LittleClipError, LargeClipError, NoComponentsError, ProfileLengthError, BadFitError):
                logger.debug(f"Skipping alpha value: {alpha}")

        # Evaluate the best profile based on reduced chi-squared.
        chi_diff = []
        alphas = []
        for alpha_key in attempts_dict.keys():
            chi_diff.append(abs(1 - attempts_dict[alpha_key]["redchisq"]))
            alphas.append(alpha_key)

        # Sometimes the profile can't be fit
        if not chi_diff:
            raise NoFitError("No suitable profile fit could be found!")

        self._best_chi = min(chi_diff)
        self._best_alpha = alphas[chi_diff.index(self._best_chi)]
        self._fit_dict = attempts_dict[self._best_alpha]

        logger.info("### Best fit results ###")
        logger.info(f"Best model found with BIC of {self._fit_dict['bic']} and reduced Chi of {self._fit_dict['redchisq']} using an alpha value of {self._best_alpha}")
        logger.info(f"W10:                   {self._fit_dict['W10']} +/- {self._fit_dict['W10_e']}")
        logger.info(f"W50:                   {self._fit_dict['W50']} +/- {self._fit_dict['W50_e']}")
        logger.info(f"Wscat:                 {self._fit_dict['Wscat']} +/- {self._fit_dict['Wscat_e']}")
        logger.info(f"Weq:                   {self._fit_dict['Weq']} +/- {self._fit_dict['Weq_e']}")
        logger.info(f"Maxima:                {self._fit_dict['maxima']}")
        logger.info(f"Maxima error:          {self._fit_dict['maxima_e']}")
        if self._fit_dict["sn"]:
            logger.info(f"S/N estimate:          {self._fit_dict['sn']} +/- {self._fit_dict['sn_e']}")


    # Plotter for the resulting fit
    def plot_fit(self):
        """Plots the best fit set of Gaussians"""
        if not self._plot_name:
            self._plot_name = "Gaussian_fit.png"
        y = self._fit_dict["profile"]
        fit = self._fit_dict["fit"]
        popt = self._fit_dict["gaussian_params"]
        maxima = self._fit_dict["maxima"]
        maxima_e= self._fit_dict["maxima_e"]
        x = np.linspace(0, len(y)-1, len(y))
        plt.figure(figsize=(30, 18))

        for j in range(0, len(popt), 3):
            z = self._multi_gauss(x, *popt[j:j+3])
            plt.plot(x, z, "--", label="Gaussian Component {}".format(int((j+3)/3)))
        if maxima:
            for i, mx in enumerate(maxima):
                plt.axvline(x=(mx + maxima_e[i]), ls=":", lw=2, color="gray")
                plt.axvline(x=(mx - maxima_e[i]), ls=":", lw=2, color="gray")

        plt.title(self._plot_name.split("/")[-1].split(".")[0], fontsize=26)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlim(0, len(y))
        plt.xlabel("Bins", fontsize=24)
        plt.ylabel("Intensity", fontsize=24)
        plt.plot(x, y, label="Original Profile", color="black")
        plt.plot(x, fit, label="Gaussian Model", color="red")
        plt.legend(loc="best", prop={'size': 18})
        plt.savefig(self._plot_name)
        plt.close()


    def _standardise_raw_profile(self):
        """Clips and normalises the raw profile"""
        y = np.array(self._raw_profile)/max(self._raw_profile)
        noise_std, self._clipped_prof = sigmaClip(y, alpha=self._alpha)
        check_clip(self._clipped_prof)
        ignore_threshold = 3 * noise_std
        y = y - np.nanmean(self._clipped_prof)
        y = y/max(y)
        self._noise_std = np.nanstd(np.array(self._clipped_prof)/max(y))
        self._std_profile = y


    def _prof_eval_gfit(self):
        """Fits multiple gaussians to a profile and subsequently finds W10, W50, Weq and maxima"""
        # Normalize, find the std
        self._standardise_raw_profile()

        # Fit gaussian
        fit, chisq, bic, popt, pcov, comp_dict, comp_idx = self._fit_gaussian()
        fit = np.array(fit)
        n_rows, _ = np.shape(pcov)
        num_gauss = n_rows/3

        # Find widths + error
        W10, W50, Weq, Wscat, W10_e, W50_e, Weq_e, Wscat_e = self._find_widths(popt, pcov)

        # Convert from bins to phase
        proflen = len(self._std_profile)
        W10 = W10/proflen
        W50 = W50/proflen
        Weq = Weq/proflen
        Wscat = Wscat/proflen
        W10_e = W10_e/proflen
        W50_e = W50_e/proflen
        Weq_e = Weq_e/proflen
        Wscat_e = Wscat_e/proflen

        _, maxima, _, maxima_e = self._find_minima_maxima_gauss(popt, pcov, len(fit))


        # Estimate SN
        if self._period:
            sn, sn_e, scattered = est_sn_from_prof(self._std_profile, self._period, alpha=self._alpha)
        else:
            sn = None
            sn_e = None
            scattered = None

        # Dump to dictionary
        fit_dict = {"W10":W10, "W10_e":W10_e, "W50":W50, "W50_e":W50_e, "Wscat":Wscat, "Wscat_e":Wscat_e,
                    "Weq":Weq, "Weq_e":Weq_e, "maxima":maxima, "maxima_e":maxima_e, "redchisq":chisq,
                    "num_gauss":num_gauss, "bic":bic, "gaussian_params":popt, "cov_mat":pcov, "comp_dict":comp_dict,
                    "comp_idx":comp_idx, "alpha":self._alpha, "profile":self._std_profile, "fit":fit, "sn":sn,
                    "sn_e":sn_e, "scattered":scattered}

        return fit_dict


    # The function that does the actual fitting
    def _fit_gaussian(self):
        """
        Fits multiple gaussian components to a pulse profile and finds the best number to use for a fit.
        Will always fit at least one gaussian per profile component.
        Profile components are defined by find_components().
        Each gaussian is defined by the following: y = amp * np.exp( -((x - ctr)/wid)**2)

        Returns:
        --------
        [fit, redchisq, best_bic, popt, pcov]: list
            fit: list
                The data containing the multi-component gaussian fit to the input profile
            redchisq: float
                The reduced chi-sqaured value of the fit
            best_bic: float
                The bayesian information criterion for the fit
            popt: list
                A list of floats where each 3 numbers describes a single gaussain and are 'ctr', 'amp' and 'wid' respectively
            pcov: numpy matrix
                The covariance matrix generated by the curve_fit function
        """
        # Chi sqaured evaluation
        def chsq(observed_values, expected_values, err):
            test_statistic=0
            for observed, expected in zip(observed_values, expected_values):
                test_statistic+=((float(observed)-float(expected))/float(err))**2
            return test_statistic

        # Find profile components
        self._fill_clipped_prof(search_scope=int(len(self._std_profile)/100))
        comp_dict, comp_idx = self._find_components()

        # Estimate gaussian parameters based on profile components
        comp_centres = []
        comp_max = []
        comp_width = []
        for i in range(self._max_N//len(comp_idx.keys())+1):
            for key in comp_idx.keys():
                comp_centres.append(np.mean(comp_idx[key]))
                comp_max.append(max(comp_dict[key])*0.5)
                comp_width.append((max(comp_idx[key])-min(comp_idx[key])))
        centre_guess = iter(comp_centres)
        width_guess=iter(comp_width)
        max_guess=iter(comp_max)

        n_comps=len(comp_dict.keys())
        logger.debug(f"Number of profile components: {n_comps} ({comp_centres[:n_comps]})")

        # Fit from 1 to max_N gaussians to the profile. Evaluate profile fit using bayesian information criterion
        x = np.linspace(0, len(self._std_profile)-1, len(self._std_profile))
        bounds_arr=[[],[]]
        guess = []
        fit_dict = {}

        for num in range(1, self._max_N):
            # Calculate bounds and guess
            guess += [next(max_guess), next(centre_guess), next(width_guess)]
            bounds_arr[0].append(0)
            bounds_arr[0].append(0)
            bounds_arr[0].append(0)
            bounds_arr[1].append(max(self.std_profile))
            bounds_arr[1].append(len(self.std_profile))
            bounds_arr[1].append(len(self.std_profile))
            bounds_tuple=(tuple(bounds_arr[0]), tuple(bounds_arr[1]))

            # Do the fitting
            maxfev_list = [10000, 50000, 100000, 500000]
            popt = None
            pcov = None
            for maxfev in maxfev_list:
                try:
                    logger.debug(f"Maxfev in curve_fit: {i}")
                    popt, pcov = curve_fit(self._multi_gauss, x, self._std_profile, bounds=bounds_tuple, p0=guess, maxfev=maxfev)
                    break # Break if fit successful
                except (RuntimeError, ValueError) as e: # No fit found in maxfev
                    pass
            if popt is None or pcov is None:
                raise BadFitError("Could not successfully execute least squares fit")

            # Work out some stuff
            fit = self._multi_gauss(x, *popt)
            chisq = chsq(self._std_profile, fit, self._noise_std)

            # Bayesian information criterion for gaussian noise
            k = 3*(num+1)
            bic = chisq + k*np.log(len(self.std_profile))
            fit_dict[str(num+1)]={"popt":[], "pcov":[], "fit":[], "chisq":[], "bic":[]}
            fit_dict[str(num+1)]["popt"] = popt
            fit_dict[str(num+1)]["pcov"] = pcov
            fit_dict[str(num+1)]["fit"] = fit
            fit_dict[str(num+1)]["redchisq"] = chisq/(len(self.std_profile)-1)
            fit_dict[str(num+1)]["bic"] = bic
            logger.debug(f"Reduced chi squared for               {num+1} components: {fit_dict[str(num+1)]['redchisq']}")
            logger.debug(f"Bayesian Information Criterion for    {num+1} components: {fit_dict[str(num+1)]['bic']}")

        # Find the best fit according to the BIC
        best_bic = np.inf
        best_fit = None
        for n_components in fit_dict.keys():
            if fit_dict[n_components]["bic"] < best_bic:
                best_bic = fit_dict[n_components]["bic"]
                best_fit = n_components
        logger.debug(f"Fit {best_fit} gaussians for a reduced chi sqaured of {fit_dict[best_fit]['redchisq']}")
        popt = fit_dict[best_fit]["popt"]
        pcov = fit_dict[best_fit]["pcov"]
        fit = fit_dict[best_fit]["fit"]
        redchisq = fit_dict[best_fit]["redchisq"]

        return [fit, redchisq, best_bic, popt, pcov, comp_dict, comp_idx]

    # SigmaClip isn't perfect. Use these next function to fix bad clips
    def _fill_clipped_prof(self, search_scope=None, nan_type=0.):
        """
        Intended for use on noisy profiles. Fills nan values that are surrounded by non-nans to avoid discontinuities
        in the profile
        """
        length = len(self._clipped_prof)
        if search_scope is None:
            # Search 5% ahead for non-nans
            search_scope = round(length*0.05)
        search_scope = np.linspace(1, search_scope, search_scope, dtype=int)

        # Loop over all values in clipped profile
        for i, val in enumerate(self._clipped_prof):
            if val == nan_type and not (i+max(search_scope)) >= length:
                #look 'search_scope' indices ahead for non-nans
                for j in sorted(search_scope, reverse=True):
                    #fill in nans
                    if self._clipped_prof[i+j]==nan_type:
                        for k in range(j):
                            self._clipped_prof[i+k]=nan_type
                        break


    # Find the individual profile components
    def _find_components(self):
        """
        Given a profile in which the noise is clipped to 0, finds the components that are clumped together.

        Returns:
        --------
        component_dict: dictionary
            dict["component_x"] contains an array of the component x
        component_idx: dictionary
            dict["component_x"] contains an array of indexes of the original profile corresponding to component x
        """
        # Find the components by looking at the clipped profile
        component_dict={}
        component_idx={}
        num_components=0
        for i, val in enumerate(self._clipped_prof):
            if np.isnan(val): # On pulse
                if not np.isnan(self._clipped_prof[i-1]) or i==0: # Off pulse or beginning of profile phase
                    num_components+=1
                    comp_key = f"component_{num_components}"
                    component_dict[comp_key]=[]
                    component_idx[comp_key]=[]
                component_dict[comp_key].append(self._std_profile[i])
                component_idx[comp_key].append(i)

        del_comps = []
        for comp_key in component_dict.keys():
            if len(component_dict[comp_key]) < self._min_comp_len or max(component_dict[comp_key]) < 0.:
                del_comps.append(comp_key)
        for i in del_comps:
            del component_dict[i]
            del component_idx[i]

        if len(component_dict.keys()) == 0:
            raise NoComponentsError("No profile components have been found")

        return component_dict, component_idx


    # Find min and max points
    def _find_minima_maxima_gauss(self, popt, pcov, x_length):
        """
        Finds all roots of a gaussian function

        Parameters:
        -----------
        popt: list
            A list of length 3N where N is the number of gaussians. This list contains the parameters amp, mean, centre respectively
        pcov: np.matrix
            The covariance matric corresponding to the parameters from popt
        x_length: int
            The length of the list used to fit the gaussian

        Returns:
        --------
        minima: list
            A list of the minimum points of the fit
        maxima: list
            A list of the maximum points of the fit
        minima_e: list
            The error in each minima
        maxima_e: list
            The error in each maxima
        """
        # Create the derivative list and spline it to find roots
        x = np.linspace(0, x_length-1, x_length)
        dy = self._multi_gauss_ddx(x, *popt)
        spline_dy = UnivariateSpline(x, dy, s=0)
        roots = spline_dy.roots()

        # Find which are max and min
        maxima = []
        minima = []
        for root in roots:
            idx = int(root + 0.5)
            if dy[idx-1] > dy[idx]:
                maxima.append(root)
            else:
                minima.append(root)

        minima_e = self._find_x_err(minima, popt, pcov)
        maxima_e = self._find_x_err(maxima, popt, pcov)

        return minima, maxima, minima_e, maxima_e


    # Find the width of profile components
    def _find_widths(self, popt, pcov):
        """
        Attempts to find the W_10, W_50 and equivalent width of a profile by using a spline approach.
        W10 and W50 errors are estimated by using: sigma_x = sigma_y/(dy/dx)
        Weq errors are estimated by finding the average difference in Weq when you add and subtract the std from the on-pulse profile

        Parameters:
        -----------
        popt: list
            The parameters that are used to create the multi-gaussian fit
        pcov: np.matrix
            The covariance matrix corresponding to the parameters from popt

        Returns:
        --------
        [W10, W50, Weq, Wscat, W10_e, W50_e, Weq_e, Wscat_e]: list
            W10: float
                The W10 width of the profile measured in number of bins
            W50: float
                The W50 width of the profile measured in number of bins
            Weq: float
                The equivalent width of the profile measured in number of bins
            Wscat: float
                The scattering width of the profile measured in number of bins
            W10_e: float
                The uncertainty in W10
            W50_e: float
                The uncertainty in W50
            Weq_e: float
                The uncertainty in Weq
            Wscar_e: float
                The unceratinty in Wscat
        """
        def error_in_x_pos(pcov, popt, x):
            J = self._jacobian_slope(x, *popt)
            JC = np.matmul(J, pcov)
            sigma_y = np.sqrt(np.matmul(JC, np.transpose(J)).item(0))
            ddx = self._multi_gauss_ddx(x, *popt)
            return sigma_y/ddx

        # Perform spline operations on the fit
        x = np.linspace(0, len(self._std_profile)-1, len(self._std_profile))
        fit = self._multi_gauss(x, *popt)
        amp_fit = max(fit) - min(fit)
        spline10 = UnivariateSpline(x, fit - np.full(len(x), 0.1*amp_fit), s=0)
        spline50 = UnivariateSpline(x, fit - np.full(len(x), 0.5*amp_fit), s=0)
        spline_s = UnivariateSpline(x, fit - np.full(len(x), 1/np.exp(1)*amp_fit), s=0)

        # Find Weq using the real profile
        std, off_pulse = sigmaClip(self._std_profile, alpha=self._alpha)
        check_clip(off_pulse)
        on_pulse=[]
        for i, data in enumerate(off_pulse):
            if np.isnan(data):
                on_pulse.append(self._std_profile[i])
        x = np.linspace(0, len(on_pulse)-1, len(on_pulse))
        spline0 = UnivariateSpline(x, on_pulse, s=0)
        integral = spline0.integral(0, len(on_pulse)-1)
        Weq = integral/max(on_pulse)

        # Find W10, W50 and Wscat
        W10_roots = spline10.roots()
        W50_roots = spline50.roots()
        Wscat_roots = spline_s.roots()
        # If any of the widths roots don't exist, this fit sucks
        if not list(W10_roots) or not list(W50_roots) or not list(Wscat_roots):
            raise BadFitError("This fit is not suitable for use")
        W10 = W10_roots[-1] - W10_roots[0]
        W50 = W50_roots[-1] - W50_roots[0]
        Wscat = Wscat_roots[-1] - Wscat_roots[0]

        # W10 root errors
        err_10_1 = error_in_x_pos(pcov, popt, W10_roots[0])
        err_10_2 = error_in_x_pos(pcov, popt, W10_roots[-1])
        W10_e = np.sqrt(err_10_1**2 + err_10_2**2)

        # W50 root errors
        err_50_1 = error_in_x_pos(pcov, popt, W50_roots[0])
        err_50_2 = error_in_x_pos(pcov, popt, W50_roots[-1])
        W50_e = np.sqrt(err_50_1**2 + err_50_2**2)

        # Wscat root errors
        err_scat_1 = error_in_x_pos(pcov, popt, Wscat_roots[0])
        err_scat_2 = error_in_x_pos(pcov, popt, Wscat_roots[-1])
        Wscat_e = np.sqrt(err_scat_1**2 + err_scat_2**2)

        # Weq errors - using covariance formula
        on_pulse_less = (on_pulse - std).clip(min=0)
        spline0 = UnivariateSpline(x, on_pulse_less, s=0)
        integral = spline0.integral(0, len(self._std_profile)-1)
        dwdint = 1/max(on_pulse)**2
        dwdmax = -integral/max(on_pulse)**2
        int_e = abs(integral/max(on_pulse - std) - integral/max(on_pulse))
        max_e = std
        Weq_e = np.sqrt( dwdint**2 * int_e**2 + dwdmax**2 * max_e**2 + 2*dwdint*dwdmax*int_e*max_e )

        return [W10, W50, Weq, Wscat, W10_e, W50_e, Weq_e, Wscat_e]


    def _find_x_err(self, x, popt, pcov):
        """
        Finds the error in the horizontal position of a gaussian fit at the point x.
        Uses the equation sigma_x = sigma_y/d2ydx2 where:
        sigma_x = error in x
        d2ydx2 = second derivative of the gaussian function at point x
        sigma_y = J*C*J_T
        J = Jacobian evalutated at point x
        C = covariance matrix of gaussian fit
        J_T = transposed jacobian

        Parameters:
        -----------
        x: list
            A list of points to evaluate the error at
        popt: list
            The parameters used to describe the gaussian fit
        pcov: numpy.matrix
            The covariance matrix corresponding to popt

        Returns:
        --------
        x_err: list
            The error evaluated at each point, x
        """
        x_err = []
        for _, point in enumerate(x):
            J = self._jacobian_slope(point, *popt)
            d2dx2 = self._multi_gauss_d2dx2(point, *popt)
            JC = np.matmul(J, pcov)
            sigma_y = np.sqrt( np.matmul(JC, np.transpose(J)).item(0) )
            x_err.append(sigma_y / abs(d2dx2))
        return x_err


    # A bunch of maths stuff
    def _integral_multi_gauss(self, *params):
        y=0
        for i in range(0, len(params), 3):
            a = params[i]
            c = params[i+2]
            y = y + a*c*np.sqrt(2*np.pi)
        return y


    def _multi_gauss(self, x, *params):
        y = np.zeros_like(x)
        for i in range(0, len(params), 3):
            a = params[i]
            b = params[i+1]
            c = params[i+2]
            y = y +  a * np.exp( -(((x-b)**2) / (2*c**2)) )
        return y


    def _multi_gauss_ddx(self, x, *params):
        # Derivative of gaussian
        y = np.zeros_like(x)
        for i in range(0, len(params), 3):
            a = params[i]
            b = params[i+1]
            c = params[i+2]
            y = y - a/c**2 * (x - b) * np.exp( -(((x-b)**2) / (2*c**2)) )
        return y


    def _multi_gauss_d2dx2(self, x, *params):
        # Double derivative of gaussian
        y = np.zeros_like(x)
        for i in range(0, len(params), 3):
            a = params[i]
            b = params[i+1]
            c = params[i+2]
            y = y + (self._multi_gauss(x, a, b, c) / c**2) * (((x - b)**2)/(c**2) - 1)
        return y


    def _partial_gauss_dda(x, a, b, c):
            return np.exp((-(b - x)**2)/(2*c**2))
    def _partial_gauss_ddb(x, a, b, c):
            return a*(x - b) * np.exp((-(b - x)**2)/(2*c**2))/c**2
    def _partial_gauss_ddc(x, a, b, c):
            return a*(x - b)**2 * np.exp((-(b - x)**2)/(2*c**2))/c**3


    def _jacobian_slope(self, x, *params):
        """
        Evaluates the Jacobian matrix of a gaussian slope at a single point, x

        Parameters:
        -----------
        x: float
            The point to evaluate
        *params: list
            A list containing three parameters per gaussian component in the order: Amp, Mean, Width

        Returns:
        --------
        J: numpy.matrix
            The Jacobian matrix
        """
        def dda(a, b, c, x):
            return -self._multi_gauss(x, a, b, c) * (x - b)/(c**2)/a
        def ddb(a, b, c, x):
            return self._multi_gauss(x, a, b, c) * (1 - (x - b)**2/(c**2))/c**2
        def ddc(a, b, c, x):
            return self._multi_gauss(x, a, b, c) * (x - b)/(c**3) * (2 - (x-b)**2/(c**2))
        J = []
        for i in range(0, len(params), 3):
            a = params[i]
            b = params[i+1]
            c = params[i+2]
            mypars = [a, b, c, x]
            J.append(dda(*mypars))
            J.append(ddb(*mypars))
            J.append(ddc(*mypars))
        J = np.asmatrix(J)
        return J