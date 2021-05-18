from vcstools.prof_utils import ProfileLengthError, BadFitError
from vcstools.prof_utils import estimate_components_onpulse, error_in_x_pos
import itertools
from scipy.stats import vonmises
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import logging
from scipy.interpolate import UnivariateSpline
import matplotlib
matplotlib.use('Agg')

# Error imports

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
    conponent_plot_name: str
        The name of the plot for the component plot that illustrates how profile components are chosen. If 
        unsupplied, will not plot.
        Default - None


    Functions for users:
    --------------------
    auto_fit:
        Runs a gaussian fit evaluation using up to max_N gaussians. The quality of the fit is decided via a chi-square evaluation.
        The best of these is saved and used to determine widths and maxima location. Part of the preprocessing of the profile will
        also determine the on-pulse regions and subsequently the noise level and signal to noise ratio.
    plot_fit:
        Plots the best chosen gaussian fit to a file whose name is gfit.plot_name. This can only be run after the fit_dict
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
            A lost of floats corresponding to the phase location of each maximum point
        maxima_e: list
            A list of floats, each correspinding to the error of the maxima of the same index. Measured in phase
        minima: list
            A lost of floats corresponding to the phase location of each minimum point
        minima_e: list
            A list of floats, each correspinding to the error of the minima of the same index. Measured in phase
        redchisq: float
            The reduced chi sqared of the fit
        n_components: int
            The number of gaussian components used in the best fit
        n_profile_components: int
            The number of distinguishable profile components
        on_pulse_estimates: list
            The output of prof_utils.estimate_components_onpulse()
        fit_params: list
            A list of length 5 + 3*(N-1) there N is num_gauss. Each set of 3 parameters corresponds to the amp, centre and width of a guassian component
            The first five are baseline, tau, amp, centre and width
        cov_mat: np.matrix
            The covariance matrix output from curve_fit
        profile: list
            The normalised and rolled input profile
        fit: list
            The best fit made into a list form
        sn: float
            The estimated signal to noise ratio, obtained from the profile
        sn_e: float
            The uncertainty in sn
        scattering: float
            The tau value of the profile fit
    """

    def __init__(self, raw_profile, max_N=6, plot_name=None, component_plot_name=None):
        # Initialise inputs
        self._raw_profile = raw_profile
        self._max_N = max_N
        self._plot_name = plot_name
        self._component_plot_name = component_plot_name

        # The return dictionary
        self._fit_dict = {}

        # Stuff to put in the dictionary
        self._fit_profile = None
        self._best_chi = None
        self._std_profile = None
        self._noise_std = None
        self._n_off_pulse = None
        self._n_gaussians = None
        self._comp_est_dict = None
        self._n_prof_components = None
        self._popt = None
        self._pcov = None
        self._W10 = None
        self._W50 = None
        self._Weq = None
        self._Wscat = None
        self._W10_e = None
        self._W50_e = None
        self._Weq_e = None
        self._Wscat_e = None
        self._minima = None
        self._maxima = None
        self._minima_e = None
        self._maxima_e = None

    # Setters and Getters

    @property
    def plot_name(self):
        return self._plot_name

    @plot_name.setter
    def plot_name(self, val):
        self._plot_name = val

    @property
    def component_plot_name(self):
        return self.component_plot_name

    @component_plot_name.setter
    def component_plot_name(self, val):
        self.component_plot_name = val

    @property
    def max_N(self):
        return self._max_N

    @max_N.setter
    def max_N(self, val):
        self._max_N = val

    @property
    def fit_dict(self):
        return self._fit_dict

    # Main function - intended for use by the user

    def auto_fit(self):
        """Fits multiple gaussian profiles and finds the best combination of N_Gaussians and alpha"""
        if len(self._raw_profile) < 100:
            raise ProfileLengthError("Profile must have length > 100")

        # Standarsdise the profile so that it's easy to fit
        self._standardise_raw_profile()
        # Estimate the on-pulse region, components and noise level
        self._comp_est_dict = estimate_components_onpulse(
            self._std_profile, plot_name=self._component_plot_name)
        self._n_prof_components = len(self._comp_est_dict["est_on_pulse"])
        self._noise_std = self._comp_est_dict["noise"]

        # Loop over the gaussian fitting function, excepting in-built errors
        try:
            self._fit_gaussian()
        except BadFitError as e:
            logger.error("No Fit Found!")

        # Do all the fancy things with the fit
        self._find_widths() # Find widths + error
        self._find_minima_maxima_gauss() # Find turning points
        self._est_sn() # Estimate SN

        # Dump to dictionary
        self._fit_dict = {}
        self._fit_dict["W10"] = self._W10
        self._fit_dict["W10_e"] = self._W10_e
        self._fit_dict["W50"] = self._W50
        self._fit_dict["W50_e"] = self._W50_e
        self._fit_dict["Wscat"] = self._Wscat
        self._fit_dict["Wscat_e"] = self._Wscat_e
        self._fit_dict["Weq"] = self._Weq
        self._fit_dict["Weq_e"] = self._Weq_e
        self._fit_dict["maxima"] = self._maxima
        self._fit_dict["maxima_e"] = self._maxima_e
        self._fit_dict["maxima"] = self._maxima
        self._fit_dict["maxima_e"] = self._maxima_e
        self._fit_dict["redchisq"] = self._best_chi
        self._fit_dict["n_components"] = self._n_gaussians
        self._fit_dict["fit_params"] = self._popt
        self._fit_dict["cov_mat"] = self._pcov
        self._fit_dict["profile"] = self._std_profile
        self._fit_dict["fit"] = self._fit_profile
        self._fit_dict["sn"] = self._sn
        self._fit_dict["sn_e"] = self._sn_e
        self._fit_dict["noise"] = self._noise_std
        self._fit_dict["scattering"] = self._popt[1] # tau
        self._fit_dict["on_pulse_estimates"] = self._comp_est_dict
        self._fit_dict["n_profile_components"] = self._n_prof_components

        # Log all the important stuff
        logger.info("### Best fit results ###")
        logger.info(
            f"Best model found with reduced Chi square of {self._fit_dict['redchisq']} using {self._n_gaussians} Gaussian components")
        logger.info(
            f"W10:                   {self._fit_dict['W10']} +/- {self._fit_dict['W10_e']}")
        logger.info(
            f"W50:                   {self._fit_dict['W50']} +/- {self._fit_dict['W50_e']}")
        logger.info(
            f"Wscat:                 {self._fit_dict['Wscat']} +/- {self._fit_dict['Wscat_e']}")
        logger.info(
            f"Weq:                   {self._fit_dict['Weq']} +/- {self._fit_dict['Weq_e']}")
        logger.info(f"Maxima:                {self._fit_dict['maxima']}")
        logger.info(f"Maxima error:          {self._fit_dict['maxima_e']}")
        logger.info(
            f"S/N estimate:          {self._fit_dict['sn']} +/- {self._fit_dict['sn_e']}")
        logger.debug(f"popt:                  {self._fit_dict['fit_params']}")

        # And we're done 

    # Plotter for the resulting fit

    def plot_fit(self):
        """Plots the best fit set of Gaussians"""
        if not self._plot_name:
            self._plot_name = "Gaussian_fit.png"
        popt = self._fit_dict["fit_params"]
        plot_x = np.linspace(0, 1, len(self._fit_profile))
        gaussian_x = np.linspace(-np.pi, np.pi,
                                 len(self._fit_profile), endpoint=False)
        plt.figure(figsize=(20, 12))

        z = self._exp_vm_gauss_conv(gaussian_x, *self._popt[0:5])
        plt_label = "Convolved Gaussian Component" if self._n_prof_components == 1 else f"Convolved Gaussian Component 1"
        plt.plot(plot_x, z, "--", label=plt_label)
        for i in range(1, self._n_prof_components):
            base = 0
            tau = self._popt[1]
            popt_idx = 5 + (i-1)*3
            z = self._exp_vm_gauss_conv(
                gaussian_x, base, tau, *self._popt[popt_idx:popt_idx + 3])
            plt.plot(plot_x, z, "--",
                     label=f"Convolved Gaussian Component {i+1}")

        start_j = 5 + (self._n_prof_components-1)*3
        for i, j in list(enumerate(range(start_j, len(self._popt), 3))):
            z = self._vm_gauss(gaussian_x, *self._popt[j:j+3])
            plt.plot(plot_x, z, "--", label=f"Gaussian Component {i}")
        if self._maxima:
            for i, mx in enumerate(self._maxima):
                plt.axvline(
                    x=(mx + self._maxima_e[i]), ls=":", lw=2, color="gray")
                plt.axvline(
                    x=(mx - self._maxima_e[i]), ls=":", lw=2, color="gray")

        plt.title(self._plot_name.split("/")[-1].split(".")[0], fontsize=26)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlim(0, 1)
        plt.xlabel("Bins", fontsize=24)
        plt.ylabel("Intensity", fontsize=24)
        plt.plot(plot_x, self._std_profile,
                 label="Original Profile", color="black")
        plt.plot(plot_x, self._fit_profile,
                 label="Gaussian Model", color="red")
        plt.legend(loc="best", prop={'size': 18})
        plt.savefig(self._plot_name, bbox_inches="tight")
        plt.close()

    def _standardise_raw_profile(self):
        """Normalises and rolls the raw profile to 0.25 in phase"""
        self._std_profile = self._raw_profile.copy()
        # Roll the profile such that the max point is at 0.25 phase
        roll_from = np.argmax(self._std_profile)
        roll_to = round(0.25*len(self._std_profile))
        roll_i = roll_to - roll_from
        self._std_profile = np.roll(self._std_profile, roll_i)
        # Normalise profile to max=1
        self._std_profile = np.array(
            self._std_profile)/max(self._std_profile)

    # The function that does the actual fitting

    def _fit_gaussian(self):
        """
        Fits multiple gaussian components to a pulse profile and finds the best number to use for a fit.
        Will always fit at least one gaussian per profile component.
        Profile components are defined by prof_utils.estimate_components_onpulse().
        Will fit a number of gaussians convoled with an exponential tail equal to the number of profile components
        Will then attempt to successively fit one more regular gaussian component until max_N is reached
        The number of gaussians fit that yields the lowes Bayesian Information Criterion is deemed the most successful
        """
        def add_bounds_guess(max_gen, width_gen, centre_gen, current_guess, current_bounds, base=False, tau=False):
            if base:
                # Base
                current_guess += [0.1]  # Guess at 10%
                current_bounds[0].append(0)
                # Noise should definitely not be above 50% max anyway
                current_bounds[1].append(0.5)
            if tau:
                # Tau
                current_guess += [1]
                current_bounds[0].append(1e-3)
                current_bounds[1].append(100)
            current_guess += [next(max_gen)*0.5**(len(current_guess)/5),
                              next(width_gen), next(centre_gen)]  # [Amp, kappa, loc]
            current_bounds[0].append(0)  # Amplitude
            current_bounds[1].append(1)
            current_bounds[0].append(1e-1)  # kappa
            current_bounds[1].append(1e5)
            current_bounds[0].append(-np.pi)  # Location
            current_bounds[1].append(np.pi)
            return current_guess, current_bounds

        # Initiate x - we will be working in 1D circular coordinates
        x = np.linspace(-np.pi, np.pi,
                        len(self._std_profile), endpoint=False)

        # Estimate gaussian parameters based on profile components
        comp_centres = []
        comp_max = []
        comp_width = []
        for pair in self._comp_est_dict["est_on_pulse"]:
            lower_bound = pair[0]
            upper_bound = pair[1]
            # Amplitude
            # Amp guess is half of max amplitude of component
            comp_max.append(
                max(self._std_profile[lower_bound:upper_bound])*0.5)
            # Gaussian width (kappa/sigma)
            width_est = (upper_bound - lower_bound) / \
                len(x)  # This is remeniscent of sigma
            width_est /= 2  # halve the width estimate
            comp_width.append(1/width_est**2)  # Convert this to kappa
            # Centre location
            # Looks like the middle of the profile component
            centre_val = np.max(
                max(self._std_profile[lower_bound:upper_bound]))
            centre_loc_in_prof = np.where(self._std_profile == centre_val)[
                0][0]  # get this in terms of profile index
            # Convert from 0 -> len(x) to -pi -> pi coordinates
            comp_centres.append((centre_loc_in_prof/len(x)) * 2*np.pi - np.pi)
        # Turn guesses into generators that cycle between profile components
        centre_gen = itertools.cycle(comp_centres)
        width_gen = itertools.cycle(comp_width)
        max_gen = itertools.cycle(comp_max)

        # Initiate bounds, guesses and fit dictionary
        bounds_arr = [[], []]
        guess = []
        fit_dict = {}

        # Add base and tau to the first guassian component for exp tail
        for i in range(0, self._n_prof_components):
            is_first = True if i == 0 else False  # Add the baseline to a single gaussian
            if i != self._n_prof_components:  # Don't add the rest of the guess to the last component, this is done in the next loop
                guess, bounds_arr = add_bounds_guess(
                    max_gen, width_gen, centre_gen, guess, bounds_arr, base=is_first, tau=is_first)

        # Fit from n_components to max_N gaussians to the profile. Evaluate profile fit using bayesian information criterion
        for num in range(self._n_prof_components, self._max_N):
            # Append to the bounds and guesses
            guess, bounds_arr = add_bounds_guess(
                max_gen, width_gen, centre_gen, guess, bounds_arr)
            # Bounds needs to be in tuple form for curve_fit
            bounds_tuple = (tuple(bounds_arr[0]), tuple(bounds_arr[1]))

            _guess = self._multi_gauss_exp_with_base(x, *guess)
            plt.figure(figsize=(12, 8))
            plt.plot(x, _guess)
            plt.plot(x, self._std_profile)
            plt.savefig(f"guess_plot_{num}.png")
            plt.close()

            # Do the fitting
            maxfev = 100000
            popt = None
            pcov = None
            try:
                popt, pcov = curve_fit(self._multi_gauss_exp_with_base, x,
                                       self._std_profile, bounds=bounds_tuple, p0=guess, maxfev=maxfev)
            except (RuntimeError, ValueError) as e:  # No fit found in maxfev
                raise BadFitError("No fit found in maxfev for this profile")

            # Work out some stuff
            if popt is not None and pcov is not None:
                fit = self._multi_gauss_exp_with_base(x, *popt)
                chisq = self._calculate_chisq(
                    self._std_profile, fit, self._noise_std)
                fit_dict[str(num+1)] = {"popt": [], "pcov": [],
                                        "fit": [], "chisq": []}
                fit_dict[str(num+1)]["popt"] = popt
                fit_dict[str(num+1)]["pcov"] = pcov
                fit_dict[str(num+1)]["fit"] = fit
                fit_dict[str(num+1)]["redchisq"] = chisq / \
                    (len(self._std_profile)-1)
                logger.debug(
                    f"Reduced chi squared for               {num+1} components: {fit_dict[str(num+1)]['redchisq']}")

        if not fit_dict:
            raise BadFitError(
                f"No suitable fit could be found for any provided values of N_gaussians")

        # Find the best fit according to reduced chi square
        best_chi_diff = np.inf
        for n_gaussians in fit_dict.keys():
            chi_diff = abs(1 - fit_dict[n_gaussians]["redchisq"])
            if chi_diff < best_chi_diff:
                best_chi_diff = chi_diff
                best_fit = n_gaussians

        self._popt = fit_dict[best_fit]["popt"]
        self._pcov = fit_dict[best_fit]["pcov"]
        self._fit_profile = fit_dict[best_fit]["fit"]
        self._best_chi = fit_dict[best_fit]["redchisq"]
        self._n_gaussians = int(best_fit)

    # Find min and max points

    def _find_minima_maxima_gauss(self):
        """
        Finds all roots of a gaussian function. Measured in Phase.
        """
        # Create the derivative list and spline it to find roots
        x = np.linspace(0, len(self._fit_profile)-1,  # Work in indexes
                        len(self._fit_profile))
        spline_prof = UnivariateSpline(x, self._fit_profile, s=0, k=4)
        dy_spline = spline_prof.derivative()
        dy_profile = dy_spline(x)
        roots = dy_spline.roots()
        plt.figure(figsize=(12, 8))
        plt.plot(self._fit_profile)
        plt.plot(np.array(dy_profile)/max(dy_profile))
        plt.savefig("stupid_shit_wont_work.png", bbox_inches="tight")
        plt.close()

        # Find which are max, min, and false
        maxima = []
        minima = []
        for root in roots:
            idx = round(root)
            # Set valid maxima to be 2 times noise
            if (self._fit_profile[idx] - self._popt[0]) > 2*self._noise_std:
                if dy_profile[idx-1] > dy_profile[idx]:
                    maxima.append(root)
                else:
                    minima.append(root)

        # Get errors
        minima_e = []
        maxima_e = []
        for m in minima:
            minima_e.append(error_in_x_pos(
                x, self._fit_profile, self._noise_std, round(m)))
        for m in maxima:
            maxima_e.append(error_in_x_pos(
                x, self._fit_profile, self._noise_std, round(m)))

        # Convert to phase, we're currently in indexes
        self._maxima = list(np.array(maxima)/len(self._fit_profile))
        self._minima = list(np.array(minima)/len(self._fit_profile))
        self._maxima_e = list(np.array(maxima_e)/len(self._fit_profile))
        self._minima_e = list(np.array(minima_e)/len(self._fit_profile))

    # Find the width of profile components

    def _find_widths(self):
        """
        Attempts to find the W_10, W_50 and equivalent width of a profile by using a spline approach.
        Weq errors are estimated by finding the average difference in Weq when you add and subtract the std from the on-pulse profile
        """

        # Perform spline operations on the fit
        x = np.linspace(0, len(self._fit_profile)-1,  # Work in indexes
                        len(self._fit_profile), endpoint=False)
        # Profile needs to be at zero baseline
        profile_for_widths = self._fit_profile - min(self._fit_profile)
        amp_fit = max(profile_for_widths) - min(profile_for_widths)
        spline10 = UnivariateSpline(
            x, profile_for_widths - np.full(len(x), 0.1*amp_fit), s=0)
        spline50 = UnivariateSpline(
            x, profile_for_widths - np.full(len(x), 0.5*amp_fit), s=0)
        spline_s = UnivariateSpline(
            x, profile_for_widths - np.full(len(x), 1/np.exp(1)*amp_fit), s=0)

        # Find Weq
        on_pulse_pairs = self._comp_est_dict["est_on_pulse"]
        on_pulse = []
        for pair in on_pulse_pairs:
            if pair[0] < pair[1]:
                on_pulse.extend(self._std_profile[pair[0]:pair[1]])
            else:
                on_pulse.extend(self._std_profile[pair[1]:])
                on_pulse.extend(self._std_profile[0:pair[0]])
        x_eq = np.linspace(0, len(on_pulse)-1, len(on_pulse))
        spline0 = UnivariateSpline(x_eq, on_pulse, s=0)
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
        err_10_1 = error_in_x_pos(
            x, self._fit_profile, self._noise_std, round(W10_roots[0]))
        err_10_2 = error_in_x_pos(
            x, self._fit_profile, self._noise_std, round(W10_roots[-1]))
        W10_e = np.sqrt(err_10_1**2 + err_10_2**2)

        # W50 root errors
        err_50_1 = error_in_x_pos(
            x, self._fit_profile, self._noise_std, round(W50_roots[0]))
        err_50_2 = error_in_x_pos(
            x, self._fit_profile, self._noise_std, round(W50_roots[-1]))
        W50_e = np.sqrt(err_50_1**2 + err_50_2**2)

        # Wscat root errors
        err_scat_1 = error_in_x_pos(
            x, self._fit_profile, self._noise_std, round(Wscat_roots[0]))
        err_scat_2 = error_in_x_pos(
            x, self._fit_profile, self._noise_std, round(Wscat_roots[-1]))
        Wscat_e = np.sqrt(err_scat_1**2 + err_scat_2**2)

        # Weq errors - using covariance formula
        on_pulse_less = (on_pulse - self._noise_std).clip(min=0)
        spline0 = UnivariateSpline(x_eq, on_pulse_less, s=0)
        integral = spline0.integral(0, len(self._std_profile)-1)
        dwdint = 1/max(on_pulse)**2
        dwdmax = -integral/max(on_pulse)**2
        int_e = abs(integral/max(on_pulse - self._noise_std) -
                    integral/max(on_pulse))
        max_e = self._noise_std
        Weq_e = np.sqrt(dwdint**2 * int_e**2 + dwdmax**2 *
                        max_e**2 + 2*dwdint*dwdmax*int_e*max_e)

        # Convert to phases
        self._W10 = W10/len(self._fit_profile)
        self._W10_e = W10_e/len(self._fit_profile)
        self._W50 = W50/len(self._fit_profile)
        self._W50_e = W50_e/len(self._fit_profile)
        self._Wscat = Wscat/len(self._fit_profile)
        self._Wscat_e = Wscat_e/len(self._fit_profile)
        self._Weq = Weq/len(self._fit_profile)
        self._Weq_e = Weq_e/len(self._fit_profile)

    def _est_sn(self):
        """
        Estimates the SN of a profile. Requires prof_utils.estimate_components_onpulse() to have been run
        on the profile and the noise stored in self._noise_std and component estimates in self._comp_est_dict
        """
        on_pulse_bins = 0
        for comp_range in self._comp_est_dict["overestimated_on"]:
            comp_beg = comp_range[0]
            comp_end = comp_range[1]
            if comp_beg < comp_end:
                on_pulse_bins += comp_end - comp_beg
            else:  # Component has wrapped around the phase bounds
                on_pulse_bins += len(self._std_profile) - comp_beg
                on_pulse_bins += comp_end
        self._n_off_pulse = len(self._std_profile) - on_pulse_bins
        self._sn = 1/self._noise_std
        # TODO: make this estimate better
        self._sn_e = 1/(self._noise_std * np.sqrt(2 * self._n_off_pulse - 2))

    
    ####################################
    ###### A bunch of maths stuff ######
    ####################################

    def _exp_tail(self, tau, x):
        return np.exp(-x/tau)/tau

    def _exp_vm_gauss_conv(self, x, base, tau, amp, kappa, loc):
        e_tail = self._exp_tail(tau/100, x)
        e_tail = e_tail/max(e_tail)
        g = self._vm_gauss(x, amp, kappa, loc)
        y = np.fft.irfft(np.fft.rfft(g) * np.fft.rfft(e_tail))  # Convolve
        y = amp*y/max(y)
        y += base
        return y

    def _vm_gauss(self, x, amp, kappa, loc):
        if kappa <= 10:  # Use vonmises
            g = vonmises.pdf(x, kappa, loc)
            g = amp * g/max(g)
        else:  # Use regular gaussian - vonmises can't handle large kappa
            sigma = np.sqrt(1/kappa)  # This is true for kappa approx > 10
            g = amp * np.exp(-(((x-loc)**2) / (2*sigma**2)))
        return g

    def _multi_gauss_exp_with_base(self, x, base, tau, *params):
        """
        Where x is a sequenetial numpy array of floats from -pi to pi of any length
        This function will employ a vonmises distribution where kappa is small and 
        convolve the first gaussian/vonmises with an exponential tail described by tau.
        base is the baseline noise
        The parameters are passed as a list of len%3 = 0 with values:
        amp: The amplitude of the gaussian
        kappa: The kappa value of the vonmises distribution. For k>10 we will use sigma=sqrt(1/k) instead
        loc: The location of the peak of the gaussian. Should be between -pi and pi
        """
        y = np.zeros_like(x)
        for i in range(0, len(params), 3):
            amp = params[i]
            kappa = params[i+1]
            loc = params[i+2]
            if i < self._n_prof_components:  # Make a gauss with exp tail
                if i != 0:  # Only use the base for one component
                    base = 0
                y += self._exp_vm_gauss_conv(x, base, tau, amp, kappa, loc)
            else:  # Make a vonmises/gaussian
                g = self._vm_gauss(x, amp, kappa, loc)
                y += g
        return y

    def _integral_multi_gauss(self, *params):
        y = 0
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
            y = y + a * np.exp(-(((x-b)**2) / (2*c**2)))
        return y

    def _multi_gauss_ddx(self, x, *params):
        # Derivative of gaussian
        y = np.zeros_like(x)
        for i in range(0, len(params), 3):
            a = params[i]
            b = params[i+1]
            c = params[i+2]
            y = y - a/c**2 * (x - b) * np.exp(-(((x-b)**2) / (2*c**2)))
        return y

    def _multi_gauss_d2dx2(self, x, *params):
        # Double derivative of gaussian
        y = np.zeros_like(x)
        for i in range(0, len(params), 3):
            a = params[i]
            b = params[i+1]
            c = params[i+2]
            y = y + (self._multi_gauss(x, a, b, c) / c**2) * \
                (((x - b)**2)/(c**2) - 1)
        return y

    def _partial_gauss_dda(x, a, b, c):
        return np.exp((-(b - x)**2)/(2*c**2))

    def _partial_gauss_ddb(x, a, b, c):
        return a*(x - b) * np.exp((-(b - x)**2)/(2*c**2))/c**2

    def _partial_gauss_ddc(x, a, b, c):
        return a*(x - b)**2 * np.exp((-(b - x)**2)/(2*c**2))/c**3

    # Chi sqaured evaluation

    def _calculate_chisq(self, observed_values, expected_values, err):
        test_statistic = 0
        for observed, expected in zip(observed_values, expected_values):
            test_statistic += ((float(observed)-float(expected))/float(err))**2
        return test_statistic