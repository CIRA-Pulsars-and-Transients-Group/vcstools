from vcstools.prof_utils import ProfileLengthError, BadFitError
from vcstools.prof_utils import estimate_components_onpulse, error_in_x_pos, profile_region_from_pairs, filled_profile_region_between_pairs, normamlise_prof, get_off_pulse
import itertools
from scipy.stats import vonmises
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import logging
from scipy.interpolate import UnivariateSpline
from scipy.special import erf
import matplotlib
matplotlib.use('Agg')

# Error imports

logger = logging.getLogger(__name__)


# Main class
class gfit:
    """
    This class is used to fit multiple Gaussians to a pulse profile.

    Parameters
    ----------
    raw_profile : `list`
        The pulse profile to be fit.
    max_N : `int`
        The maximum number of gaussians to attempt to fit.
        |br| Default: 10.
    plot_name : `str`
        The name of the output plot. Can be set with gfit.plot_name. If unsupplied, will use a generic name.
        |br| Default: None.
    conponent_plot_name : `str`
        The name of the plot for the component plot that illustrates how profile components are chosen. If
        unsupplied, will not plot. Only applicable if on_pulse_range unsupplied.
        |br| Default: None.
    scattering_threshold : `float`
        The threshold for which any tau (scattering pulse width) value greater will be deemed scattered (in phase).
        |br| Default: 0.7.
    on_pulse_ranges : `list`
        A list of two-lists/tuples that describes the on pulse region in phase.
        e.g. [[0.1, 0.2], [0.6, 0.7]]
        If not supplied, will attempt to find on-pulse region
        |br| Default: None.

    Methods
    -------
    auto_fit:
        Runs a gaussian fit evaluation using up to max_N gaussians. The quality of the fit is decided via a chi-square evaluation.
        The best of these is saved and used to determine widths and maxima location. Part of the preprocessing of the profile will
        also determine the on-pulse regions and subsequently the noise level and signal to noise ratio.
    plot_fit:
        Plots the best chosen gaussian fit to a file whose name is gfit.plot_name. This can only be run after the fit_dict
        dictionary has been filled (presumably by auto_gfit).

    Returns
    -------
    fit_dict : `dict`
        contains the following keys:

        W10 : `float`
            The W10 width of the profile measured in phase.
        W10_e : `float`
            The uncertainty in the W10.
        W50 : `float`
            The W50 width of the profile measured in phase.
        W50_e : `float`
            The uncertainty in the W50.
        Weq : `float`
            The equivalent width of the profile measured in phase.
        Weq_e : `float`
            The uncertainty in the equivalent width.
        Wscat : `float`
            The scattering width of the profile measured in phaseprofile_region_from_pairs.
            The number of distinguishable profile components.
        on_pulse_estimates : `list`
            The output of prof_utils.estimate_components_onpulse(). Will be the rolled on_pulse_range input instead if supplied.
        fit_params : `list`
            A list of length 5 + 3*(N-1) there N is num_gauss. Each set of 3 parameters corresponds to the amp, centre and width of a guassian component.
            The first five are baseline, tau, amp, centre and width.
        cov_mat : np.matrix
            The covariance matrix output from curve_fit.
        profile : `list`
            The normalised and rolled input profile.
        fit : `list`
            The best fit made into a list form.
        sn : `float`
            The estimated signal to noise ratio, obtained from the profile.
        sn_e : `float`
            The uncertainty in sn.
        scattering : `float`
            The tau value of the profile fit.
        scattered : `boolean`
            Whether or not the final profile's tau value is greater than the sattering threshold.
    """

    def __init__(self, raw_profile, max_N=10, plot_name=None, component_plot_name=None, scattering_threshold=0.7, on_pulse_ranges=None):
        # Initialise inputs
        self._raw_profile = raw_profile
        self._max_N = max_N
        self._plot_name = plot_name
        self._component_plot_name = component_plot_name
        self._scattering_threshold = scattering_threshold
        self._on_pulse_ranges = on_pulse_ranges

        # The return dictionary
        self._fit_dict = {}

        # Stuff to put in the dictionary
        self._fit_profile = None
        self._best_chi = None
        self._std_profile = None
        self._on_pulse_prof = None
        self._noise_std = None
        self._noise_mean = None
        self._n_off_pulse = None
        self._n_gaussians = None
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
        return self._component_plot_name

    @component_plot_name.setter
    def component_plot_name(self, val):
        self._component_plot_name = val

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
        """
        Fits multiple gaussian profiles and finds the best combination of N_Gaussians and alpha.

        Uses
        """
        #if len(self._raw_profile) < 100:
        #    raise ProfileLengthError("Profile must have length > 100")
        if len(self._raw_profile) < 64:
            raise ProfileLengthError("Profile must have length > 64")

        if self._on_pulse_ranges:  # Use the on-pulse region to get noise estimates
            self._standardise_raw_profile(roll_phase=None)

            # Convert phases to bins
            on_pulse_ranges = []
            for phase_range in self._on_pulse_ranges:
                lower = int(phase_range[0] * len(self._std_profile))
                upper = int(phase_range[1] * len(self._std_profile))
                on_pulse_ranges.append([lower, upper])
            self._on_pulse_ranges = on_pulse_ranges

            # Get the off-pulse region and calculate noise
            off_pulse_range = get_off_pulse(self._on_pulse_ranges)
            off_pulse_region = profile_region_from_pairs(self._std_profile, off_pulse_range)
            self._noise_std = np.nanstd(off_pulse_region)
            self._noise_mean = np.nanmean(off_pulse_region)
            self._n_prof_components = len(self._on_pulse_ranges)

            # Roll the profile while preserving the known on-pulse region
            profile_bins = np.linspace(
                0, len(self._std_profile)-1, len(self._std_profile))
            roll_i = self._roll_to_ideal_phase(off_pulse_range)
            rolled_profile_bins = np.roll(profile_bins, -roll_i)
            on_pulse_ranges = []
            for phase_range in self._on_pulse_ranges:
                lower = rolled_profile_bins[phase_range[0]]
                upper = rolled_profile_bins[phase_range[1]]
                on_pulse_ranges.append([int(lower), int(upper)])
            self._on_pulse_ranges = on_pulse_ranges
            #quit()

        else:  # We need to estimate what the on-pulse region is
            self._standardise_raw_profile()
            comp_est_dict = estimate_components_onpulse(self._std_profile)
            self._n_prof_components = len(comp_est_dict["overest_on_pulse"])

            # Roll the profile such that an off-pulse region is on the phase edge
            self._roll_to_ideal_phase(comp_est_dict["overest_off_pulse"])

            # This has ruined our component estimates so we need get them back
            comp_est_dict = estimate_components_onpulse(
                self._std_profile, plot_name=self._component_plot_name)
            self._on_pulse_ranges = comp_est_dict["overest_on_pulse"]
            self._n_prof_components = len(comp_est_dict["overest_on_pulse"])
            self._noise_std = comp_est_dict["norm_noise_std"]
            self._noise_mean = comp_est_dict["norm_noise_mean"]

        # Create the on-pulse profile where the noise is set to the mean noise value
        self._on_pulse_prof = filled_profile_region_between_pairs(
            self._std_profile, self._on_pulse_ranges, fill_value=self._noise_mean)

        self._on_pulse_bool = []
        for opp in self._on_pulse_prof:
            if opp == self._noise_mean:
                self._on_pulse_bool.append(False)
            else:
                self._on_pulse_bool.append(True)

        # Attempt to fit gaussians to the profile
        self._fit_gaussian()  # Initial fit
        self._est_noise_sn()  # Initial Estimate SN estimate from the fit
        # Force no scattering if the fit isn't scattered now
        no_scat = True if self._scattered and self._popt[1] <= self._scattering_threshold else False
        # Redo the fit with updated scattering and noise information
        self._fit_gaussian(force_no_scattering=no_scat)

        # Do all the fancy things with the fit
        self._find_widths()  # Find widths + error
        self._find_minima_maxima_gauss()  # Find turning points
        self._est_noise_sn()  # Estimate SN

        # Update scattering information
        self._scattered = True if self._Wscat >= self._scattering_threshold else False

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
        self._fit_dict["minima"] = self._minima
        self._fit_dict["minima_e"] = self._minima_e
        self._fit_dict["redchisq"] = self._best_chi
        self._fit_dict["n_components"] = self._n_gaussians
        self._fit_dict["fit_params"] = self._popt
        self._fit_dict["cov_mat"] = self._pcov
        self._fit_dict["profile"] = self._std_profile
        self._fit_dict["on_pulse_bool"] = self._on_pulse_bool
        self._fit_dict["fit"] = self._fit_profile
        self._fit_dict["sn"] = self._sn
        self._fit_dict["sn_e"] = self._sn_e
        self._fit_dict["noise_std"] = self._noise_std
        self._fit_dict["noise_mean"] = self._noise_mean
        self._fit_dict["scattering"] = self._popt[1]  # tau
        self._fit_dict["scattered"] = self._scattered
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
        logger.info(
            f"S/N:                   {self._fit_dict['sn']} +/- {self._fit_dict['sn_e']}")
        logger.info(
            f"Noise estiamte:        {self._noise_std}")
        logger.info(
            f"Scattering Timescale:  {self._popt[1]}")
        logger.info(f"popt:                  {self._fit_dict['fit_params']}")

        # And we're done

    # Plotter for the resulting fit

    def plot_fit(self):
        """Plots the best fit set of Gaussians.

        Default output is "Gaussian_fit.png".
        """
        if not self._plot_name:
            self._plot_name = "Gaussian_fit.png"
        popt = self._fit_dict["fit_params"]
        # We need an x range for the plot - which will be phase and for
        # Recreating the Gaussians - which will be -pi to pi
        plot_x = np.linspace(0, 1, len(self._fit_profile), endpoint=False)
        gaussian_x = np.linspace(-np.pi, np.pi,
                                 len(self._fit_profile), endpoint=False)
        fig = plt.figure(figsize=(20, 12))

        # Plot the convolved gaussian
        z = self._exp_vm_gauss_conv(gaussian_x, *self._popt[0:5])
        plt_label = "Convolved Gaussian Component" if self._n_prof_components == 1 else f"Convolved Gaussian Component 1"
        plt.plot(plot_x, z, "--", label=plt_label)

        # Plot any other convoled gaussians - one for each component past the first
        for i in range(1, self._n_prof_components):
            base = 0
            tau = self._popt[1]
            popt_idx = 5 + (i-1)*3
            z = self._exp_vm_gauss_conv(
                gaussian_x, base, tau, *self._popt[popt_idx:popt_idx + 3])
            plt.plot(plot_x, z, "--",
                     label=f"Convolved Gaussian Component {i+1}")

        # Plot the other gaussian components
        start_j = 5 + (self._n_prof_components-1)*3
        for i, j in list(enumerate(range(start_j, len(self._popt), 3))):
            z = self._vm_gauss(gaussian_x, *self._popt[j:j+3])
            plt.plot(plot_x, z, "--", label=f"Gaussian Component {i+1}")
        if self._maxima:
            for mx in self._maxima:
                plt.axvline(x=mx, ls=":", lw=2, color="gray")

        # Error bar - for noise and SN
        plt.text(0.03, 0.94,
                 f"S/N:   {round(self._sn, 2)} +/- \n          {round(self._sn_e, 3)}", fontsize=14)
        plt.text(
            0.03, 0.91 - self._noise_std/2, f"Noise: {round(self._noise_std, 5)}", fontsize=14)
        plt.errorbar(0.02, 0.91 - self._noise_std/2, yerr=self._noise_std, color="gray",
                     capsize=10, markersize=0, label="Noise")

        # All the other plotting stuff
        formats = tuple(list(fig.canvas.get_supported_filetypes().keys()))
        for f in formats:
            if self._plot_name.endswith(f):
                title = self.plot_name[:-len(f)-1]
        if not self.plot_name.endswith(formats):
            title = self._plot_name
            self._plot_name += ".png"

        plt.title(title, fontsize=26)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlim(0, 1)
        plt.xlabel("Phase", fontsize=24)
        plt.ylabel("Intensity", fontsize=24)
        plt.plot(plot_x, self._std_profile,
                 label="Original Profile", color="black")
        plt.plot(plot_x, self._fit_profile,
                 label="Gaussian Model", color="red")
        plt.legend(loc="best", prop={'size': 14})
        plt.savefig(self._plot_name, bbox_inches="tight")
        plt.close()


    def _standardise_raw_profile(self, roll_phase=0.25):
        """Normalises and rolls the raw profile to 0.25 in phase"""
        self._std_profile = self._raw_profile.copy()
        # Roll the profile such that the max point is at 0.25 phase (default)
        if roll_phase:
            roll_from = np.argmax(self._std_profile)
            roll_to = round(roll_phase*len(self._std_profile))
            roll_i = roll_to - roll_from
            self._std_profile = np.roll(self._std_profile, roll_i)
        # Normalise profile to max=1, min=0
        self._std_profile = normamlise_prof(self._std_profile)

    def _roll_to_ideal_phase(self, off_pulse_ranges):
        """Finds the largest off-pulse region from the off pulse ranges and rolls the profile to that index"""
        # Find the largest off-pulse range
        off_pulse_lens = []
        for off_range in off_pulse_ranges:
            start_off = off_range[0]
            end_off = off_range[1]
            if start_off > end_off:  # Wrapped around phase end
                off_pulse_lens.append(
                    len(self._std_profile) - start_off + end_off)
            else:
                off_pulse_lens.append(end_off - start_off)

        longest_range, _ = sorted(zip(off_pulse_ranges, off_pulse_lens), reverse=True)[
            0]  # Grabs the longest one
        start_off, end_off = longest_range
        # Find the centre index of this range
        if start_off > end_off:  # off pulse is wrapped over phase end
            off_pulse_len = (len(self._std_profile) - start_off + end_off)
            centre_of_off_pulse = start_off + off_pulse_len/2
            if centre_of_off_pulse > len(self._std_profile):
                centre_of_off_pulse = centre_of_off_pulse - \
                    len(self._std_profile)
            roll_from = int(centre_of_off_pulse)
        else:  # Easy case - off pulse is not wrapped
            roll_from = int((end_off - start_off)/2 + start_off)

        # Do the rolling
        roll_to = len(self._std_profile)
        roll_i = roll_to - roll_from
        self._std_profile = np.roll(self._std_profile, roll_i)
        return roll_i

    # The function that does the actual fitting

    def _fit_gaussian(self, force_no_scattering=False):
        """
        Fits multiple gaussian components to a pulse profile and finds the best number to use for a fit.
        Will always fit at least one gaussian per profile component.
        Profile components are defined by prof_utils.estimate_components_onpulse().
        Will fit a number of gaussians convoled with an exponential tail equal to the number of profile components
        Will then attempt to successively fit one more regular gaussian component until max_N is reached
        The number of gaussians fit that yields the lowes Bayesian Information Criterion is deemed the most successful
        """
        # Initiate x - we will be working in 1D circular coordinates
        x = np.linspace(-np.pi, np.pi,
                        len(self._std_profile), endpoint=False)

        # Estimate gaussian parameters based on profile components
        comp_locs = []
        comp_amp = []
        comp_width = []
        components = []
        for pair in self._on_pulse_ranges:
            lower_bound = pair[0]
            upper_bound = pair[1]
            # Amplitude
            # Amp guess is half of max amplitude of component
            on_pulse_region = profile_region_from_pairs(
                self._std_profile, [pair])
            comp_amp.append(max(on_pulse_region)*0.5)

            # Gaussian width (kappa/sigma)
            width_est = len(profile_region_from_pairs(
                self._std_profile, [pair]))/2  # halve the width est because it always looks too big
            width_est = 2*np.pi * width_est/len(self._std_profile)
            width_est = 1/width_est**2  # Convert this to kappa
            comp_width.append(width_est)

            # The entire component profile will be useful to know about
            components.append(profile_region_from_pairs(
                np.linspace(-np.pi, np.pi, len(self._std_profile), endpoint=False), [pair]))

        # Centre location - loop over on-pulse. This estimate only makes sense if the edge of the phase is
        # off-pulse. Which should always be the case by this point.
        for on_range in self._on_pulse_ranges:
            middle = (on_range[0] + on_range[1])/2
            # Convert m from 0 -> len(x) to -pi -> pi coordinates
            comp_locs.append((middle/len(x)) * 2*np.pi - np.pi)

        # Turn guesses into generators that cycle between profile components
        amp_gen = itertools.cycle(comp_amp)
        width_gen = itertools.cycle(comp_width)
        loc_gen = itertools.cycle(comp_locs)
        components_gen = itertools.cycle(components)

        # Check if the profile is scattered - this will determine which profile to use for curve_fit()
        if not force_no_scattering:
            scattering_bounds = [[], []]
            scattering_guess = []
            amp_gen_sc = itertools.cycle(comp_amp)
            width_gen_sc = itertools.cycle(comp_width)
            loc_gen_sc = itertools.cycle(comp_locs)
            components_gen_sc = itertools.cycle(components)
            for i in range(self._n_prof_components):
                # Add the baseline and tau to a single gaussian
                is_first = True if i == 0 else False
                scattering_guess, scattering_bounds = self._add_bounds_guess(
                    amp_gen_sc, width_gen_sc, loc_gen_sc, components_gen_sc, scattering_guess, scattering_bounds, base=is_first, tau=is_first)
            # Bounds needs to be in tuple form for curve_fit
            scattering_bounds_tuple = (
                tuple(scattering_bounds[0]), tuple(scattering_bounds[1]))
            try:
                popt, _ = curve_fit(self._multi_gauss_exp_with_base, x, self._std_profile,
                                    bounds=scattering_bounds_tuple, p0=scattering_guess, maxfev=10000)
                # TODO: make this estimate better
                self._scattered = True if popt[1] >= self._scattering_threshold else False
            except (RuntimeError, ValueError) as e:
                self._scattered = True
            # The profile to fit - scattered profiles will use the entire standardised profile as we expect that the
            # on-pulse region is not accurate
            fitting_profile = self._std_profile if self._scattered else self._on_pulse_prof
        else:
            fitting_profile = self._on_pulse_prof

        # Do the fittingloop
        fit_dict = self._fit_loop(
            fitting_profile, amp_gen, width_gen, loc_gen, components_gen)

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

        # Record the best things in the object
        self._popt = fit_dict[best_fit]["popt"]
        self._pcov = fit_dict[best_fit]["pcov"]
        self._fit_profile = fit_dict[best_fit]["fit"]
        self._best_chi = fit_dict[best_fit]["redchisq"]
        self._n_gaussians = int(best_fit)

    def _fit_loop(self, profile, amp_gen, width_gen, loc_gen, components_gen, reduce_chi_threshold=0.05):
        """
        A loop that attempts to fit multiple gaussians, some with exponential tails, to a profile.
        If scattered=True,
        """
        # Initiate x - we will be working in 1D circular coordinates
        x = np.linspace(-np.pi, np.pi,
                        len(self._std_profile), endpoint=False)

        # Initiate bounds, guesses and fit dictionary
        bounds_arr = [[], []]
        guess = []
        fit_dict = {}

        # Fit from n_components to max_N gaussians to the profile. Evaluate profile fit using bayesian information criterion
        for num in range(self._max_N):
            # Add the baseline and tau to a single gaussian
            is_first = True if num == 0 else False
            guess, bounds_arr = self._add_bounds_guess(
                amp_gen, width_gen, loc_gen, components_gen, guess, bounds_arr, base=is_first, tau=is_first)
            # Bounds needs to be in tuple form for curve_fit
            bounds_tuple = (tuple(bounds_arr[0]), tuple(bounds_arr[1]))

            if num+1 < self._n_prof_components:  # Skip until we have at least num = profile_components
                continue

            # Do the fitting
            popt = None
            pcov = None
            try:
                popt, pcov = curve_fit(self._multi_gauss_exp_with_base, x,
                                       profile, bounds=bounds_tuple, p0=guess, maxfev=10000)
            except (RuntimeError, ValueError) as e:  # No fit found in maxfev
                logger.warn(e)
                for i, g in enumerate(guess):
                    logger.debug(
                        f"Guess: {g} | bounds: {bounds_tuple[0][i]} -> {bounds_tuple[1][i]}")
                # continue
                break

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
                    (len(self._std_profile)-len(popt))
                logger.debug(
                    f"Reduced chi squared for               {num+1} components: {fit_dict[str(num+1)]['redchisq']}")

                # function exit condition
                if abs(fit_dict[str(num+1)]["redchisq"]-1) <= reduce_chi_threshold:
                    return fit_dict
        return fit_dict

    def _add_bounds_guess(self, amp_gen, width_gen, loc_gen, component_gen, current_guess, current_bounds, base=False, tau=False):
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
            current_bounds[1].append(200)
        # Add the guesses
        amp_guess = next(amp_gen)
        width_guess = next(width_gen)
        loc_guess = next(loc_gen)
        component_guess = next(component_gen)
        # Add the bounds
        amp_lo = 3*self._noise_std if self._noise_std < 0.2 else 0.6
        amp_hi = 1.0
        # width bounds is hard - need to reverse engineer how we got here in the first place
        # kappa_sigma (width_est) = (2pi * og_width / prof_len)^-2
        # og_width = (kappa_sigma^-1/2 * prof_len) / 2pi
        og_width = (width_guess**(-1/2) * len(self._std_profile)) / (2*np.pi)
        # Profiles are not allowed to be much wider than their component estimates
        width_lo = og_width * 3
        width_lo = (2*np.pi * width_lo / len(self._std_profile))**-2
        width_hi = 1e6  # Components are allowed to be very sharp
        # Locations should be confined to their respective component
        loc_lo = component_guess[0]
        loc_hi = component_guess[-1]
        if loc_lo > loc_hi:
            loc_lo = -np.pi
            loc_hi = np.pi

        # Fix guesses if they're out of bounds
        amp_guess = amp_lo if amp_guess < amp_lo else amp_guess
        amp_guess = amp_hi if amp_guess > amp_hi else amp_guess
        width_guess = width_lo if width_guess < width_lo else width_guess
        width_guess = width_hi if width_guess > width_hi else width_guess
        loc_guess = loc_lo if loc_guess < loc_lo else loc_guess
        loc_guess = loc_hi if loc_guess > loc_hi else loc_guess
        # [Amp, kappa, loc]
        current_guess += [amp_guess, width_guess, loc_guess]
        # Add the bounds
        current_bounds[0].append(amp_lo)  # Amplitude
        current_bounds[1].append(amp_hi)
        current_bounds[0].append(width_lo)  # kappa
        current_bounds[1].append(width_hi)
        current_bounds[0].append(loc_lo)  # Location
        current_bounds[1].append(loc_hi)
        return current_guess, current_bounds

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
        d2y_spline = UnivariateSpline(
            x, self._fit_profile, s=0, k=5).derivative().derivative()
        d2y_profile = d2y_spline(x)
        roots = dy_spline.roots()
        logger.debug(roots)
        roots = [int(round(i)) for i in roots]
        logger.debug(roots)
        roots = list(set(roots)) # Remove duplicates
        logger.debug(roots)

        # Find which are max, min, and false
        maxima = []
        minima = []
        for root in roots:
            # Set valid maxima to be 2 times noise
            if (self._fit_profile[root] - self._popt[0]) > 2*self._noise_std:
                if d2y_profile[root] < 0:
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
        on_pulse = []
        for pair in self._on_pulse_ranges:
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
        W10_roots = [int(round(i)) for i in W10_roots]
        W50_roots = spline50.roots()
        W50_roots = W50_roots = [int(round(i)) for i in W50_roots]
        Wscat_roots = spline_s.roots()
        Wscat_roots = [int(round(i)) for i in Wscat_roots]
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

    def _est_noise_sn(self):
        """
        Estimates the SN of a profile using the fit profile to determine an off-pulse region.
        """
        sorted_profile = [i for _, i in sorted(
            zip(self._fit_profile, self._std_profile), key=lambda pair: pair[0])]
        ten_percent = round(0.1*len(self._std_profile))
        self._noise_std = np.std(sorted_profile[:ten_percent])
        self._noise_mean = np.mean(sorted_profile[:ten_percent])
        self._sn = 1/self._noise_std
        # TODO: make this estimate better
        self._sn_e = 1/(self._noise_std * np.sqrt(2 * ten_percent - 2))

    ####################################
    ###### A bunch of maths stuff ######
    ####################################

    def _exp_tail(self, tau, x):
        return np.exp(-x/tau)/tau

    def _exp_vm_gauss_conv(self, x, base, tau, amp, kappa, loc):
        e_tail = self._exp_tail(tau/100, x)
        e_tail = e_tail/max(e_tail)
        g = self._vm_gauss(x, amp, kappa, loc)
        y = np.real(np.fft.ifft(np.fft.fft(g) *
                    np.fft.fft(e_tail)))  # Convolve
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
            if i < 3*self._n_prof_components:  # Make a gauss with exp tail
                if i != 0:  # Only use the base for one component
                    base = 0
                conv = self._exp_vm_gauss_conv(x, base, tau, amp, kappa, loc)
                y += conv
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
