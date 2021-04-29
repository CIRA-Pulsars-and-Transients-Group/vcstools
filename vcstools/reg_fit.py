import numpy as np
import logging
from scipy.interpolate import UnivariateSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from vcstools.prof_utils import est_sn_from_prof, sigmaClip, check_clip
from vcstools.stickel import Stickel
# Error imports
from vcstools.prof_utils import LittleClipError, LargeClipError, NoComponentsError, ProfileLengthError, NoFitError, BadFitError

logger = logging.getLogger(__name__)

# Main class
class reg_fit:
    """
    This class is used to fit a pulse profile using a regularisation technique by Stickel 2010

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
    clip_type: str:
        The verbosity of clipping. Choose between 'regular', 'noisy' or 'verbose'.
        Default: 'regular'
    """

    def __init__(self, raw_profile, plot_name=None, min_comp_len=None):
        self._raw_profile = raw_profile
        #self._lambda_range = lambda_range
        self._plot_name = plot_name
        self._min_comp_len = min_comp_len
        #self._clip_type = clip_type

        self._fit_dict = {} # Use this to store the solution
        self._best_chi = None
        self._best_alpha = None
        self._best_lambda = None

        # Private stuff
        self._current_std_prof = []
        self._current_alpha = None
        self._current_lambda = None
        self._current_chi = None
        self._current_noise_std = None
        self._current_noise_mean = None
        self._clipped_prof = []
        self._current_off_pulse_prof = []
        self._current_n_off_pulse = None
        self._current_reg_prof = None

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
    def lambda_range(self):
        return self.lambda_range
    @lambda_range.setter
    def lambda_range(self, val):
        self.lambda_range = val

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
    def current_std_profile(self):
        return self._current_std_prof
    @current_std_profile.setter
    def current_std_profile(self, val):
        self._current_std_prof = val

    @property
    def current_alpha(self):
        return self._current_alpha
    @current_alpha.setter
    def current_alpha(self, val):
        self._current_alpha = val

    @property
    def current_noise_std(self):
        return self._current_noise_std
    @current_noise_std.setter
    def current_noise_std(self, val):
        self._current_noise_std = val

    @property
    def current_clipped_prof(self):
        return self._current_clipped_prof
    @current_clipped_prof.setter
    def current_clipped_prof(self, val):
        self._current_clipped_prof = val

    @property
    def current_off_pulse_prof(self):
        return self._current_off_pulse_prof
    @current_clipped_prof.setter
    def current_off_pulse_prof(self, val):
        self._current_off_pulse_prof = val


    def auto_reg_fit(self):
        # One solution per lambda/alpha combination
        lambda_range = np.geomspace(1e-10, 1e-3, 8)
        #alpha_range = np.linspace(1, 5, 33)
        alpha_range = np.linspace(1, 5, 17)
        #alpha_range = np.linspace(1, 5, 9)
        # Set the 'best' trackers
        self._best_alpha = None
        self._best_lambda = None
        self._best_chi = np.inf # Reduced Chi squared is our metric 
        
        for a in alpha_range:
            self._current_alpha = a
            print(f"alpha: {a}")
            # Standardise the profile
            try:
                self._standardise_raw_profile()
            except(LittleClipError, LargeClipError, NoComponentsError, ProfileLengthError, BadFitError):
                logger.debug(f"Skipping alpha: {a}")
                print(f"Skipping alpha: {a}")
                continue # This is due to the alpha value so we can skip all lambdas
            for l in lambda_range:
                self._current_lambda = l
                # Try with this a and l combination
                self._reg_fit_attempt()
                # Check the reduced chi squared value
                print(f"lambda: {l} | chi: {self._current_chi}")
                if abs(1-self._current_chi) < abs(1-self._best_chi):
                    self._best_alpha = a
                    self._best_lambda = l
                    self._best_chi = self._current_chi

        # Now that we have determined the best alpha and lambda combination, we can do some more analysis
        self._current_alpha = self._best_alpha
        self._current_lambda = self._best_lambda
        print(f"Using alpha = {self._best_alpha} | lambda = {self._best_lambda} | chi = {self._best_chi}")
        self._standardise_raw_profile()
        self._reg_fit_attempt()
        plt.figure(figsize=(20, 12))
        plt.plot(self._current_std_prof, label="std_prof")
        plt.plot(self._current_clipped_prof, label="clipped_prof")
        plt.plot(self._current_reg_prof, label="reg_prof")
        plt.legend()
        plt.savefig(f"{self._plot_name}")
        plt.close()

            
    def _reg_fit_attempt(self):
        """Fits a profile using a regularisation technique"""
        # Chi sqaured evaluation
        def chi_sq(observed_values, expected_values, err):
            test_statistic=0
            for observed, expected in zip(observed_values, expected_values):
                test_statistic+=((float(observed)-float(expected))/float(err))**2
            return test_statistic
        
        # Create a regularised version of the clipped profile
        # Stickel takes a 2-D array
        x = np.linspace(0, 1, len(self._current_clipped_prof))
        xy_prof = []
        for i, val in enumerate(self._current_clipped_prof):
            xy_prof.append([x[i], val])
        xy_prof = np.array(xy_prof)
        stickel = Stickel(xy_prof)
        # Regularise
        stickel.smooth_y(self._current_lambda)
        self._current_reg_prof = stickel.yhat

        # Evaluate fit using reduced chi squared
        self._current_chi = chi_sq(self._current_std_prof, self._current_reg_prof, self._current_noise_std)/(len(self._current_std_prof)-1)
        #self._current_chi = chi_sq(self._current_std_prof, self._current_reg_prof, self._current_noise_mean)/(len(self._current_std_prof)-1)


    def _standardise_raw_profile(self):
        """Clips, normalises and fills in the raw profile"""
        self._current_std_prof = self._raw_profile.copy()
        # Roll the profile such that the max point is at 0.25 phase
        roll_from = np.argmax(self._current_std_prof)
        roll_to = round(0.25*len(self._current_std_prof))
        roll_i = roll_to - roll_from
        self._current_std_prof = np.roll(self._current_std_prof, roll_i)
        # Normalise profile to max=1
        self._current_std_prof = np.array(self._current_std_prof)/max(self._current_std_prof)
        # Make a signal-clipped profile to determine which bins are noise vs on-pulse
        _, self._current_off_pulse_prof = sigmaClip(self._current_std_prof, alpha=self._current_alpha)
        # Check the clip to fill in gaps and then remove standalone on-pulse regions
        check_clip(self._current_off_pulse_prof)
        self._fill_off_pulse_prof()
        self._clean_small_noise()
        self._remove_noisy_components()
        # Remove the noise such that the on-pulse profile is centered about 0   
        self._current_noise_std = np.nanstd(np.array(self._current_off_pulse_prof))     
        self._current_noise_mean = np.nanmean(np.array(self._current_off_pulse_prof))
        #self._current_std_prof = self._current_std_prof - self._current_noise_std
        self._current_std_prof = self._current_std_prof - self._current_noise_mean
        # Renormalise
        self._current_std_prof = self._current_std_prof/max(self._current_std_prof)
        # Remove noise from std_profile - make new array for this
        self._current_clipped_prof = self._current_std_prof.copy()
        for i, _ in enumerate(self._current_off_pulse_prof):
            if not np.isnan(self._current_off_pulse_prof[i]):
                self._current_clipped_prof[i] = 0
        """
        plt.figure(figsize=(20, 12))
        plt.plot(self._current_std_prof, label="std")
        plt.plot(self._current_clipped_prof, label="clipped")
        plt.savefig("test.png")
        plt.close()
        quit()   
        """

    # SigmaClip isn't perfect. Use these next function to fix bad clips
    def _fill_off_pulse_prof(self, search_scope=None):
        """
        Intended for use on noisy profiles. Fills nan values that are surrounded by non-nans to avoid discontinuities
        in the profile
        """
        length = len(self._current_off_pulse_prof)
        if search_scope is None:
            # Search 5% ahead for non-nans
            search_scope = int(np.ceil(np.sqrt(length*0.1)))
        search_scope = np.linspace(search_scope, 1, search_scope, dtype=int)
        # Search scope = [x, x-1, x-2, ..., 1]

        # Loop over all values in clipped profile
        for i, val in enumerate(self._current_off_pulse_prof):
            # If current index is on-pulse
            fill_indexes = []
            if np.isnan(val):
                found_gaps = []
                # Loop over search scope
                for search_i in search_scope: 
                    # Make a copy of the profile and roll back 
                    rolled_prof = np.roll(self._current_off_pulse_prof, -search_i)
                    # If we've found a gap and an on-pulse near it, break
                    if np.isnan(rolled_prof[i]) and found_gaps:
                        break
                    elif not np.isnan(rolled_prof[i]):
                        found_gaps.append(i + search_i)
                fill_indexes += found_gaps

        fill_indexes = list(set(fill_indexes)) # Remove duplicates
        fill_indexes = [i%len(self._current_off_pulse_prof) for i in fill_indexes] # Wrap around 
        # Fill in the on-pulse
        for i in fill_indexes:
            self._current_off_pulse_prof[i] = np.nan


    def _clean_small_noise(self):
        """Removes pieces of signal that haven't been flagged by sigmaClip"""
        min_len = int(np.ceil(np.sqrt(len(self._current_off_pulse_prof)*0.05))) # Sqrt of 5% profile length rounded up
        my_prof = self._current_off_pulse_prof.copy()
        for i, val in enumerate(my_prof):
            if np.isnan(val): # On pulse
                # Check the left and right of the current bin to see how large the component is
                left = 0
                right = 0
                try_left = True
                try_right = True
                test_val = np.nan
                for j in range(min_len):
                    if np.isnan(np.roll(my_prof, j)[i]) and try_left: # check if j spaces to the left is on-pulse
                        left += 1
                    else: try_left = False
                    if np.isnan(np.roll(my_prof, -j)[i]) and try_right: # check if j spaces to the right is on-pulse
                        right += 1
                    else: try_right = False
                if left + right + 1 < min_len: # If this component is too small, mark it as off-pulse
                    my_prof[i] = self._current_std_prof[i]    
        self._current_off_pulse_prof = my_prof # Set the off-pulse to this modified profile
    

    def _remove_noisy_components(self):
        """Removes the components that look too much like noise"""
        # First find the components
        comp_dict, comp_idxs = self._find_components()

        # Define what looks like noise
        noise_std = np.nanstd(self._current_off_pulse_prof)
        noise_mean = np.nanmean(self._current_off_pulse_prof)
        std_tolerance = 3*noise_std
        mean_tolerance = 3*noise_mean

        # Check if each component looks like noise
        components_to_remove = []
        for key in comp_dict.keys():
            comp = comp_dict[key]
            comp_std = np.std(comp)
            comp_mean = np.mean(comp)
            if abs(comp_std - noise_std) <= std_tolerance and abs(comp_mean - noise_mean) <= mean_tolerance:
                components_to_remove.append(key)
        
        # Remove components
        for comp_key in components_to_remove:
            for i in comp_idxs[comp_key]:
                self._current_off_pulse_prof[i] = self._current_std_prof[i]


    def _find_components(self):
        """
        Given a profile in which the on-pulse clipped to nans, finds the components that are clumped together.

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
        for i, val in enumerate(self._current_off_pulse_prof):
            if np.isnan(val): # On pulse
                if not np.isnan(self._current_off_pulse_prof[i-1]) or i==0: # Off pulse or beginning of profile phase
                    num_components+=1
                    comp_key = f"component_{num_components}"
                    component_dict[comp_key]=[]
                    component_idx[comp_key]=[]
                component_dict[comp_key].append(self._current_std_prof[i])
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


