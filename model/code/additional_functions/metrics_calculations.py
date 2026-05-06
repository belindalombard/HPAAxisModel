

import numpy as np
from scipy.integrate import simpson
from numpy import trapezoid as trapz


'''Function to calculate the duration of the stressor'''
def calculate_stressor_duration(crh_original, crh_stressor, time_of_stressor, threshold = 0.01): 
    cut_original = crh_original[time_of_stressor:]
    cut_stressor = crh_stressor[time_of_stressor:]
    cut_original = np.array(cut_original)
    cut_stressor = np.array(cut_stressor)

    positions = np.where((cut_stressor-cut_original) < threshold)[0]

    return positions[0] if len(positions) > 0 else None


def calculate_amplitude_change(original_signal, stressed_signal, time_of_stressor):
    '''Calculate change in amplitude (max) of signal after introducing stressor'''
    during_stressor = stressed_signal[time_of_stressor:] - original_signal[time_of_stressor:]
    amplitude_change = np.abs(during_stressor).max()
    return amplitude_change


def calculate_phase_shift(original_signal, stressed_signal, time_of_stressor):

    original_peaks = np.where(np.diff(np.sign(np.diff(original_signal))) < 0)[0]
    stressed_peaks = np.where(np.diff(np.sign(np.diff(stressed_signal))) < 0)[0]
 
    return original_peaks, stressed_peaks


def difference_in_auc_trap(original_signal, stressed_signal, time_of_stressor): 
    area_original = trapz(original_signal, dx=1)
    area_stressed = trapz(stressed_signal, dx=1)
    return area_original-area_stressed

def difference_in_auc_simpson(original_signal, stressed_signal, time_of_stressor): 
    area_original = simpson(original_signal, dx=1)
    area_stressed = simpson(stressed_signal, dx=1)
    return area_original-area_stressed

def run_metrics(stressor_parameters, crh_values_original, crh_values, simulated_values_original, simulated_values): 
    results_dictionary = {}
    time_in_scope = int(stressor_parameters['time_in_scope'])
    results_dictionary['stressor_duration'] = calculate_stressor_duration(crh_values_original
                                                                          , crh_values
                                                                          , time_in_scope)
    results_dictionary['acth_amplitude'] = calculate_amplitude_change(simulated_values_original[:, 0], simulated_values[:, 0], time_in_scope)
    results_dictionary['cort_amplitude'] = calculate_amplitude_change(simulated_values_original[:, 1], simulated_values[:, 1], time_in_scope)
    
    acth_peaks, acth_peaks_shifted= calculate_phase_shift(simulated_values_original[:, 0]
                                                                , simulated_values[:, 0]
                                                                , time_in_scope)
    cort_peaks, cort_peaks_shifted = calculate_phase_shift(simulated_values_original[:, 1]
                                                                  , simulated_values[:, 1]
                                                                  , time_in_scope)
    
    results_dictionary['acth_auc_diff_simpson'] = difference_in_auc_simpson(simulated_values_original[:, 0], simulated_values[:, 0], time_in_scope)

    results_dictionary['acth_auc_diff_trap'] = difference_in_auc_trap(simulated_values_original[:, 0], simulated_values[:, 0], time_in_scope)

    results_dictionary['cort_auc_diff_simpson'] = difference_in_auc_simpson(simulated_values_original[:, 1], simulated_values[:, 1], time_in_scope)
    results_dictionary['cort_auc_diff_trap'] = difference_in_auc_trap(simulated_values_original[:, 1], simulated_values[:, 1], time_in_scope)

    return results_dictionary, acth_peaks, acth_peaks_shifted, cort_peaks, cort_peaks_shifted