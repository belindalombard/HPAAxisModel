import pints 
import numpy as np
import scipy
from scipy.signal import hilbert
from datetime import datetime

class ErrorMeasureTwoSignal(pints.ErrorMeasure):
    def __init__(self, model, times, values):
        self.model = model
        self.times = times
        self.values = values
        self.last_simulated_values = None
        self.current_parameters = None
        self.current_error = 0

    def n_parameters(self):
        return self.model.n_parameters()

    def __call__(self, parameters):
        self.current_parameters = parameters
        self.last_simulated_values = self.model.simulate(parameters, self.times)

        normalised_results = [
            self.normalise_array(self.last_simulated_values[:, i], self.values[:, i])
            for i in range(self.last_simulated_values.shape[1])
        ]
        normalised_sim = np.array([normalised_results[0][0], normalised_results[1][0]])
        normalised_obs = np.array([normalised_results[0][1], normalised_results[1][1]])
        base_error = float(np.sum((normalised_sim - normalised_obs) ** 2))

        # Pure MSE error, no penalty
        error = base_error
        self.current_error = error
        return error

    def normalise_array(self, simulated, observed):
        max_sim = np.max(simulated)
        max_obs = np.max(observed)
        max_value = max(max_sim, max_obs)
        if max_value == 0:
            return np.zeros_like(simulated, dtype=float), np.zeros_like(observed, dtype=float)
        return simulated / max_value, observed / max_value

class PhaseAlignedMSE(pints.ErrorMeasure):
    """
    Phase-invariant MSE wrapper:
    - Simulate once with the current parameters.
    - Find best circular shift between simulated and observed (default ACTH).
    - Roll the simulated series by that shift.
    - Compute the same normalised SSE as ErrorMeasureTwoSignal.
    """
    def __init__(self, model, times, values, align_channel=0):
        self.model = model
        self.times = times
        self.values = values
        self.align_channel = int(align_channel) # 0=ACTH, 1=Cortisol
        self.last_simulated_values = None
        self.current_parameters = None
        self.current_error = 0

    def n_parameters(self):
        return self.model.n_parameters()

    def _best_circular_shift(self, sim, obs):
        # Demean and scale to comparable magnitude before correlation
        sim = np.asarray(sim, dtype=float)
        obs = np.asarray(obs, dtype=float)
        s = sim - sim.mean()
        o = obs - obs.mean()
        s /= (np.max(np.abs(s)) + 1e-12)
        o /= (np.max(np.abs(o)) + 1e-12)

        n = len(s)
        # Circular cross-correlation via FFT: r[k] = sum_n s[n] * o[n+k]
        corr = np.fft.irfft(np.fft.rfft(s) * np.conj(np.fft.rfft(o)), n=n)
        k = int(np.argmax(corr))
        if k > n // 2:
            k -= n  # convert to signed lag in (-n/2, n/2]
        return k

    def _best_circular_shift_both(self, sim2d, obs2d):
        """Find a single circular shift that maximises combined correlation across channels."""
        sim2d = np.asarray(sim2d, dtype=float)
        obs2d = np.asarray(obs2d, dtype=float)
        n = sim2d.shape[0]
        n_ch = min(sim2d.shape[1], obs2d.shape[1])
        total_corr = None
        for i in range(n_ch):
            s = sim2d[:, i]
            o = obs2d[:, i]
            s = s - s.mean()
            o = o - o.mean()
            s /= (np.max(np.abs(s)) + 1e-12)
            o /= (np.max(np.abs(o)) + 1e-12)
            corr_i = np.fft.irfft(np.fft.rfft(s) * np.conj(np.fft.rfft(o)), n=n)
            total_corr = corr_i if total_corr is None else (total_corr + corr_i)
        k = int(np.argmax(total_corr))
        if k > n // 2:
            k -= n
        return k

    def _normalise_array(self, simulated, observed):
        max_sim = float(np.max(simulated)) if len(simulated) else 0.0
        max_obs = float(np.max(observed)) if len(observed) else 0.0
        max_value = max(max_sim, max_obs)
        if max_value == 0.0:
            return np.zeros_like(simulated, dtype=float), np.zeros_like(observed, dtype=float)
        return simulated / max_value, observed / max_value

    def __call__(self, parameters):
        self.current_parameters = parameters
        sim = self.model.simulate(parameters, self.times)
        self.last_simulated_values = sim
        if sim.ndim == 2 and self.values.ndim == 2 and sim.shape[1] >= 2 and self.values.shape[1] >= 2:
            k = self._best_circular_shift_both(sim, self.values)
        else:
            k = self._best_circular_shift(sim[:, self.align_channel], self.values[:, self.align_channel])
        sim_shifted = np.roll(sim, -k, axis=0)
        normalised_results = [
            self._normalise_array(sim_shifted[:, i], self.values[:, i])
            for i in range(sim_shifted.shape[1])
        ]
        normalised_sim = np.array([normalised_results[0][0], normalised_results[1][0]])
        normalised_obs = np.array([normalised_results[0][1], normalised_results[1][1]])
        error = float(np.sum((normalised_sim - normalised_obs) ** 2))

        self.current_error = error
        return error



class EnvelopeError(pints.ErrorMeasure):
    def __init__(self, model, times, values, window_size=60, polyorder=3):
        self.model = model
        self.times = times
        self.values = values
        self.window_size = window_size
        self.polyorder = polyorder
        self.last_simulated_values = None
        self.current_error = 0
        self.current_parameters = None

    def n_parameters(self):
        return self.model.n_parameters()

    def smooth(self, signal):
        return scipy.signal.savgol_filter(signal, self.window_size, self.polyorder)

    def __call__(self, parameters):
        sim_values = self.model.simulate(parameters, self.times)
        self.last_simulated_values = sim_values

        error = 0
        for i in range(self.values.shape[1]):
            sim = sim_values[:, i]
            obs = self.values[:, i]

            max_val = max(np.max(sim), np.max(obs))
            sim = sim.astype(float) / max_val
            obs = obs.astype(float) / max_val

            sim_smooth = self.smooth(sim)
            obs_smooth = self.smooth(obs)

            error += np.mean((sim_smooth - obs_smooth) ** 2)

        self.current_error = error
        self.current_parameters = parameters
        return error

class PhaseError(pints.ErrorMeasure):
    def __init__(self, model, times, values):
        self.model = model
        self.times = times
        self.values = values
        self.last_simulated_values = None
        self.current_error = 0
        self.current_parameters = None


    def n_parameters(self):
        return self.model.n_parameters()

    def compute_phase(self, signal):
        analytic = hilbert(signal)
        return np.angle(analytic)

    def __call__(self, parameters):
        sim_values = self.model.simulate(parameters, self.times)
        self.last_simulated_values = sim_values

        error = 0
        for i in range(self.values.shape[1]):
            sim = sim_values[:, i]
            obs = self.values[:, i]

            max_val = max(np.max(sim), np.max(obs))
            sim = sim.astype(float) / max_val
            obs = obs.astype(float) / max_val

            phase_sim = self.compute_phase(sim)
            phase_obs = self.compute_phase(obs)

            phase_diff = np.angle(np.exp(1j * (phase_sim - phase_obs)))
            error += np.mean(1 - np.cos(phase_diff))

        self.current_error = error
        self.current_parameters = parameters

        return error
    

