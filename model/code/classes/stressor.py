import numpy as np
import matplotlib.pyplot as plt

class Stressor:
    def __init__(self, stressor_params=None, decay_shape=0.01):
        defaults = {
            "magnitude": 0,
            "start": 0,
            "duration": 0,         # used as rise duration for HS
            "rho": None,           # rise steepness (if None -> auto-compute)
            "beta": 0.8,           # decay-shape exponent
            "plateau": 0.0,        # minutes to hold at M after the rise
            
        }
        self.stressor_params = {**defaults, **(stressor_params or {})}

        self.magnitude = self.stressor_params["magnitude"]
        self.start = self.stressor_params["start"]
        self.duration = self.stressor_params["duration"]
        self.increase_proportion = 0.1
        self.time_to_peak = 10
        # kappa (decay rate): pull from params if provided, else use constructor default
        self.decay_shape = float(self.stressor_params.get("decay_shape", decay_shape))  # kappa

        self.rho = self.stressor_params.get("rho")   # if None, will be set to default
        self.beta = float(self.stressor_params.get("beta", 0.8))
        self.plateau = float(self.stressor_params.get("plateau", 0.0))
        if self.rho is None:
            # Default to a moderately steep rise; caller can override via 'rho'
            self.rho = 4.6

    def stressor(self, t):
        if t < self.start:
            return 0

        peak_time = self.start + self.time_to_peak
        end_time = peak_time + self.duration

        if self.start <= t < peak_time:
            return self.magnitude * (t - self.start) / self.time_to_peak

        if peak_time <= t < end_time:
            decay = -np.log(0.01) / self.duration  
            time_since_peak = t - peak_time
            return self.magnitude * np.exp(-decay * time_since_peak)

        return 0

    def stressor_heart_surgery(self, t):
        """CABG/heart surgery stressor following the paper equation:

        - Rise: normalized exponential approach to magnitude over 'duration'
                f(t) = M * (1 - exp(-rho * s)) / (1 - exp(-rho)), with s=(t-start)/duration
        - Plateau: hold at M for 'plateau' minutes
        - Decay: stretched exponential with parameters kappa=decay_shape and beta
        """
        if t < self.start:
            return 0.0

        # RISE: Exponential approach to M over 'duration', normalized to reach M at rise_end
        rise_end = self.start + self.duration
        if self.start <= t < rise_end:
            s = (t - self.start) / max(self.duration, 1e-12)
            if self.rho <= 1e-8:
                # Near-zero rho -> linear ramp fallback
                return self.magnitude * s
            # Normalized exponential: M * (1 - e^{-rho * s}) / (1 - e^{-rho})
            num = 1.0 - np.exp(-self.rho * s)
            den = 1.0 - np.exp(-self.rho)
            return self.magnitude * (num / den)

        # PLATEAU: Hold at magnitude for 'plateau' minutes
        plateau_end = rise_end + max(self.plateau, 0.0)
        if rise_end <= t < plateau_end:
            return self.magnitude

        # DECAY: Exponential decay from magnitude
        if t >= plateau_end:
            time_since_plateau = t - plateau_end
            kappa = max(float(self.decay_shape), 0.0)  
            beta = max(self.beta, 1e-6)              
            return self.magnitude * np.exp(-kappa * (time_since_plateau ** beta))

        return 0.0