

PARAMETER_BOUNDARIES = {
    'k_a': (1, 25),
    'k_c': (8, 150),
    'alpha': (1, 100),
    'delay': (1, 40),
    'lambda_s': (0, 20),
    'lambda_a': (0, 20),
    't_s': (-1440, 1440),
    'm_a': (2, 5),
    'm_c': (2, 5),
    'sigma': (0, 1),
    'gamma_a': (0.0231, 0.0693),
    'gamma_c': (0.00578, 0.00991),
    # Sleep model specific
    'mu_wake': (0.0003, 0.003),
    'mu_sleep': (0.002, 0.010),
    'phase_shift': (300, 700),
    'omega': (0, 0.01),
    'hysteresis': (0.1, 0.4),
    'A_c': (0, 2),
    'crh_reduction': (-1, 1)
}

