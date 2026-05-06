import pints
from ddeint import ddeint
import math
import numpy as np
from model.code.classes.stressor import Stressor
from model.code.config.parameter_boundaries import PARAMETER_BOUNDARIES

#Define pints model class
class Model(pints.ForwardModel):
    def __init__(self
                 , fixed_params
                 , suggested_params
                 , num_days = 4
                 , signal ='Both'
                 , signal_range = (7,13)
                 , reject=True
                 , stressor_parameters = {} 
                 , length_model =1440
                 , stressor_type = ''
                 , days_to_keep = 1
        ):
        self.signal = signal
        self.num_days = num_days
        self.days_to_keep = days_to_keep
        self.fixed_params = fixed_params
        self.suggested_params_dict = suggested_params
        self.n_parameters_value = len(suggested_params) 
        self.signal_range = signal_range
        self.parameters = {**fixed_params, **suggested_params}
        self.reject = reject
        self.length_model = length_model
        self.stressor_type = stressor_type
        
        self.parameter_boundaries = PARAMETER_BOUNDARIES.copy()
        
        if stressor_parameters != {}: 
            self.stressor_object = Stressor(stressor_parameters)
        else: 
            self.stressor_object = None 

    def crh(self, t, t_s=None, lambda_a=None, lambda_s=None, sigma=None, T_c=1440):
            if t_s is None: 
                t_s = self.parameters['t_s']
            if lambda_a is None: 
                lambda_a = self.parameters['lambda_a']
            if lambda_s is None: 
                lambda_s = self.parameters['lambda_s']
            if sigma is None: 
                sigma = self.parameters['sigma']

            crh  =  lambda_a * math.e**(lambda_s*math.cos(2*math.pi*((t-t_s)/T_c)+sigma*math.cos(2*math.pi*((t-t_s)/T_c))))

            if self.stressor_object != None:
                if self.stressor_type == 'HS': 
                    stressor_value = self.stressor_object.stressor_heart_surgery(t)
                else: 
                    stressor_value = self.stressor_object.stressor(t)
                crh += stressor_value
                
            return crh

    def simulate(self, parameters, times):
        param_keys = list(self.suggested_params_dict.keys())
        for i, key in enumerate(param_keys):
            self.parameters[key] = parameters[i]

        h = self.parameters['h']
        k_a = self.parameters['k_a']
        k_c = self.parameters['k_c']
        alpha = self.parameters['alpha']
        delay = self.parameters['delay']
        t_s = self.parameters['t_s']
        gamma_a = self.parameters['gamma_a']
        gamma_c = self.parameters['gamma_c']
        sigma = self.parameters['sigma']
        m_a = self.parameters['m_a']
        m_c = self.parameters['m_c']
        lambda_a =self.parameters['lambda_a']
        lambda_s =self.parameters['lambda_s']

        def model(Y, t):
            A, C= Y(t)
            C_delay = Y(t - delay)[1]

            dAdt = -gamma_a*A + h*((k_c**m_a)*self.crh(t, t_s, lambda_a, lambda_s, sigma))/(k_c**m_a+C_delay**m_a)
            dCdt = -gamma_c*C + alpha*((A**m_c)/(k_a**m_c + A**m_c))

            return [dAdt, dCdt]

        def initial_conditions(t):
            return [5, 400]  
        
        result = ddeint(model, initial_conditions, times)

        result = result[self.length_model*(self.num_days-self.days_to_keep):]

        if self.signal == 'Cortisol': 
            result = result[:, 1]

        if self.signal == 'ACTH': 
            result = result[:, 0]

        if result.ndim == 1:
            result = result[:, np.newaxis]

        if (self.reject == True):
            if (self.reject_parameter_combination(result)):
                return np.full((self.length_model, result.ndim), 5000)

        return result
    
    def n_outputs(self):
        return 2

    def n_times(self): 
        return self.length_model*self.num_days

    def n_parameters(self):
        return self.n_parameters_value

    def suggested_parameters(self): 
        return list(self.suggested_params_dict.values())
    

    def get_and_create_boundaries(self):
        lowerbounds = []
        upperbounds = []

        for item in self.suggested_params_dict.keys():
            lowerbound, upperbound = self.parameter_boundaries.get(item, (None, None))
            if lowerbound is None or upperbound is None: 
                print(f'{item} has no defined bounds. Setting to default (0, 1000)')
                lowerbound, upperbound = (0, 1000)
            
            lowerbounds.append(lowerbound)
            upperbounds.append(upperbound)
        

        return pints.RectangularBoundaries(lowerbounds, upperbounds)


    # Function to reject parameter combination if number of peaks are outside a plausible range
    def reject_parameter_combination(self, result): 
        from scipy import signal as scipy_signal
        
        for i in range(result.shape[1]):
            signals, _ = scipy_signal.find_peaks(result[:, i])
            number_of_signals = len(signals)

            lower_bound, upper_bound = self.signal_range

            if not (lower_bound <= number_of_signals <= upper_bound): 
                print(f'Number of ultradian pulses are unrealistic')
                print(f'Number of signal peaks: {number_of_signals}')
                print(f'We are expecting it to be between {lower_bound} and {upper_bound}')
                print(f'Rejecting parameter combination')
                return True
        
        return False



class ModelConstCRH(pints.ForwardModel):
    def __init__(self, fixed_params, suggested_params, num_days = 4, signal ='Both', signal_range = (7,13), days_to_keep = 1, length_model = 1440):
        self.signal = signal
        self.num_days = num_days
        self.days_to_keep = days_to_keep
        self.fixed_params = fixed_params
        self.suggested_params_dict = suggested_params
        self.n_parameters_value = len(suggested_params) 
        self.signal_range = signal_range
        self.length_model = length_model
        self.parameters = {**fixed_params, **suggested_params}

    def simulate(self, parameters, times):

        param_keys = list(self.suggested_params_dict.keys())
        for i, key in enumerate(param_keys):
            self.parameters[key] = parameters[i]

        h = self.parameters['h']
        k_c = self.parameters['k_c']
        k_a = self.parameters['k_a']
        alpha = self.parameters['alpha']
        delay = self.parameters['delay']
        t_s = self.parameters['t_s']
        gamma_a = self.parameters['gamma_a']
        gamma_c = self.parameters['gamma_c']
        sigma = self.parameters['sigma']
        m_a = self.parameters['m_a']
        m_c = self.parameters['m_c']
        lambda_a =self.parameters['lambda_a']
        lambda_s =self.parameters['lambda_s']

        def model(Y, t):
            A, C= Y(t)
            
            A_delay = Y(t - delay)[0]
            #C_delay1 = Y(t - delay)[1]

            dAdt = -gamma_a*A + h*((k_c**m_a)*lambda_a)/(k_c**m_a+C**m_a)
            dCdt = -gamma_c*C + alpha*((A_delay**m_c)/(k_a**m_c + A_delay**m_c))

            return [dAdt, dCdt]

        def initial_conditions(t):
            return [5, 400]  
        
        result = ddeint(model, initial_conditions, times)

        result = result[self.length_model*(self.num_days-self.days_to_keep):]

        if self.signal == 'Cortisol': 
            result = result[:, 1]

        if self.signal == 'ACTH': 
            result = result[:, 0]

        if result.ndim == 1:
            result = result[:, np.newaxis]

        return result
    
    def n_outputs(self):
        return 2

    def n_times(self): 
        return self.length_model*self.num_days

    def n_parameters(self):
        return self.n_parameters_value

    def suggested_parameters(self): 
        return list(self.suggested_params_dict.values())

