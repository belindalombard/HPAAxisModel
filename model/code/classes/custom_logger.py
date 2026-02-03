import os
import matplotlib.pyplot as plt
import json



class CustomLogger:
    """Logger for optimization progress, errors, and plots."""
    def __init__(self, error_measure, num_days, signal, output_dir=None, save_plot=True, participant_number=None, error_penalty=None, plot_on_improvement=True):
        self.iteration = 0
        self.participant_number = participant_number
        self.save_plot = save_plot
        self.error_measure = error_measure
        self.error_penalty = error_penalty
        self.num_days = num_days
        self.signal = signal
        self.output_dir = output_dir
        self.current_best_fit = float('inf')
        self.plot_on_improvement = plot_on_improvement

        if self.output_dir:
            abs_output_dir = os.path.abspath(self.output_dir)
            os.makedirs(abs_output_dir, exist_ok=True)
            self.config_save_file = os.path.join(self.output_dir, "best_fit_config.json")
            self.output_file = os.path.join(self.output_dir, "iterations.output")
            self.best_fit_file = os.path.join(self.output_dir, "current_best_fit")
            with open(self.output_file, 'w') as f:
                f.write("")

    def __call__(self, iteration, optimiser):
        """Log progress and errors at each optimizer iteration."""
        self.iteration += 1
        xk, fxk = None, None
        try:
            xk = optimiser.xbest()
            fxk = optimiser.fbest()
        except AttributeError:
            pass

        penalty_val = None
        if self.error_penalty is not None:
            # Do not recompute penalty; reuse last computed value from the penalty error
            penalty_val = getattr(self.error_penalty, 'current_error', None)

        self._write_best_fit_file(penalty_val)
        self._write_output_file(xk, fxk, penalty_val)

        improved = self.error_measure.current_error < self.current_best_fit
        if improved:
            self.current_best_fit = self.error_measure.current_error
            if self.output_dir:
                self._save_config_file(penalty_val)

        # Plot only on improvement (if enabled) or every iteration when disabled
        if self.output_dir and self.save_plot and (not self.plot_on_improvement or improved):
            if self.signal == 'Both':
                self.plot_simulation_vs_actual_both(self.error_measure.last_simulated_values)
            else:
                self.plot_simulation_vs_actual_single(self.error_measure.last_simulated_values)

    def _write_best_fit_file(self, penalty_val):
        if not self.output_dir:
            return
        with open(self.best_fit_file, 'w') as f:
            f.write(f"Iteration {self.iteration}\n")
            f.write(f"Parameters: {self.error_measure.current_parameters}\n")
            f.write(f"Error: {self.error_measure.current_error}\n")
            if penalty_val is not None:
                f.write(f"Penalty error: {penalty_val}\n")

    def _write_output_file(self, xk, fxk, penalty_val):
        if not self.output_dir:
            return
        with open(self.output_file, 'a') as f:
            f.write(f"Iteration {self.iteration}\n")
            f.write(f"Optimiser best: {xk}, Error: {fxk}\n")
            f.write(f"ErrorMeasure best: {self.error_measure.current_parameters}, Error: {self.error_measure.current_error}\n")
            if penalty_val is not None:
                f.write(f"Penalty error: {penalty_val}\n")
            f.write("\n")


    def plot_simulation_vs_actual_both(self, simulated_values):
        """Plot both signals (ACTH, CORT) vs simulated values."""
        if simulated_values is None:
            print("[CustomLogger] Skipping plot: simulated_values is None")
            return
        fig, ax = plt.subplots(2, figsize=(5, 4))
        t = self.error_measure.times[((self.num_days - 1) * 1440):]
        ax[0].plot(t, self.error_measure.values[:, 0], 'b-', label='Actual ACTH')
        ax[0].plot(t, simulated_values[:, 0], 'b--', label='Simulated ACTH')
        ax[0].set_ylabel('ACTH')
        ax[0].legend()
        ax[1].plot(t, self.error_measure.values[:, 1], 'g-', label='Actual CORT')
        ax[1].plot(t, simulated_values[:, 1], 'g--', label='Simulated CORT')
        ax[1].set_ylabel('CORT')
        ax[1].set_xlabel('Time')
        ax[1].legend()
        fig.suptitle(f'Iteration {self.iteration}')
        fig.tight_layout()
        img_path = os.path.join(self.output_dir, f'Iteration_{self.iteration}.png')
        plt.savefig(img_path, dpi=72)
        plt.close(fig)


    def plot_simulation_vs_actual_single(self, simulated_values):
        """Plot a single signal vs simulated values."""
        if simulated_values is None:
            print("[CustomLogger] Skipping plot: simulated_values is None")
            return
        fig, ax = plt.subplots(1, figsize=(5, 2.5))
        t = self.error_measure.times[((self.num_days - 1) * 1440):]
        ax.plot(t, self.error_measure.values[:, 0], 'b-', label='Actual')
        ax.plot(t, simulated_values[:, 0], 'b--', label='Simulated')
        ax.set_ylabel('Signal')
        ax.set_xlabel('Time')
        ax.legend()
        fig.suptitle(f'Iteration {self.iteration}')
        fig.tight_layout()
        img_path = os.path.join(self.output_dir, f'Iteration_{self.iteration}.png')
        plt.savefig(img_path, dpi=72)
        plt.close(fig)


    def _save_config_file(self, penalty_val=None):
        """Save current and best fit parameters and errors to a JSON file."""
        if self.error_measure.current_parameters is None:
            return
        best_params = self.error_measure.current_parameters if self.error_measure.current_error == self.current_best_fit else getattr(self, 'best_parameters', self.error_measure.current_parameters)
        best_penalty = None
        if self.error_penalty is not None and best_params is not None:
            try:
                best_penalty = self.error_penalty(best_params)
            except Exception as e:
                print(f"[CustomLogger] Penalty error evaluation for best fit failed: {e}")
        config_data = {
            "current_fit": {
                "parameters": dict(zip(
                    self.error_measure.model.suggested_params_dict.keys(),
                    self.error_measure.current_parameters
                )),
                "fixed_params": self.error_measure.model.fixed_params,
                "signal": self.signal,
                "participant": self.participant_number,
                "num_days": self.num_days,
                "error": self.error_measure.current_error
            },
            "best_fit": {
                "parameters": dict(zip(
                    self.error_measure.model.suggested_params_dict.keys(),
                    best_params
                )),
                "fixed_params": self.error_measure.model.fixed_params,
                "signal": self.signal,
                "participant": self.participant_number,
                "num_days": self.num_days,
                "error": self.current_best_fit
            }
        }
        if best_penalty is not None:
            config_data["best_fit"]["penalty_error"] = best_penalty
        if self.error_measure.current_error == self.current_best_fit:
            self.best_parameters = self.error_measure.current_parameters.copy()
        with open(self.config_save_file, 'w') as f:
            json.dump(config_data, f, indent=2)
