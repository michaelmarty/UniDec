import numpy as np

# --- Calibration data (replace with yours) ---
cal_RV = np.array([13.81, 16.82, 18.09, 19.05, 20.79])    # retention volumes for standards
cal_D  = np.array([17, 10.52, 5.62, 3.5, 1.7])    # known Stokes diameters for standards

# Fit linear calibration D = m*RV + b and get covariance of coeffs
coeffs, cov = np.polyfit(cal_RV, cal_D, 1, cov=True)
m, b = coeffs
Var_m = cov[0,0]
Var_b = cov[1,1]
Cov_mb = cov[0,1]

# --- Your sample triplicate results (replace with yours) ---
RV_means = np.array([16.99, 15.28, 15.89, 16.15, 15.74, 16.18])   # mean retention volume for each sample (triplicate mean)
RV_sds   = np.array([0.09, 0.04, 0.06, 0.12, 0.04, 0.01])   # SD of the triplicate retention volumes
# If you prefer SEM for error-of-mean, do: RV_sems = RV_sds / np.sqrt(3)

# Compute predicted diameters
D_means = m * RV_means + b

# Compute variance from fit uncertainty for each sample RV
Var_fit = (RV_means**2) * Var_m + Var_b + 2 * RV_means * Cov_mb

# Compute variance due to measurement scatter (choose SD or SEM)
Var_meas = (m**2) * (RV_sds**2)        # use RV_sds for SD-based bars
# If you want SEM-based bars (uncertainty of mean), uncomment:
# Var_meas = (m**2) * (RV_sds**2 / 3.0)

Var_total = Var_fit + Var_meas
D_err = np.sqrt(Var_total)

# Now plot
import matplotlib.pyplot as plt
plt.errorbar(D_means, RV_means, xerr=D_err, fmt='o')  # xerr in Stokes diameter
plt.xlabel('Stokes diameter')
plt.ylabel('Retention volume')
plt.gca().invert_yaxis()   # if your RV axis is typically inverted for SEC
plt.show()
