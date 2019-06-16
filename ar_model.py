import statsmodels.tsa.api as smtsa
import numpy as np
# Number of samples
n=600
# Generate AR(1) dataset
ar = np.r_[1, -0.6]
ma = np.r_[1, 0]
ar1_data = smtsa.arma_generate_sample(ar=ar, ma=ma, nsample=n)
# Plot the AR1 process
plotds(ar1_data)
