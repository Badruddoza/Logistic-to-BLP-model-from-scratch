import statsmodels.tsa.api as smtsa
# Number of samples
n=600
# Generate AR(1) dataset
ar = np.r_[1, -0.6]
ma = np.r_[1, 0]
ar1_data = smtsa.arma_generate_sample(ar=ar, ma=ma, nsample=n)
plotds(ar1_data)
