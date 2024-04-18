
import numpy as np
import pyfeng as pf

# sigma: model volatility or variance at t=0.
# vov: volatility of volatility
# rho: correlation between price and volatility
# mr: mean-reversion speed (kappa)
# theta: long-term mean of volatility or variance. If None, same as sigma
# intr: interest rate (domestic interest rate)
# divr: dividend/convenience yield (foreign interest rate)
# is_fwd: if True, treat `spot` as forward price. False by default.
strike = np.arange(3200, 3800, 50)
#strike =np.arange(2900, 3500, 50)
sigma, vov, mr, rho, texp, spot = 0.3, 1, 0.5, -0.9, 20, 3431.1099

m = pf.HestonMcAndersen2008(sigma, vov=vov, mr=mr, rho=rho)
m.set_num_params(n_path=1e5, dt=1 / 8, rn_seed=123456)
print('true price:', m.price(strike, spot, texp))

# %%
m = pf.HestonFft(sigma, vov=vov, mr=mr, rho=rho)
print('true price:', m.price(strike, spot, texp))

# %%

# Heaton-Dupire
from heston_mc import HestonMcAndersen2008

m = HestonMcAndersen2008(sigma, vov=vov, mr=mr, rho=rho)
m.set_num_params(n_path=1e5, dt=1 / 8, rn_seed=123456)
m.price(strike, spot, texp)
print('end')

# %%
# %%
