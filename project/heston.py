import numpy as np

def generate_heston_paths(S, T, r, kappa, theta, rho, vov, 
                          steps, Npaths, return_vol=False, v_0  = None):
    v_0 = theta if v_0 is None else v_0
    dt = T/steps
    size = (Npaths, steps)
    prices = np.zeros(size)
    sigs = np.zeros(size)
    # S_t = S
    X_t = np.log(S) * np.ones((Npaths,))
    v_t = v_0
    # 设置随机数种子
    # np.random.seed(20) 
    for t in range(steps):
        WT1 = np.random.normal(size=Npaths) * np.sqrt(dt)
        Z = np.random.normal(size=Npaths) * np.sqrt(dt)
        WT2 = rho*WT1 + np.sqrt(1-rho**2)*Z
        # S_t = S_t*(np.exp( (r- 0.5*v_t)*dt+ np.sqrt(v_t) *WT1 ) ) 
        X_t = X_t + (r - 0.5*v_t) * dt + np.sqrt(v_t) * WT1
        S_t = np.exp(X_t)
        # v_t = np.abs(v_t + kappa*(theta-v_t)*dt + vov*np.sqrt(v_t)*WT2)
        v_t = v_t + kappa*(theta-v_t)*dt + vov*np.sqrt(v_t)*WT2 + 1/4 * vov**2 * (WT2**2 - dt)
        v_t[v_t < 0] = 0
        prices[:, t] = S_t
        sigs[:, t] = v_t
    
    if return_vol:
        return prices, sigs
    
    return prices

