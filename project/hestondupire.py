import numpy as np
from tqdm import tqdm

def generate_hestondupire_paths(S, T, r, kappa, theta, rho, vov, 
                          steps, Npaths, return_vol=False, v_0  = None):
    v_0 = theta if v_0 is None else v_0
    dt = T/steps
    size = (Npaths, steps)
    prices = np.zeros(size)
    sigs = np.zeros(size)
    # S_t = S * np.ones((Npaths,))
    X_t = np.log(S) * np.ones((Npaths,))
    v_t = v_0
    E = v_0 * np.ones((Npaths,))
    L = np.zeros((Npaths, ))
    meanL = []
    bin = 20
    n = Npaths // bin
    TT = [20, 40, 60, 120,180, 250]
    p = np.load(r".\para_svi.npy")
    # 设置随机数种子
    # np.random.seed(20) 
    # for t in tqdm(range(1, steps)):
    for t in tqdm(range(1, steps)):
        WT1 = np.random.normal(size=Npaths) * np.sqrt(dt)
        Z = np.random.normal(size=Npaths) * np.sqrt(dt)
        WT2 = rho*WT1 + np.sqrt(1-rho**2)*Z
        for j in range(Npaths):
            L[j] = LocVolCalibrator((t-1)*dt, X_t[j], p, TT, S, r) ** 2 / E[j]

        # S_t = S_t*(np.exp( (r- 0.5*v_t* L)*dt+ np.sqrt(v_t* L) *WT1 ) ) 
        X_t = X_t + (r - 0.5*v_t* L) * dt + np.sqrt(v_t* L) * WT2
        
        # v_t = np.abs(v_t + kappa*(theta-v_t)*dt + vov*np.sqrt(v_t)*WT2)
        v_t = v_t + kappa*(theta-v_t)*dt + vov*np.sqrt(v_t)*WT1 + 1/4 * vov**2 * (WT1**2 - dt)
        v_t[v_t < 0] = 0
        # print(np.mean(v_t), np.std(v_t))
        prices[:, t] = np.exp(X_t)
        sigs[:, t] = v_t

        index = np.argsort(X_t)
        for k in range(bin):
            ee = np.sum(v_t[index[k * n : (k +1)* n]]) #分段求和
            E[index[k * n : (k +1) * n]] = bin / Npaths * ee
   
    if return_vol:
        return prices, sigs
    
    return prices

def LocVolCalibrator(t, s, p, TT, S0, r):
        lT = len(TT)
        if abs(t - 1) < 10**(-6):
            t = 1 - 10**(-5)
        t_real = t

        x = s-np.log(S0)
        if S0==0:
            print('!!')

        if t_real < TT[0]:
            i = 0
            taon = 1
            taon_1 = 0
        elif t_real >= TT[lT - 1]:
            i = lT - 2
            taon = 0
            taon_1 = 1
        else:
            for i in range(lT - 1):
                if t_real >= TT[i] and t_real < TT[i + 1]:
                    break
            taon = (TT[i + 1] - t_real) / (TT[i + 1] - TT[i])
            taon_1 = (t_real - TT[i]) / (TT[i + 1] - TT[i])

        sign = np.sqrt(max(p[i, 0] + p[i, 1] * (p[i, 2] * (x - p[i, 3]) + np.sqrt(max((x - p[i, 3])**2 + p[i, 4],0))), 0)) / np.sqrt(TT[i] / 250)
        sign_1 = np.sqrt(max(p[i + 1, 0] + p[i + 1, 1] * (p[i + 1, 2] * (x - p[i + 1, 3]) + np.sqrt(max((x - p[i + 1, 3])**2 + p[i + 1, 4], 0))), 0)) / np.sqrt(TT[i + 1] / 250)
        #SVIVol = np.sqrt(max(p[0] + p[1] * (p[2] * (k[j] - p[3]) + np.sqrt(max((k[j] - p[3]) ** 2 + p[4], 0))), 0)) / np.sqrt(TT[i] / 365)  # max is to avoid error
        sign = max(sign, 1e-6) #!!!!!!!!!!!!!
        sign_1 = max(sign_1, 1e-6) #!!!!!!!!!!!!!!!!1
        sig = taon * sign + taon_1 * sign_1

        sig_T = (sign_1 - sign) / (TT[i + 1] - TT[i])

        fn = (p[i + 1, 1] / (2 * sign)) * (p[i, 2] + (x - p[i, 3]) / np.sqrt(max((x - p[i, 3])**2 + p[i, 4], 1e-6))) #!!!!!!!!!!!!!!!!!!!
        fn_1 = (p[i + 1, 1] / (2 * sign_1)) * (p[i + 1, 2] + (x - p[i + 1, 3]) / np.sqrt(max((x - p[i + 1, 3])**2 + p[i + 1, 4], 1e-6))) #!!!!!!!!!!!!!!!!!!!!!!!!
        gn = (p[i, 1] * p[i, 4] / (2 * ((x - p[i, 3])**2 + p[i, 4])**(3 / 2)) - fn**2) / sign
        gn_1 = (p[i + 1, 1] * p[i + 1, 4] / (2 * ((x - p[i + 1, 3])**2 + p[i + 1, 4])**(3 / 2)) - fn_1**2) / sign_1
        sig_x = taon * fn + taon_1 * fn_1
        sig_xx = taon * gn + taon_1 * gn_1

        d1 = ((r + (sig**2) / 2) * t_real / 250 - x) 
        SigmaLV = np.sqrt(max((sig**2 + 2 * t_real / 250 * sig * (sig_T + r * sig_x)) / ((1 + d1/sig*sig_x)**2 + sig * t_real / 250 * (sig_xx - sig_x - d1/sig  * (sig_x)**2)), 10**(-6)))
        if not np.isfinite(SigmaLV):
            print('!!!!!!!!!!')
        return SigmaLV



# strike = np.arange(2900, 3500, 50)
# sigma, vov, mr, rho, texp, spot = 0.5, 1, 0.5, -0.9, 20, 3431.1099
# r = 0.0
# paths =100000
# steps = 160

# prices_neg  = generate_hestondupire_paths(spot, texp, r, mr, sigma,
#                                      rho=rho, vov=vov, steps=steps, Npaths=paths,
#                                     return_vol=False)[:,-1]  

# calls = [] 
# for K in strike:
#     P = np.mean(np.maximum(prices_neg - K,0))*np.exp(-r*texp)
#     calls.append(P)
# print(calls)