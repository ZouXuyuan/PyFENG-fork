# <center>European Option Pricing Based on Heston-Dupire</center>
### <center>Group 8: Gong Jiaxin, Zou Xuyuan</center>
# 1. Theory
## 1.1 Heston-Dupire Model
<font size=2 >The Heston-Dupire immediate local volatility model is a model that combines the local volatility obtained from the Dupire formula in its unparameterized form with the Heston stochastic volatility model.
</font>

![本地路径](.\picture\2.2small.png)

## 1.2 Leverage function
### 1.2.1 Derivation of leverage function
<font size=2> Let the price of the European call option be　</font>


![本地路径](.\picture\2.6small.png)



<font size=2> Differentiating the the above equation, and using Fubini's theorem, we get, </font>

![本地路径](.\picture\2.7small.png)

<font size=2> This function is not differentiable at x =  c  and cannot be solved directly by Itoˆ 's lemma. However, this problem can be solved by the following Tanaka-Meyer formula. </font>

![本地路径](.\picture\2.8small.png)

<font size=2> To further simplify the calculation, the well-known conclusion from Feng (2010) can be used.  If S<sub>t</sub>, V<sub>t</sub>  obey the stochastic local volatility model in 1.1, then the following equation holds for the price of a European call option.</font>

![本地路径](.\picture\2.10small.png)

<font size=2> Then we get the leverage function</font>

![本地路径](.\picture\2.12small.png)

### 1.2.2 Calculation of leverage function
### Numerator
![本地路径](.\picture\nu.png)
### Denominator


![本地路径](.\picture\de.png)

<font size=2> To get σ<sub>I</sub>  We first use the SVI (Stochastic Volatility Inspired) function to fit in the K-direction , followed by linear interpolation in the T-direction. And then we obtain the implied volatility surface by combining the two directions.
</font>


![本地路径](.\picture\svi.png)

![本地路径](.\picture\minmize.png)

![本地路径](.\picture\linear.png)

<font size=2> To complete the calculation of the rest of leverage function, we take the numerical derivatives by difference.
</font>

![本地路径](.\picture\szqd.png)


# 2. Practice
## 2.1 Fitting volatility based on SVI
<font size=2> We select the CSI 300 put option data on December 29, 2023, and fit the implied volatility surface using the SVI method.
</font>


![本地路径](.\picture\subplot.jpg)

![本地路径](.\picture\surface.jpg)

## 2.2 Model result
* <font size=2> strike = np.arange(2900, 3500, 50) </font>
* <font size=2>sigma, vov, mr, rho, texp, spot = 0.3, 1, 0.5, -0.9, 20, 3431.1099 </font>
### 2.2.1 Heston model

* <font size=2> Lewis AL (2000) Option valuation under stochastic volatility: with Mathematica code. Finance Press</font>
![本地路径](.\picture\conditional_mc.png)

* <font size=2> Conditional MC for Heston model based on QE discretization scheme by Andersen (2008)</font>

![本地路径](.\picture\fft.png)

### 2.2.2 Heston-Dupire
* <font size=2> Var calculation is based on QE method, and we use Euler method to obtain ST </font>

![本地路径](.\picture\eula.png)

![本地路径](.\picture\res.png)