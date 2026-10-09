# Freeze-out conditions

Here, several different freeze-out conditions are defined.

## 1. Comparison of expansion and weak-interaction timescales
Electron fraction $Y_\mathrm{e}$ evolves with

$$\frac{dY_\mathrm{e}}{dt} = -\lambda_\mathrm{ec} Y_p + \lambda_\mathrm{pc} Y_n,$$

where $\lambda_\mathrm{ec}$ and $\lambda_\mathrm{pc}$ are the reaction rates of electron and positron captures, 
and $Y_p$ and $Y_n$ are the number fractions of free protons and neutrons. Here, the weak-interaction processes
of nuclei and the neutrino absorption process are ignored.

We define the weak-interaction timescale as

$$t_\mathrm{weak} = \frac{1}{\lambda_\mathrm{ec} + \lambda_\mathrm{pc}}.$$

Under sufficiently high temperature, we can ignore the presence of nuclei,
and we can assume $Y_p \approx Y_\mathrm{e}$ and $Y_n \approx 1-Y_\mathrm{e}$. 
The evolution equation of $Y_\mathrm{e}$ then becomes 

$$\frac{dY_\mathrm{e}}{dt} = -(\lambda_\mathrm{ec} + \lambda_\mathrm{pc})(Y_\mathrm{e} - Y_\mathrm{e}^\mathrm{(eq)}) = -\frac{Y_\mathrm{e} - Y_\mathrm{e}^\mathrm{(eq)}}{t_\mathrm{weak}},$$

where we define the equilibrium electron fraction as

$$ Y_\mathrm{e}^\mathrm{(eq)} = \frac{\lambda_\mathrm{pc}}{\lambda_\mathrm{ec} + \lambda_\mathrm{pc}}.$$

Thus, $t_\mathrm{weak}$ can be regarded as the timescale on which the electron fraction approaches the equilibrium value.

We define the expansion time as

$$t_\mathrm{exp} = r/v^r$$.

Under the above definitions, the first definition for the freeze-out time is the latest time along the Lagrangian evolution at which $t_\mathrm{exp} = t_\mathrm{weak}$.
This defines the freeze-out locally.

## 2. Net change of $Y_\mathrm{e}$
We can also define the freeze-out using the history of the weak-interaction reactions that each Lagrangian particle undergoes.
For example, we define the freeze-out time as

$$\Delta Y_\mathrm{e}|_\mathrm{thr} = \int_t^{t_\mathrm{fin}} \frac{dY_\mathrm{e}}{dt}(\rho(t'), T(t'), Y_\mathrm{e}(t'))dt.$$

Here, we note that all terms in the right-hand side of $dY_\mathrm{e}/dt$ can be evaluated for given $(\rho, T, Y_\mathrm{e})$ under the assumption of NSE.
$\Delta Y_\mathrm{e}|_\mathrm{thr}$ is the arbitral threshold value for the definition of the freeze-out. Default value is 0.01.

## 3. Total variation
The criterion with $\Delta Y_\mathrm{e}|_\mathrm{thr}$ misdetermines the freeze-out when the increase and decrease in $Y_\mathrm{e}$ cancel.
To avoid this, we can define the second criterion as

$$\Delta N|_\mathrm{var,thr} = \int_t^{t_\mathrm{fin}} \big|-\lambda_\mathrm{ec} Y_p + \lambda_\mathrm{pc} Y_n\big|dt.$$

The right-hand side in the above expression is the total number of $Y_\mathrm{e}$ change after the time $t$.

## 4. Gross reaction number
The last and the most stringent criterion can be defined as

$$\Delta N|_\mathrm{gr,thr} = \int_t^{t_\mathrm{fin}} (\lambda_\mathrm{ec} Y_p + \lambda_\mathrm{pc} Y_n)dt.$$

The right-hand side in the above expression is the number of total weak-interaction processes after the time $t$.

Unlike the total-variation criterion, this quantity remains finite even when
electron and positron captures nearly balance each other and therefore produce
little net change in $Y_\mathrm{e}$.

Thus,

- the net-change criterion measures the remaining net change in $Y_\mathrm{e}$;
- the total-variation criterion measures the remaining variation of $Y_\mathrm{e}$ without temporal cancellation;
- the gross-reaction criterion measures the remaining weak-interaction activity itself.
  
