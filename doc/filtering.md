
# Filtering and Smoothing

## Filtering

From [@DurbinKoopman2012, Sec 4.3]

Let $\vec{a}_t = \E(\vec{\alpha}_t | y_{1, \dots, t - 1})$ be the expected value
$\vec{P}_t = \Var(\vec{\alpha}_t | y_{1, \dots, t - 1})$ be the variance
of the state in $t + 1$ given data up to time $t$.
To calculate $\vec{\alpha}_{t + 1}$ and $\mat{P}_{t + 1}$ given the arrival of new
data at time $t$,
$$
\begin{aligned}[t]
\vec{v}_t &= \vec{y}_t - \mat{Z}_t \vec{a}_t - \vec{d}_t, \\
\mat{F}_t &= \mat{Z}_t \mat{P}_t \mat{Z}_t\T + \mat{H}_t, \\
\mat{K}_t &= \mat{T}_t \mat{P}_t \mat{Z}_t\T \mat{F}_t^{-1} \\
\vec{a}_{t + 1} &= \mat{T}_t \vec{a}_t + \mat{K}_t \vec{v}_t + \vec{c}_t \\
\mat{P}_{t + 1} &= \mat{T}_t \mat{P}_t (\mat{T}_t - \mat{K}_t \mat{Z}_t)\T + \mat{R}_t \mat{Q}_t \mat{R}_t\T .
\end{aligned}
$$
The vector $\vec{v}_t$ are the *one-step ahead forecast errors$,
and the matrix $\mat{K}_t$ is called the *Kalman gain*.

The filter can also be written to estimate the *filtered states*,
where $\vec{a}_{t|t} = \E(\vec{\alpha}_t | y_{1, \dots, t})$ is the expected value
and $\vec{P}_{t|t} = \Var(\vec{\alpha}_t | y_{1, \dots, t})$ is the variance of
the state $\vec{\alpha}_t$ given information up to *and including* $\vec{y}_t$.
The filter written this way is,
$$
\begin{aligned}[t]
\vec{v}_t &= \vec{y}_t - \mat{Z}_t \vec{a}_t - \vec{d}_t, \\
\mat{F}_t &= \mat{Z}_t \mat{P}_t \mat{Z}_t\T + \mat{H}_t, \\
\vec{a}_{t|t} &= \vec{a}_t + \mat{P}_t \mat{Z}_t\T \mat{F}_t^{-1} v_t , \\
\mat{P}_{t|t} &= \mat{P}_t - \mat{P}_t \mat{Z}_t\T \mat{F}_t^{-1} \mat{Z}_t \mat{P}_t , \\
\vec{a}_{t + 1} &= \mat{T}_{t} \vec{a}_{t|t} + \vec{c}_t, \\
\mat{P}_{t + 1} & = \mat{T}_t \mat{P}_{t|t} \mat{T}_t\T + \mat{R}_t \mat{Q}_t \mat{R}_t\T .
\end{aligned}
$$



matrix/vector         dimension
--------------------- --------------
$\vec{v}_t$           $p \times 1$
$\vec{a}_t$           $m \times 1$
$\vec{a}_{t|t}$       $m \times 1$
$\mat{F}_t$           $p \times p$
$\mat{K}_t$           $m \times p$
$\mat{P}_t$           $m \times m$
$\mat{P}_{t|T}$       $m \times m$
$\vec{x}_t$           $m \times 1$
$\mat{L}_t$           $m \times m$
--------------------- --------------

Table: Dimensions of matrices and vectors in the SSM

See [@DurbinKoopman2012, Sec 4.3.4]:
For a time-invariant state space model, the Kalman recursion for $\mat{P}_{t + 1}$ converges to a constant matrix $\bar{\mat{P}}$,
$$
\bar{\mat{P}} = \mat{T} \bar{\mat{P}} \mat{T}\T - \mat{T} \bar{\mat{P}} \mat{Z}\T \bar{\mat{F}}^{-1} \mat{Z} \bar{\mat{P}} \mat{T}\T + \mat{R} \mat{Q} \mat{R}\T ,
$$
where $\bar{\mat{F}} = \mat{Z} \bar{\mat{P}} \mat{Z}\T + \mat{H}$.

See [@DurbinKoopman2012, Sec 4.3.5]:
The *state estimation error* is,
$$
\vec{x}_t = \vec{\alpha}_t - \vec{a}_t,
$$
where $\Var(\vec{x}_t) = \mat{P}_t$.
The $v_t$ are sometimes called *innovations*, since they are the part of $\vec{y}_t$ not predicted from the past.
The innovation analog of the state space model is
$$
\begin{aligned}[t]
\vec{v}_t &= \mat{Z}_t \vec{x}_t + \vec{\varepsilon}_t ,  \\
\vec{x}_{t + 1} &= \mat{L} \vec{x}_{t} + \mat{R}_t \vec{\eta}_t - \mat{K}_t \vec{\varepsilon}_t , \\
\mat{K}_t &= \mat{T}_t \mat{P}_t \mat{Z}_t\T \mat{F}_t^{-1} , \\
\mat{L}_t &= \mat{T}_t - \mat{K}_t \mat{Z}_t ,
\mat{P}_{t + 1} &= \mat{T}_t \mat{P}_t \mat{L}_t\T +  \mat{R}_t \mat{Q}_t \mat{R}_T\T  .
\end{aligned}
$$
These recursions allow for a simpler derivation of $\mat{P}_{t + 1}$, and are useful for
the smoothing recursions.
Moreover, the one-step ahead forecast errors are indendendent, which allows for a simple derivation of the log-likelihood.

Alternative methods **TODO**

- square-root filtering
- precision filters
- sequential filtering

## Smoothing

While filtering calculates the conditional densities the states and disturbances given data prior to or up to the current time,
smoothing calculates the conditional densities states and disturbances given the entire series of observations, $\vec{y}_{1:n}$.

*State smoothing* calculates the conditional mean, $\hat{\vec{\alpha}}_t = \E(\vec{\alpha}_t | \vec{y}_{1:n})$, and variance, $\mat{V}_t = \Var(\vec{\alpha}_t | \vec{y}_{1:n})$, of the states.
Observation disturbance smoothing calculates the conditional mean, $\hat{\vec{\varepsilon}}_t = \E(\vec{\varepsilon}_t | \vec{y}_{1:n})$, and variance, $\Var(\vec{\varepsilon}_t | \vec{y}_{1:n})$, of the state disturbances.
Likewise, state disturbance smoothing calculates the conditional mean, $\hat{\vec{\eta}}_t = \E(\vec{\eta}_t | \vec{y}_{1:n})$, and variance, $\Var(\vec{\eta}_t | \vec{y}_{1:n})$, of the state disturbances.


Vector/Matrix                Dimension
---------------------------- -----------------------
$\vec{r}_t$                   $m \times 1$
$\vec{\vec{\alpha}}_t$        $m \times 1$
$\vec{u}_t$                   $p \times 1$
$\hat{\vec{\varepsilon}}_t$   $p \times 1$
$\hat{\vec{\eta}}_t$          $r \times 1$
$\mat{N}_t$                   $m \times m$
$\mat{V}_t$                   $m \times m$
$\mat{D}_t$                   $p \times p$
---------------------------- -----------------------

Table: Dimensions of vectors and matrices used in smoothing recursions



### State Smoothing

Smoothing calculates conditional density of the states given all observations, $p(\vec{\alpha} | \vec{y}_{1:n})$.  Let $\hat{\vec{\alpha}} = \E(\vec{\alpha}_t | \vec{y}_{1:n})$ be the mean and $\mat{V}_t = \Var(\vec{\alpha} | \vec{y}_{1:n})$ be the variance of this density.
The following recursions can be used to calculate these densities [@DurbinKoopman2012, Sec 4.4.4],
$$
\begin{aligned}[t]
\vec{r}_{t - 1} &= \mat{Z}_t\T \mat{F}_t^{-1} \vec{v}_t + \mat{L}_t\T \vec{r}_t , &
\mat{N}_{t - 1} &= \mat{Z}_t\T \mat{F}_t^{-1} \mat{Z}_t + \mat{L}_t\T \mat{N}_t \mat{L}_t, \\
\hat{\vec{\alpha}}_t &= \vec{a}_t + \mat{P}_t \vec{r}_{t - 1} , &
\mat{V}_t &= \mat{P}_t - \mat{P}_t \mat{N}_{t - 1} \mat{P}_t ,
\end{aligned}
$$
for $t = n, \dots, 1$, with $\vec{r}_n = \vec{0}$, and $\mat{N}_n = \mat{0}$.

During the filtering pass $\vec{v}_t$, $\mat{F}_t$, $\mat{K}_t$, and $\mat{P}_t$ for $t = 1, \dots, n$ need to be stored.
Alternatively, $\vec{a}_t$ and $\mat{P}_t$ only can be stored, and $\vec{v}_t$, $\mat{F}_t$, $\mat{K}_t$ recalculated on the fly.
However, since the dimensions of $\vec{v}_t$, $\mat{F}_t$, $\mat{K}_t$ are usually small relative to $\vec{a}_t$ and $\mat{P}_t$ is is usually worth storing them.



### Disturbance smoothing



Disturbance smoothing calculates the density of the state and observation disturbances ($\vec{\eta}_t$ and $\vec{\varepsilon}_t$) given the full series of observations $\vec{y}_{1:n}$.
Let $\hat{\vec{\varepsilon}}_t = \E(\vec{\varepsilon} | \vec{y}_{1:n})$ be the mean and $\Var(\vec{\varepsilon}_t | \vec{y}_{1:n})$ be the variance of the smoothed density of the observation disturbances at time $t$, $p(\vec{\varepsilon}_t | \vec{y}_{1:n})$.
Likewise, let $\hat{\vec{\eta}} = \E(\vec{\eta}_t | \vec{y}_{1:n})$ be the mean and $\Var(\vec{\eta}_t | \vec{y}_{1:n})$  be the variance of the smoothed density of the state disturbances at time $t$, $p(\vec{\eta}_t | \vec{y}_{1:n})$.
The following recursions can be used to calculate these values [@DurbinKoopman2012, Eq 4.69]:
$$
\begin{aligned}[t]
\hat{\vec{\varepsilon}}_t &= \mat{H}_t (\mat{F}^{-1} \vec{v}_t - \mat{K}_t\T \vec{r}_t) , &
\Var(\vec{\varepsilon}_t | \mat{Y}_n) &= \mat{H}_t - \mat{H}_t (\mat{F}_t^{-1} + \mat{K}_t\T \mat{N}_t \mat{K}_t) \mat{H}_t , \\
\hat{\vec{\eta}}_t &= \mat{Q}_t \mat{R}_t\T \vec{r}_t , &
\Var(\vec{\eta}_t | \mat{Y}_n) &= \mat{Q}_t - \mat{Q}_t \mat{R}_t\T \mat{N}_t \mat{R}_t \mat{Q}_t , \\
\vec{r}_{t - 1} &= \mat{Z}_t\T \mat{F}_t^{-1} \vec{v}_t + \mat{L}_t\T \vec{r}_t , &
\mat{N}_{t - 1} &= \mat{Z}_t\T \mat{F}_t^{-1} \mat{Z}_t + \mat{L}_t\T \mat{N}_t \mat{L}_t
\end{aligned}
$$
Alternatively, these equations can be rewritten as [@DurbinKoopman2012, Sec 4.5.3]:
$$
\begin{aligned}[t]
\hat{\vec{\varepsilon}}_t &= \mat{H}_t \vec{u}_t , &
\Var(\vec{\varepsilon}_t | \mat{Y}_n) &= \mat{H}_t - \mat{H}_t \mat{D}_t \mat{H}_t , \\
\hat{\vec{\eta}}_t &= \mat{Q}_t \mat{R}_t\T \vec{r}_t , &
\Var(\vec{\eta}_t | \mat{Y}_n) &= \mat{Q}_t - \mat{Q}_t \mat{R}_t\T \mat{N}_t \mat{R}_t \mat{Q}_t , \\
\vec{u}_t &= \mat{F}^{-1} \vec{v}_t - \mat{K}_t\T \vec{r}_t , &
\mat{D}_t &= \mat{F}_t^{-1} + \mat{K}_t\T \mat{N}_t \mat{K}_t , \\
\vec{r}_{t - 1} &= \mat{Z}_t\T \vec{u}_t + \mat{T}_t\T \vec{r}_t , &
\mat{N}_{t - 1} &= \mat{Z}_t\T \mat{D}_t \mat{Z}_t + \mat{T}_t\T \mat{N}_t \mat{T}_t - \mat{Z}_t\T \mat{K}_t\T \mat{N}_t \mat{T}_t - \mat{T}_t\T \mat{N}_t \mat{K}_t \mat{Z}_t .
\end{aligned}
$$
This reformulation can be computationally useful since it relies on the system matrices $\mat{Z}_t$ and $\mat{T}_t$ which are often sparse.
The disturbance smoothing recursions require only $\vec{v}_t$, $\mat{f}_t$, and $\mat{K}_t$ which are calculated with a forward pass of the Kalman filter.
Unlike the state smoother, the disturbance smoothers do not require either the mean ($\vec{a}_t$) or variance ($\mat{P}_t$) of the predicted state density.



### Fast state smoothing



If the variances of the states do not need to be calculated, then a faster smoothing algorithm can be used (Koopman 1993).
The fast state smoother is defined as [@DurbinKoopman2012, Sec 4.6.2],
$$
\begin{aligned}[t]
\hat{\vec{\alpha}}_t &= \mat{T}_t \hat{\vec{\alpha}}_t + \mat{R}_t \mat{Q}_t \mat{R}_t\T \vec{r}_t , && t = 2, \dots, n \\
\hat{\vec{\alpha}}_1 &= \vec{a}_1 + \mat{P}_1 \vec{r}_0 .
\end{aligned}
$$
The values of $\vec{r}_t$ come from the recursions in the disturbance smoother.



## Simulation smoothers

Simulation smoothing draws samples of the states, $p(\vec{\alpha}_1, \dots, \vec{\alpha}_n | \vec{y}_{1:n})$, or disturbances, $p(\vec{\varepsilon}_1, \dots, \vec{\varepsilon}_n | \vec{y}_{1:n})$ and $p(\vec{\eta}_1, \dots, \vec{\eta}_n | \vec{y}_{1:n})$.[^simsmo]

[simsmo]: See @McCauslandMillerPelletier2011a for a comparison of the computational efficiency of various approaches.

### Mean correction simulation smoother

The mean-correction simulation smoother was introduced in @DurbinKoopman2002 . See @DurbinKoopman2012 (Sec 4.9) for an exposition of it. It requires only the previously described filters and smoothers, and generating samples from multivariate distributions.

#### Disturbances

1. Run a filter and disturbance smoother to calculate $\hat{\vec{\varepsilon}}_{1:n}$ and $\hat{\vec{\eta}}_{1:(n - 1)}$
2. Draw samples from the unconditional distribution of the disturbances,
    $$
    \begin{aligned}[t]
    \vec{\eta}^+_t &\sim N(0, \mat{H}_t) && t = 1, \dots, n - 1 \\
    \vec{\varepsilon}^+_t &\sim N(0, \mat{Q}_t) && t = 1, \dots, n
    \end{aligned}
    $$
3. Simulate observations from the system using the simulated disturbances,
    $$
    \begin{aligned}[t]
    \vec{y}^+_t &= \vec{d}_t + \mat{Z}_t \vec{\alpha}_t + \vec{\varepsilon}^+_t, \\
    \vec{\alpha}_{t + 1} &= \vec{c}_t + \mat{T}_t \vec{\alpha}_t + \mat{R}_t \vec{\eta}^+_t, \\
    \end{aligned}
    $$
    where $\vec{\alpha}_1 \sim N(\vec{a}_1, \mat{P}_1)$.
4. Run a filter and disturbance smoother on the simulated observations $\vec{y}^+$ to calculate $\hat{\vec{\varepsilon}}_t^+ = \E(\vec{\varepsilon}_t | \vec{y}^+_{1:n})$ and $\hat{\vec{\eta}}_t^+ = \E(\vec{\eta}_t | \vec{y}^+_{1:n})$.
5. A sample from $p(\hat{\vec{\eta}}_{1:(n - 1)}, \hat{\vec{\varepsilon}}_{1:n} | \vec{y}_{1:n})$ is
    $$
    \begin{aligned}[t]
    \tilde{\vec{\eta}}_t &= \vec{\eta}^+_t - \hat{\vec{\eta}}^+_t + \hat{\vec{\eta}}_t , \\
    \tilde{\vec{\varepsilon}}_t &= \vec{\varepsilon}^+_t - \hat{\vec{\varepsilon}}^+_t + \hat{\vec{\varepsilon}}_t .
    \end{aligned}
    $$


#### States

1. Run a filter and disturbance smoother to calculate the mean of the states conditional on the full series of observations, $\hat{\vec{\alpha}}_{1:n} = \E(\vec{\alpha}_{1:n} | \vec{y}_{1:n})$.
2. Draw samples from the unconditional distribution of the disturbances,
    $$
    \begin{aligned}[t]
    \vec{\eta}^+_t &\sim N(0, \mat{H}_t) && t = 1, \dots, n - 1 \\
    \vec{\varepsilon}^+_t &\sim N(0, \mat{Q}_t) && t = 1, \dots, n
    \end{aligned}
    $$
3. Simulate states and observations from the system using the simulated disturbances,
    $$
    \begin{aligned}[t]
    \vec{y}^+_t &= \vec{d}_t + \mat{Z}_t \vec{\alpha}_t + \vec{\varepsilon}^+_t, \\
    \vec{\alpha}^+_{t + 1} &= \vec{c}_t + \mat{T}_t \vec{\alpha}_t + \mat{R}_t \vec{\eta}^+_t, \\
    \end{aligned}
    $$
    where $\vec{\alpha}^+_1 \sim N(\vec{a}_1, \mat{P}_1)$.
4. Run a filter and smoother on the simulated observations $\vec{y}^+$ to calculate $\hat{\vec{\alpha}}_t^+ = \E(\vec{\alpha}_t | \vec{y}^+_{1:n})$.
5. A sample from $p(\hat{\vec{\alpha}}_{1:n} | \vec{y}_{1:n})$ is
    $$
    \begin{aligned}[t]
    \tilde{\vec{\alpha}}_t &= \vec{\alpha}^+_t - \hat{\vec{\alpha}}^+_t + \hat{\vec{\alpha}}_t .
    \end{aligned}
    $$

One convenient feature of this method is that since only the conditional means of the states are required, the fast state smoother can be used, since the variances of the states are not required.

### de Jong-Shephard method


These recursions were developed in @DeJongShephard1995 .
Although the the mean-correction simulation smoother will work in most cases, there are a few in which it will not work.

**TODO**

### Forward-filter backwards-smoother (FFBS)

This was the simulation method developed in @CarterKohn1994 and @Fruehwirth-Schnatter1994 .

**TODO**

## Missing observations

When all observations at time $t$ are missing, the filtering recursions become [@DurbinKoopman2012, Sec 4.10],
$$
\begin{aligned}[t]
\vec{a}_{t|t} &= \vec{a}_t , \\
\mat{P}_{t|t} &= \mat{P}_t , \\
\vec{a}_{t + 1} &= \mat{T}_t \vec{a}_t + \vec{c}_t \\
\mat{P}_{t + 1} &= \mat{T}_t \mat{P}_t \mat{T}_t\T + \mat{R}_t \mat{Q}_t \mat{R}_t\T
\end{aligned}
$$

This is equivalent to setting $\mat{Z}_t = \mat{0}$ (implying also that $\mat{K}_t = \mat{0}$) in the filtering equations.



For smoothing, also replace $\mat{Z}_t = \mat{0}$,
$$
\begin{aligned}[t]
\vec{r}_{t - 1} &= \mat{T}_t\T \vec{r}_t , \\
\mat{N}_{t - 1} &= \mat{T}_t\T \mat{N}_t \mat{T}_t,
\end{aligned}
$$

When some, but not all observations are missing, replace the observation equation by,
$$
\begin{aligned}[t]
\vec{y}^*_t &= \mat{Z}^*_t \vec{\alpha}_t + \vec{\varepsilon}_t^*, & \vec{\varepsilon}_t^* &\sim N(\vec{0}, \mat{H}_t^*),
\end{aligned}
$$
where,
$$
\begin{aligned}[t]
\vec{y}^*_t &= \mat{W}_t \vec{y}_t, \\
\mat{Z}^* &= \mat{W}_t \mat{Z}_t , \\
\vec{\varepsilon}_t &= \mat{W}_t \vec{\varepsilon}_t , \\
\mat{H}^*_t &= \mat{W}_t \mat{H}_t \mat{W}_t\T ,
\end{aligned}
$$
and $\mat{W}_t$ is a selection matrix to select non-missing values.
In smoothing the missing elements are estimated by the appropriate elements of $\mat{Z}_t \hat{\vec{alpha}}_t$, where $\hat{\vec{\alpha}}_t$ is the smoothed state.

If $y_{t,j}$ is missing, then setting the relevant entries in the forecast precision matrix, $F^{-1}_{t,j,.} = \vec{0}$ and $F^{-1}_{t,.,j} = \vec{0}$, and Kalman gain matrix, $K_{t,.,j} = \vec{0}$, will handle missing values without having to directly pass that information to the smoother.
When $K_{t,.,j} = 0$, it indicates that the observation $y_{t,j}$ provided no information.
However, it may be computationally more efficient if the values of the locations of the missing observations are known.

However, the simulation smoothers using the mean correction method need to have the missing values indicated, so filtering of the unconditional simulations accounts for those missing values. 
The missing values could be inferred from the values of $K_t$ or can be provided directly.


## Forecasting matrices

Forecasting future observations is the same as treating the future observations as missing [@DurbinKoopman2012, Sec 4.11],
$$
\begin{aligned}[t]
\bar{\vec{y}}_{n + j} &= \mat{Z}_{n + j} \bar{\vec{a}}_{n + j} \\
\bar{\mat{F}}_{n + j} &= \mat{Z}_{n + j} \bar{\mat{P}}_{n + j} \mat{Z}_{n + j}\T + \mat{H}_{n + j} .
\end{aligned}
$$
