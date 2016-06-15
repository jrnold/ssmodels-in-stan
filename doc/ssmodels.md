## Models

A linear Gaussian state space model (SSM)[^dlm] is defined by

$$
\begin{aligned}[t]
\vec{y}_t &= \mat{Z}_t \vec{\alpha}_t + \vec{\varepsilon}_t,  &
\vec{\varepsilon}_t & \sim N(0, \mat{H}),
\end{aligned}
$$

The *state equation* is defined as,
$$
\begin{aligned}[t]
\vec{\alpha}_{t + 1} &= \mat{T}_t \vec{\alpha}_t + \mat{R}_t \vec{\eta}_t,  &
\vec{\eta}_t & \sim N(0, \mat{Q}),
\end{aligned}
$$

The initial state, $\alpha_1$ is normally distributed such that,
$$
\alpha_1 \sim N(\vec{a}_1, \vec{P}_1)
$$

[^dlm]: This is also called a dynamic linear model (DLM).

## Integrated Sampler

See
$$
p(\vec{y} | \mat{H}, \mat{\Psi}) = \int p(\vec{y} | \vec{\alpha}, \mat{H}) p(\vec{\alph} | \mat{\Psi}) d\,\vec{\alpha},
$$

The log-likelihood of a state-space function can be calculated analytically and expressed up to an integrating constant marginally of the state vector $\vec{\alpha}$,
$$
\lot p(\vec{y} | \mat{H}, \mat{\Psi}) = \text{const} - 0.5 \left( \sum_{t = 1}^n \sum_{i = 1}^p \ln f_{t, i} - v_t^2 f^{-1}_{t,i} \right) .
$$
Thus, parameters of the state space model can be sampled as,
$$
p(\mat{H}, \mat{\Psi} | \vec{y}) \propto p(\vec{y} | \mat{H}, \mat{\Psi}) p(\mat{H}) p(\mat{\Psi}) .
$$

##
