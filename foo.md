$$
L_T = \Pr\left(X^{(T)} = x^{(T)}\right)  = \delta P(x_1) \Gamma P(x_2) \cdots \Gamma P(x_T) 1'
$$

where $\delta$ is the initial distribution (distribution of $C_1$), and $P(x)$ is the $m \times m$ diagonal matrix where the $i$th diagonal is the state-dependent probability $p_i(x)$.
Then $L_T = \alpha_T 1'$
$$
\alpha_1 = \delta P(x_1)
$$
and
$$
\begin{aligned}[t]
\alpha_t = \alpha_{t-1} \Gamma P(x_t) & \text{for $t = 2, 3, \dots, T$.}
\end{aligned}
$$
If the chain is stationary $\delta = \delta \Gamma$,
$$
\alpha_0 = \delta
$$
and
\begin{aligned}[t]
\alpha_t = \alpha_{t-1} \Gamma P(x_t) & \text{for $t = 1, 2, 3, \dots, T$.}
\end{aligned}
