///
/// # Common Smoother Functions
///

/**
Update $\vec{r}_t$ in smoothing recursions

@param vector r An $m \times 1$ vector with $\vec{r}_{t-1}$
@param matrix Z A $p \times m$ vector with $\mat{Z}_{t}$
@param vector v A $p \times 1$ vector of the forecast errors, $\vec{v}_t$.
@param matrix Finv A $p \times p$ matrix of the forecast precision, $\mat{F}^{-1}_t$.
@param matrix L An $m \times m$ matrix with $\mat{L}_t$.
@return matrix An $m \times 1$ vector with $\vec{r}_t$.

In smoothing recursions, the vector $\vec{r}_t$ is updated with,
$$
\vec{r}_t = \mat{Z}' \mat{F}^{-1}_t \vec{v}_t + \mat{L}' \vec{r}_{t - 1} .
$$

See [@DurbinKoopman2012, p. 91]
*/
vector ssm_smooth_update_r(vector r, matrix Z, vector v, matrix Finv,
                           matrix L) {
  vector[num_elements(r)] r_new;
  r_new = Z ' * Finv * v + L ' * r;
  return r_new;
}

/**
Update $\mat{N}_t$ in smoothing recursions

@param vector N An $m \times 1$ vector with $\vec{N}_{t-1}$
@param matrix Z A $p \times m$ vector with $\mat{Z}_{t}$
@param matrix Finv A $p \times p$ matrix of the forecast precision, $\mat{F}^{-1}_t$.
@param matrix L An $m \times m$ matrix with $\mat{L}_t$.
@return matrix An $m \times m$ matrix with $\vec{N}_t$.

In smoothing recursions, the matrix $\vec{N}_t$ is updated with,
$$
\mat{N}_t = \mat{Z})_t' \mat{F}^{-1}_t \mat{Z}_t + \mat{L}_t' \mat{N}_t \mat{L}_t .
$$

See [@DurbinKoopman2012, p. 91]
*/
matrix ssm_smooth_update_N(matrix N, matrix Z, matrix Finv, matrix L) {
  matrix[rows(N), cols(N)] N_new;
  N_new = quad_form(Finv, Z) + quad_form(N, L);
  return N_new;
}
