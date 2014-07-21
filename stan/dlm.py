#!/usr/bin/env python3
"""
Functions to generate Stan code for Kalman filtering and smoothing / simulation
"""
import string
import random
import itertools


# Raymond Hettinger
# http://code.activestate.com/recipes/273085-sample-with-replacement/
def sample_wr(population, k):
    "Chooses k random elements (with replacement) from a population"
    n = len(population)
    _random, _int = random.random, int  # speed hack
    return [population[_int(_random() * n)] for i in itertools.repeat(None, k)]


def random_prefix(k = 5, prefix = "", postfix = ""):
    """ Generate a random prefix for Stan variable """
    LETTERS = string.ascii_letters
    ALL = LETTERS + string.digits + "_"
    ret = random.choice(LETTERS)
    if k > 1:
        ret = ret + ''.join(sample_wr(ALL, k - 1))
    return prefix + ret + postfix


def indentlines(x, indent = 0):
    """ Indent string """
    return '\n'.join([" " * max(int(indent), 0) + y for y in x.split("\n")])


def kalman_seq(y = "y",
               y_obs = "y_obs",
               n = "n", r = "r", p = "p",
               F = "F", G = "G",
               V = "V", W = "W",
               b = "b", g = "g",
               loglik_obs = "loglik_obs",
               m0 = "m0", C0 = "C0",
               a = "a", R = "R",
               f = "f", Q = "Q",
               m = "m", C = "C",
               Qfill = "1e7",
               tv = False,
               prefix = "kalman_seq_",
               indent = 0):
    """ Generate Stan code for Kalman Filter (sequential, missing values)
    """

    TEMPLATE = """
{
  real $err;
  vector[$p] $K;
  matrix[$p, $p] $J;
  vector[$p] $m_tmp;
  matrix[$p, $p] $C_tmp;
  vector[$p] $Fj;
  matrix[$p] $I;
  $m[1] <- $m0;
  $C[1] <- $C0;
  $I <- diag_matrix(rep_vector($p, 1.0));
  for ($t in 1:$n) {
    $a[$t] <- $g + $G * $m[$t];
    $R[$t] <- quad_form($C[$t], $G ') + $W;
    $m_tmp <- $a[$t];
    $C_tmp <- $R[$t];
    for ($j in 1:$r) {
      if (step($y_obs[$t, $j])) {
        $Fj <- row($F, $j) ';
        $f[$t, $j] <- $b[$j] + dot_product($Fj, $m_tmp);
        $Q[$t, $j] <- $Fj ' * $C_tmp * $Fj + $V[$j]; 
        $err <- $y[$t, $j] - $f[$t, $j];
        $K <- $C_tmp * $Fj / $Q[$t, $j];
        $m_tmp <- $m_tmp + $K * $err;
        $J <- ($I - $K * $Fj ');
        $C_tmp <- quad_form($C_tmp, $J ') + $K ' * $K * $V[$j];
        $loglik_obs[$t, $j] <- -0.5 * (log(2 * pi())
                                       + log($Q[$t, $j])
                                       + pow($err, 2) / $Q[$t, $j]);
      } else {
        $f[$t, $j] <- 0.0;
        $Q[$t, $j] <- $Qfill;
        $loglik_obs[$t, $j] <- 0.0;
      }
    }
    $m[$t] <- $m_tmp;
    $C[$t] <- $C_tmp;
  }
}
"""

    def mangle(x):
        return "%s%s" % (prefix, x)

    def iftv(x):
        if tv:
            return "%s[%s]" % (x, mangle("t"))
        else:
            return x

    # time variable
    data = {
        't' : mangle("t"),
        'j' : mangle("j"),
        'err' : mangle("err"),
        'K' : mangle("K"),
        'J' : mangle("J"),
        'm_tmp' : mangle("m_tmp"),
        'C_tmp' : mangle("C_tmp"),
        'I' : mangle("I"),
        'Fj' : mangle("Fj"),
        'y' : y, 'y_obs' : y_obs,
        'n' : n, 'r' : r, 'p' : p,
        'F' : iftv(F), 'G' : iftv(G), 
        'V' : iftv(V), 'W' : iftv(W),
        'b' : iftv(b), 'g' : iftv(g),
        'loglik_obs' : loglik_obs,
        'm0' : m0, 'C0' : C0,
        'a' : a, 'R' : R,
        'f' : f, 'Q' : Q,
        'm' : m, 'C' : C,
        'Qfill' : Qfill
    }
    return indentlines(string.Template(TEMPLATE).substitute(**data), indent)


def kalman_seq_nonmiss(y = "y",
                       y_obs = "y_obs",
                       n = "n", r = "r", p = "p",
                       F = "F", G = "G",
                       V = "V", W = "W",
                       b = "b", g = "g",
                       loglik_obs = "loglik_obs",
                       m0 = "m0", C0 = "C0",
                       a = "a", R = "R",
                       f = "f", Q = "Q",
                       m = "m", C = "C",
                       Qfill = "1e7",
                       tv = False,
                       prefix = "kalman_seq_nonmiss_",
                       indent = 0):

    """ Generate Stan Code for Kalman Filtering (sequential, no-missing values)
    """

    TEMPLATE = """
{
  real $err;
  vector[$p] $K;
  matrix[$p, $p] $J;
  vector[$p] $m_tmp;
  matrix[$p, $p] $C_tmp;
  vector[$p] $Fj;
  matrix[$p] $I;
  $m[1] <- $m0;
  $C[1] <- $C0;
  $I <- diag_matrix(rep_vector($p, 1.0));
  for ($t in 1:$n) {
    $a[$t] <- $g + $G * $m[$t];
    $R[$t] <- quad_form($C[$t], $G ') + $W;
    $m_tmp <- $a[$t];
    $C_tmp <- $R[$t];
    for ($j in 1:$r) {
      $Fj <- row($F, $j) ';
      $f[$t, $j] <- $b[$j] + dot_product($Fj, $m_tmp);
      $Q[$t, $j] <- $Fj ' * $C_tmp * $Fj + $V[$j]; 
      $err <- $y[$t, $j] - $f[$t, $j];
      $K <- $C_tmp * $Fj / $Q[$t, $j];
      $m_tmp <- $m_tmp + $K * $err;
      $J <- ($I - $K * $Fj ');
      $C_tmp <- quad_form($C_tmp, $J ') + $K ' * $K * $V[$j];
      $loglik_obs[$t, $j] <- - 0.5 * (log(2 * pi())
                                      + log($Q[$t, $j])
                                      + pow($err, 2) / $Q[$t, $j]);
    }
    $m[$t] <- $m_tmp;
    $C[$t] <- $C_tmp;
  }
}
"""
    
    def mangle(x):
        return "%s%s" % (prefix, x)

    def iftv(x):
        if tv:
            return "%s[%s]" % (x, mangle("t"))
        else:
            return x

    # time variable
    data = {
        't' : mangle("t"),
        'j' : mangle("j"),
        'err' : mangle("err"),
        'K' : mangle("K"),
        'J' : mangle("J"),
        'm_tmp' : mangle("m_tmp"),
        'C_tmp' : mangle("C_tmp"),
        'I' : mangle("I"),
        'Fj' : mangle("Fj"),
        'y' : y,
        'n' : n, 'r' : r, 'p' : p,
        'F' : iftv(F), 'G' : iftv(G), 
        'V' : iftv(V), 'W' : iftv(W),
        'b' : iftv(b), 'g' : iftv(g),
        'loglik_obs' : loglik_obs,
        'm0' : m0, 'C0' : C0,
        'a' : a, 'R' : R,
        'f' : f, 'Q' : Q,
        'm' : m, 'C' : C
    }
    return indentlines(string.Template(TEMPLATE).substitute(**data), indent)


def kalman_batch_nomiss(y = "y",
                        n = "n", r = "r", p = "p",
                        F = "F", G = "G",
                        V = "V", W = "W",
                        b = "b", g = "g",
                        loglik_obs = "loglik_obs",
                        m0 = "m0", C0 = "C0",
                        a = "a", R = "R",
                        f = "f", Q = "Q",
                        m = "m", C = "C",
                        tv = False,
                        prefix = "kalman_batch_nomiss_",
                        indent = 0):

    """ Generate Stan Code for Kalman Filtering (batch, no-missing values)
    """

    TEMPLATE = """
{
  matrix[$p, $p] $I;
  $I <- diag_matrix(rep_vector($p, 1.0));
  $m[1] <- $m0;
  $C[1] <- $C0;
  for ($t in 1:$n) {
    vector[$r] $err;
    matrix[$p, $r] $K;
    matrix[$r, $r] $Q_tmp;      
    matrix[$r, $r] $Qinv;
    matrix[$p, $p] $J;
    matrix[$p, $p] $R_tmp;
    $a[$t] <- $g + $G * $m[$t];
    $R_tmp <- quad_form($C[$t], $G ') + $W;
    $R[$t] <- 0.5 * ($R_tmp + $R_tmp ');
    $f[$t] <- $b + $F * $a[$t];
    $Q_tmp <- quad_form($R[$t], $F ') + $V;
    $Q[$t] <- 0.5 * ($Q_tmp + $Q_tmp ');
    $err <- $y[$t] - $f[$t];
    $Qinv <- inverse_spd($Q[$t]);
    $K <- $R[$t] * $F ' * $Qinv;
    $m[$t + 1] <- $a[$t] + $K * $err;
    $J <- ($I - $K * $F);
    $C[$t + 1] <- quad_form($R[$t], $J ') + quad_form($V, $K ');
    $loglik_obs[$t] <- - 0.5 * ($r * log(2 * pi())
                                + log_determinant($Q_tmp)
                                + quad_form($Qinv, $err));
  }
}
"""
    
    def mangle(x):
        return "%s%s" % (prefix, x)

    def iftv(x):
        if tv:
            return "%s[%s]" % (x, mangle("t"))
        else:
            return x

    # time variable
    data = {
        'loglik_obs' : loglik_obs,
        'y' : y, 
        'n' : n, 'r' : r, 'p' : p,
        'F' : iftv(F), 'G' : iftv(G), 
        'V' : iftv(V), 'W' : iftv(W),
        'b' : iftv(b), 'g' : iftv(g),
        'm0' : m0, 'C0' : C0,
        'a' : a, 'R' : R,
        'f' : f, 'Q' : Q,
        'm' : m, 'C' : C,
        't' : mangle("t"),
        'I' : mangle("I"),
        'err' : mangle("err"),
        'K' : mangle("K"),
        'Q_tmp' : mangle("Q_tmp"),
        'Qinv' : mangle("Qinv"),
        'J' : mangle("J"),
        'R_tmp' : mangle("R_tmp")
    }
    return indentlines(string.Template(TEMPLATE).substitute(**data), indent)


def kalman_batch(y = "y",
                 nobs = "nobs",
                 obs_ind = "obs_ind",
                 n = "n", r = "r", p = "p",
                 F = "F", G = "G",
                 V = "V", W = "W",
                 b = "b", g = "g",
                 m0 = "m0", C0 = "C0",
                 loglik_obs = "loglik_obs",
                 a = "a", R = "R",
                 f = "f", Q = "Q",
                 m = "m", C = "C",
                 Qmiss = "1e7",
                 tv = False,
                 prefix = "kalman_batch_",
                 indent = 0):
    """ Generate Stan code for Kalman Filter (batch, missing values)
    """

    TEMPLATE = """
{
  matrix[$p, $p] $I;
  $I <- diag_matrix(rep_vector(1, p));
  // set initial states
  $m[1] <- $m0;
  $C[1] <- $C0;

  // loop through observations    
  for ($t in 1:$n) {
    vector[$nobs[$t]] $y_tmp;
    matrix[$nobs[$t], $p] $F_tmp;
    vector[$nobs[$t]] $b_tmp;
    vector[$nobs[$t], $nobs[$t]] $V_tmp;
    vector[$nobs[$t]] $f_tmp;      
    matrix[$nobs[$t], $nobs[$t]] $Q_tmp;
    vector[$nobs[$t]] $err;
    matrix[$p, $nobs[$t]] $K;
    matrix[$nobs[$t], $nobs[$t]] $Qinv;
    matrix[$p, $p] $J;

    // one step ahead predictive distribution of \theta_t | y_{1:(t-1)}
    $a[$t] <- $g + $G * $m[$t];
    $R[$t] <- quad_form($C[t], $G ') + $W;
    $f[$t] <- rep_vector(0.0, $r);
    $Q[t] <- diag_matrix(rep_vector($Qmiss, $r));
    if ($nobs[$t] == 0) {
      // all observations missing
      $m[$t + 1] <- $a[$t];
      $R[$t + 1] <- $R[$t];
      $loglik_obs[$t] <- 0.0;
    } else {
      if ($nobs[$t] == $r) {
        // no observations missing
        $y_tmp <- $y;
        $b_tmp <- $b;
        $F_tmp <- $F;
        $V_tmp <- $V;
      } else {
        // at least one observation missing
        for ($i in 1:$nobs[$t]) {
          $y_tmp[$i] <- y[$obs_ind[$i]];
          $b_tmp[$i] <- y[$obs_ind[$i]];
          for ($j in 1:$p) {
            $F_tmp[$i, $j] <- $F[$obs_ind[$i], $j];
          }
          for ($j in 1:$nobs[$t]) {
            $V_tmp[$i, $j] <- $V[$obs_ind[$i], $obs_ind[$j]];
          }
        }
      }
      // one step ahead predictive distribution of y_t | y_{1:(t-1)}
      $f_tmp <- $b_tmp + $F_tmp * $a[$t];
      $Q_tmp <- quad_form($R[$t], $F_tmp ') + $V_tmp;
      // forecast error
      $err <- $y_tmp[$t] - $f_tmp[$t];
      // Kalman gain
      $Qinv <- inverse_spd($Q[$t]);
      $K <- R[t] * $F_tmp ' * $Qinv;
      // posterior distribution of \theta_t | y_{1:t}
      $m[$t + 1] <- $a[$t] + $K * $err;
      // matrix used in Joseph stabilized form
      $J <- ($I - $K * $F_tmp);
      $C[$t + 1] <- quad_form($R[$t], $J ') + quad_form($V_tmp, $K ');
      // log likelihood
      $loglik_obs[$t] <- - 0.5 * ($nobs[$t] * log(2 * pi())
                                 + log_determinant($Q_tmp)
                                 + quad_form($Qinv, $err));
      if ($nobs[$i] == $r) {
        $f[$t] <- $f_tmp;
        $Q[$t] <- $Q_tmp;
      } else {
        for (i in 1:$nobs[$t]) {
          $f[$obs_ind[$i]] <- $f_tmp[$i];
          for ($j in 1:$nobs[t]) {
            $Q[$obs_ind[$i], $obs_ind[$j]] <- $Q_tmp[$i, $j];
          }
        }
      }
    }
  }
}
"""
    
    def mangle(x):
        return "%s%s" % (prefix, x)

    def iftv(x):
        if tv:
            return "%s[%s]" % (x, mangle("t"))
        else:
            return x

    # time variable
    data = {
        'y': y,
        'nobs' : nobs,
        'obs_ind' : obs_ind,
        'n' : n, 'r' : r, 'p' : p,
        'F' : iftv(F), 'G' : iftv(G),
        'V' : iftv(V), 'W' : iftv(W),
        'b' : iftv(b), 'g' : iftv(g),
        'm0' : m0, 'C0' : C0,
        'loglik_obs' : loglik_obs,
        'a' : a, 'R' : R,
        'f' : f, 'Q' : Q,
        'm' : m, 'C' : C,
        'Qmiss' : Qmiss,
        # temporary
        'I' : mangle("I"),
        'y_tmp' : mangle("y_tmp"),
        'F_tmp' : mangle("F_tmp"),        
        'b_tmp' : mangle("b_tmp"),
        'V_tmp' : mangle("v_tmp"),
        'f_tmp' : mangle("f_tmp"),
        'Q_tmp' : mangle("Q_tmp"),
        'K' : mangle("K"),
        'err' : mangle("err"),
        'Qinv' : mangle("Qinv"),
        'J' : mangle("J"),
        't' : mangle("t"),        
        'i' : mangle("i"),
        'j' : mangle("j"),
    }
    return indentlines(string.Template(TEMPLATE).\
                       substitute(**data), indent)


def kalman_sim(theta = "theta",
               n = "n", p = "p",
               G = "G",
               m = "m", C = "C",
               a = "a", R = "R",
               tv = False,
               prefix = "kalman_sim_",
               indent = 0):
    """ Generate Stan Code for Backward Sampling
    """

    TEMPLATE = """
theta[$n + 1] <- multi_normal_rng(m[$n + 1], C[$n + 1]);
for ($i in 1:$n) {
  int $t;
  vector[$p] $h;
  matrix[$p, $p] $H;
  matrix[$p, $p] $Rinv;
  $t <- $n - $i + 1;
  $Rinv <- inverse($R[$t]);
  // sample 
  $h <- $m[$t] + $C[$t] * $G ' * $Rinv * ($theta[$t + 1] - $a[$t]);
  $H <- $C[$t] - $C[$t] * $G ' * $Rinv * $G * C[$t];
  $theta[$t] <- multi_normal_rng($h, 0.5 * ($H + $H '));
}
"""

    def mangle(x):
        return "%s%s" % (prefix, x)

    def iftv(x):
        if tv:
            return "%s[%s]" % (x, mangle("t"))
        else:
            return x
    
    data = {
        't' : mangle("t"),
        'i' : mangle("i"),
        'h' : mangle("h"),
        'H' : mangle("H"),
        'Rinv' : mangle("Rinv"),
        'theta' : theta,
        'n' : n, 'p': p,
        'G' : iftv(G),
        'a' : a, 'R' : R, 'm' : m, 'C' : C,
    }
    return indentlines(string.Template(TEMPLATE).substitute(**data), indent)    


def kalman_smooth(s = "s", S = "S",
                  n = "n", p = "p",
                  G = "G",
                  m = "m", C = "C",
                  a = "a", R = "R",
                  tv = False,
                  prefix = "kalman_smooth_",
                  indent = 0):
    """ Generate Stan Code for smoothed state estimates
    """

    TEMPLATE = """
$s[$n + 1] <- $m[$n + 1]; 
$S[$n + 1] <- $C[$n + 1];
for ($i in 1:$n) {
  int $t;
  matrix[$p, $p] $Rinv;
  matrix[$p, $p] $S_tmp
  $t <- $n - $i + 1;
  $Rinv <- inverse($R[$t]);
  $s[$t] <- $m[$t] + $C[$t] * $G ' * $Rinv * ($s[$t + 1] - $a[$t]);
  $S_tmp <- $C[$t] - $C[$t] * $G ' * $Rinv * ($R[$t] - $S[$t + 1]) * $Rinv * $G * $C[$t];
  $S[$t] <- 0.5 * ($S_tmp ' + $S_tmp);
}
"""
    def mangle(x):
        return "%s%s" % (prefix, x)

    def iftv(x):
        if tv:
            return "%s[%s]" % (x, mangle("t"))
        else:
            return x
    
    data = {
        't' : mangle("t"),
        'i' : mangle("i"),
        'Rinv' : mangle("Rinv"),
        'S_tmp' : S_tmp,
        'n' : n, 'p': p,
        'G' : iftv(G),
        'a' : a, 'R' : R, 'm' : m, 'C' : C,
    }
    return indentlines(string.Template(TEMPLATE).substitute(**data), indent)    

