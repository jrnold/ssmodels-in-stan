functions {
  #include ssm.stan
}
data {
  int n;
  int x[n];
}
model {}
generated quantities {
  int output;
  output = int_sum_false(x);
}
