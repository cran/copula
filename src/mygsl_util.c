double gsl_fdiv (const double x, const double y)
{
  return x / y;
}

double gsl_nan (void)
{
  return gsl_fdiv (0.0, 0.0);
}
