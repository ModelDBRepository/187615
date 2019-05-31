function dvs = membpotderivs(time,vrec)
  N = length(time);
  tdiff = time(2:N)-time(1:N-1);
  vdiff = vrec(2:N)-vrec(1:N-1);
  mderiv = vdiff./tdiff;
  dvs = 0.5*(mderiv(2:N-1)+mderiv(1:N-2));
