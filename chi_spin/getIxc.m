%# calculate fxc needed for interaction spin susceptibility
%# fxc(r,r') = d^2 Exc / dm(r)dm(r') = delta(r-r') Ixc(n(r))
%# Ixc(n) = d^2 (n e_xc(n,m)) /dm^2
%# Ixc(n) = 1/n ( alpha_c(n) + expp )

%# run setup.m
#setup;

%# parameter for -ac from Perdew Wang paper
%# last column of TABLE I (atomic units)
A  = 0.016887;
a1 = 0.11125;
b1 = 10.357;
b2 = 3.6231;
b3 = 0.88026;
b4 = 0.49671;

%# load in density in real space from espresso
n  = load("cd.dat");
#ninv = 1./n;
#ninv(ninv > 75) = 75;
#n = 1./ninv;

rs = (3/4/pi./n).^(1/3);
rrs= sqrt(rs);

%# calculate spin stiffness
ac   = log(1+1./(2*A*( b1*rrs + b2*rrs.^2 + b3*rrs.^3 + b4*rrs.^4)));
ac .*= 2*A*(1+a1.*rs);

%# calculate 2nd derivative of exchange energy
expp  = -3/4/pi./rs * (9*pi/4)^(1/3);
expp *= 4/9;

%# final result for I(n)
Ixc = 2./n .* (ac + expp);

#x0 = 100;
#x = 1./n;
#Ixc = Ixc.*(1 - exp(-x0^2./x.^2));

%# Rydberg: multiply by 2
Ixc *= 2;
%# Fourier transform I(n(r)) to Fourier space
#Ixc = cJ(Ixc);

