more off;
setup;

Nfreq = 12;
nmtxdat = load("nmtx.dat");
%# factor of two - see Rousseau paper
chi0 = load("fullchi.dat")/2;
vcdat  = load("vcoul");

Nq = length(nmtxdat);
chi0 = chi0(:,1) + I*chi0(:,2);

%# calculate Ixc:
getIxc;

chi0h = zeros(Nfreq,Nq);
chih  = zeros(Nfreq,Nq);

for iq = 1:Nq;

  printf("doing %d of %d qpoints \n",iq,Nq);
  nmtx  = nmtxdat(iq);
  ind1 = sum(nmtxdat(1:iq-1));
  ind2 = sum(nmtxdat(1:iq)) ;
  gvecs = vcdat( ind1+1 : ind2 , 4:6 );  
  indx = FFTbox(gvecs,S);

  ind1 = sum(nmtxdat(1:iq-1).^2);
  ind2 = sum(nmtxdat(1:iq).^2) ;
  chi0w = chi0( ind1*Nfreq+1 : ind2*Nfreq ) ;
  chi0w = reshape(chi0w,Nfreq,nmtx^2);

  for ifreq = 1:Nfreq;
    printf("ifeq = %d \n",ifreq);
    chi0f = reshape( chi0w(ifreq,:),nmtx,nmtx);

    epsmat = zeros(nmtx,nmtx);
    for ncol = 1:nmtx;
      printf("col %d \n",ncol);
      %# transform each col of chi0 to real space
      %# multiply by Ixc and transform back
      row_chi0 = zeros(lenS,1);
      row_chi0(indx) = chi0f(:,ncol);
      row_chi0   = cJ( row_chi0);
      epsmat(:,ncol) = cI( row_chi0 .* Ixc)(indx);
      
    endfor;
    
    epsmat = eye(nmtx) - epsmat;
    chi = inv(epsmat) * chi0f;
    
    chi0h(ifreq,iq) = chi0f(1,1);
    chih(ifreq,iq)  = chi(1,1);
    
  endfor;
endfor;
 
more on;