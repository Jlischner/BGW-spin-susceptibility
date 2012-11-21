%# out=fft3(dat,N,s) - computes 3d fft (dimensions in N) of sign s,
%# out(l,m,n)=sum_{a,b,c} exp(2 pi s i *(a*l/N(1)+b*m/N(2)+c*n/N(3))*in(a,b,c).
%# Notes: 1) fortran/matlab ordering assumed, ordering in mem is out(1,1,1),
%#           out(2,1,1), ..., out(N(1),1,1), out(1,2,1), ...
%#        2) a=fft3(b,N,1) => b=fft3(a,N,-1)/prod(N); ie., fft(dat,N,1) and
%#           fft(dat,N,-1) are inverses except for normalization by
%#           N(1)*N(2)*N(3)
function out=fft3(dat,N,s)
%  tic;
  if s==1
    out=reshape(ifftn(reshape(dat,N(1),N(2),N(3)))*prod(N),size(dat));
  else
    out=reshape(fftn(reshape(dat,N(1),N(2),N(3))),size(dat));
  end

%  Ntot=prod(size(dat));
%  MFLOPS=5*Ntot*log2(Ntot)/1e6/toc
endfunction
