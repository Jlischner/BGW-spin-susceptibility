function indx = FFTbox(gkvectors,S);

  Ngk = size( gkvectors)(1);
  lenS = prod(S);
  ivec = gkvectors + ones(Ngk,1)*[1 1 1];
  
  ivec(:,1) += S(1)*(ivec(:,1)<=0) ;
  ivec(:,2) += S(2)*(ivec(:,2)<=0) ;
  ivec(:,3) += S(3)*(ivec(:,3)<=0) ;
  indx = S(2)*S(1)*( ivec(:,3) -1) + S(1)*( ivec(:,2)-1) + ivec(:,1);

endfunction;
