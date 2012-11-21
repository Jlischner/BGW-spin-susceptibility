function out = cI(in);
global gbl_S;

out=zeros(size(in));
for col=1:size(in,2);
out(:,col)=fft3(in(:,col),gbl_S,1);
endfor;

endfunction;
