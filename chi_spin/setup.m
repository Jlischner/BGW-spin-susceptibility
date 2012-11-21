global gbl_S; 
global gbl_R; 
global gbl_G;
global gbl_G2;
global gbl_r;
global gbl_dV;
global gbl_Vcell;


S=[ 60; 60; 360]
R=[ 1.00000     0.0000      0.0000;
    0.00000     1.0000      0.0000;
    0.00000     0.0000      5.8624 ]

R *= 6.96787;

Vcell=abs(det(R));
dV=abs(det(R))/prod(S);
lenS = prod(S);

ms=[0:prod(S)-1];
ms=ms.';
m1=rem(ms,S(1));
m2=rem(floor(ms/S(1)),S(2));
m3=rem(floor(ms/(S(1)*S(2))),S(3));

M=[m1, m2, m3];
n1=m1-(m1>S(1)/2)*S(1);
n2=m2-(m2>S(2)/2)*S(2);
n3=m3-(m3>S(3)/2)*S(3);
N=[n1, n2, n3];

r=M*inverse(diag(S))*R.';
G=2.*pi*N*inverse(R);
G2=sum(G.^2,2);

gbl_S=S;
gbl_R=R;
gbl_G=G;
gbl_G2=G2;
gbl_r=r;
gbl_dV=dV;
gbl_Vcell=Vcell;
