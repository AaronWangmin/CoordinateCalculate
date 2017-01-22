%-----------------                          L2Gausee              -------------
% transfer the Geodetic coordinates the Gauss plane rectangular coordinates
% args   : 
%          B,L,H		    I	Longitude(Lng),Latitude(Lat), Gedetic Height
%		   a,f      		I   semi-major axis,flattening
%          L0               I	central meridian
% return : X,Y,Z 		    O   spatial rectangular coordinate 
% notes  : 
%------------------------------------------------------------------------------



function [x,y,h]  = BL2Gauss(B,L,H,a,f,L0)

b = (1-f) * a;
sqrt_e = (a^2 - b^2) /a^2;

delta_L = L - L0;
N = a /(1 -sqrt_e *( sin(B) )^2) ^0.5;
t = tan(B);
sqrt_e1 = (a^2 - b^2) /b^2;
sqrt_eta = sqrt_e1 * ( cos(B) )^2;

C1 = 1 + 3/4 * sqrt_e + 45/64 * sqrt_e^2 + 175/256 * sqrt_e^3 + 11025/16384 * sqrt_e^4;
C2 = 3/4 * sqrt_e + 15/16 * sqrt_e^2 + 525/512 * sqrt_e^3 + 2205/2048 * sqrt_e^4;
C3 = 15/64 * sqrt_e^2 + 105/256 * sqrt_e^3 + 2205/4096 * sqrt_e^4;
C4 = 35/512 * sqrt_e^3 + 315/2048 * sqrt_e^4;
C5 = 315/131072 * sqrt_e^4;

X = a * (1 - sqrt_e) *( C1 * B - 1/2 * C2 * sin( 2*B ) + 1/4 * C3 * sin( 4* B) - 1/6 *C4 * sin( 6*B ) + C5 * sin( 8*B ) );

x = X + 1/2 * N * sin( B ) * cos ( B ) * ( delta_L )^2 * ( 1 + 1/12 * ( delta_L )^2 * ( cos (B) )^2 * ( 5 - t^2 + 9 * sqrt_eta + 4 * sqrt_eta^2 ) + 1/360 * delta_L^4 * ( cos (B) )^4 * ( 61 - 58 * t^2 + t^4)  );

y = 500000 + N * cos( B ) * delta_L * ( 1 + 1/6 * delta_L^2 * ( cos( B) )^2 * ( 1 - t^2 + sqrt_eta ) + 1/120 * delta_L^4 * ( cos( B ))^4 * ( 5 - 18 * t^2 + t^4 - 14 * sqrt_eta - 58 * sqrt_eta * t^2)  ); 

h = H;




