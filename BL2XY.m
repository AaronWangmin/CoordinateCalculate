%-----------------                           BL2XY                -------------
% transfer the spatial rectangular coordinates to the Geodetic coordinates
% args   : 
%          B,L,H		    I	Longitude(Lng),Latitude(Lat),
%		   a,f      		I   semi-major axis,flattening
% return : X,Y,Z 		    O   spatial rectangular coordinate 
% notes  : 
%------------------------------------------------------------------------------

function [X,Y,Z] = BL2XY(B,L,H,a,f)

b = (1-f) * a;
sqrt_e = (a^2 - b^2) /a^2;
N = a /(1-sqrt_e *( sin(B) )^2) ^0.5;

X = (N+H) *cos(B) *cos(L);
Y = (N+H) *cos(B) *sin(L);
Z = ( N *(1- sqrt_e) + H) *sin(B);


