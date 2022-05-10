function [ SCoeff ] = linf( x, y )
% try to find the coefficient SCoeff such that, x is closet to SCoeff * y by least
% square approach

% % the boundaries of Scale Coefficient
% a = 0.10;
% b = 10.00;
% 
% % the step of scale coefficient
% dx = 0.1;
% SCoeff = a:dx:b;




% create the funciton handle
f = @ (a) norm ( ( abs (  a(1) * x + a(2) - y ) ) );
a0=[1 0];
SCoeff = fminsearch ( f, a0 );





end

