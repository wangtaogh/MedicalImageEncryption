function [T, K]=FOChen(parameters, orders, TSim, Y0)
%
% Numerical Solution of the Fractional-Order Chen's System
%
%   D^q1 x(t) = a(y(t)-x(t))
%   D^q2 y(t) = dx(t) - x(t)z(t) + cy(t)
%   D^q3 z(t) = x(t)y(t) - bz(t)
%
% function [T, Y] = FOChen(parameters, orders, TSim, Y0)
%
% Input:    parameters - model parameters [a, b, c, d]
%           orders - derivatives orders [q1, q2, q3]
%           TSim - simulation time (0 - TSim) in sec
%           Y0 - initial conditions [Y0(1), Y0(2), Y0(3)]
%
% Output:   T - simulation time (0 : Tstep : TSim)
%           Y - solution of the system (x=Y(1), y=Y(2), z=Y(3))
%
% Author:  (c) Ivo Petras (ivo.petras@tuke.sk), 2010.
%

% time step:
h=0.005; 
% number of calculated mesh points:
n=round(TSim/h);
%orders of derivatives, respectively:
q1=orders(1); q2=orders(2); q3=orders(3);
% constants of Chen's system:
a=parameters(1); b=parameters(2); 
c=parameters(3); d=c-a;
% binomial coefficients calculation:
cp1=1; cp2=1; cp3=1;
for j=1:n
    c1(j)=(1-(1+q1)/j)*cp1;
    c2(j)=(1-(1+q2)/j)*cp2;
    c3(j)=(1-(1+q3)/j)*cp3;
    cp1=c1(j); cp2=c2(j); cp3=c3(j);
end
% initial conditions setting:
x_chen(1)=Y0(1); y_chen(1)=Y0(2); z_chen(1)=Y0(3);
% calculation of phase portraits /numerical solution/:
for i=2:n
    x_chen(i)=(a*(y_chen(i-1)-x_chen(i-1)))*h^q1 - memo(x_chen, c1, i);
    y_chen(i)=(-d*x_chen(i)-x_chen(i)*z_chen(i-1)+c*y_chen(i-1))*h^q2 - memo(y_chen, c2, i);
    z_chen(i)=(x_chen(i)*y_chen(i)-b*z_chen(i-1))*h^q3 - memo(z_chen, c3, i);
end
for j=1:n
    K(j,1)=x_chen(j);
    K(j,2)=y_chen(j);
    K(j,3)=z_chen(j);
end
T=h:h:TSim;
%

