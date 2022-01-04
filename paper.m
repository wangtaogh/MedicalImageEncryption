image1 = imread("8.gif");
imsize = size(image1);
%imshow(image1);
% Size of the image nxn
% n = 512;
image2 = imresize(image1, [512 512]);
imshow(image2)
%imshow(image2(1,3:end))
rows = [1:512];
columns = [1:512];

pixelIntensity = reshape(image2(:,:,1),512,512);

% For the Henon map to be chaotic, values of a and b : 
a = 1.4;
b = 0.3;

% secret key for Henon Chaotic Map 
% x0 = input("Secret Key for the Henon Chaotic Map in the X irection: ");
% y0 = input("Secret Key for the Henon Chaotic Map in the Y direction: ");
x0 = 0.00002;
y0 = 0.27;
x(1) = 1-a*x0^2+y0;
y(1) = b*x0;
for i = 2:512 
   x(i) = 1-a*x(i-1)^2+y(i-1);
   y(i) = b*x(i-1);
end
[outX,indX] = sort(x);
[outY,indY] = sort(y);
% for n = 1:64
%     b = 1+(n-1)*8;
%     for n = 1:64                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
%         a = 1+(n-1)*8;
%         for i = b:b+7
%             for j=a:a+7
%                 pixelIntensity(i,j) = pixelIntensity(indX(i),indY(j));
%                 image2(i,j,:) = image2(indX(i),indY(j));
%             end
%         end
%     end
% end
for n = 1:8
    b = 1+(n-1)*8;
    for n = 1:8                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
        a = 1+(n-1)*8;
        for i = b:b+63
            for j=a:a+63
                pixelIntensity(i,j) = pixelIntensity(i+(indX(n)-1)*64,j+(indY(n)-1)*64);
                image2(i,j,:) = image2(i+(indX(n)-1)*64,j+(indY(n)-1)*64);
            end
        end
    end
end

imshow(image2)



% This simulation illustrates a fast implementation of three dimensional
% Brownian motion, the output is the Euclidean distance between initial
% and final positions.
% To calculate the mean value of T runs, run the following code in the 
% Command window :
%
% >>T=100;
% >> for n=1:T
% >>Brownianmotion;
% >>close;
% >>D(n)=d;
% >>end
% >>figure; plot(D),title(' Distance ');
% >>dd=mean(D)
% (c) Youssef Khmou, Applied mathematics, may 2015.
N=1000;
x=cumsum(randn(1,N));
y=cumsum(randn(1,N));
z=cumsum(randn(1,N));
% define starting point of brownian particle - secret key 
x(1)=0;
y(1)=0;
z(1)=0;
r0=[x(1) y(1) z(1)];
% final radius
rf=[x(end) y(end) z(end)];
t = rf-r0;
for i = 2:1536256
    x=cumsum(randn(1,N));
    y=cumsum(randn(1,N));
    z=cumsum(randn(1,N));
    % define starting point of brownian particle - secret key 
    x(1)=0;
    y(1)=0;
    z(1)=0;
    % h=figure;
    % view(3);
    % set(gca,'GridLineStyle','--')
    % hold on;
    % plot3(x,y,z,'LineWidth',1.5)
    % axis([min(x) max(x) min(y) max(y) min(z) max(z)]);
    % xlabel('x');
    % ylabel('y');
    % zlabel('z');
    % initial radius
    m = r0;
    n = rf;
    r0=[x(1) y(1) z(1)];
    % final radius
    rf=[x(end) y(end) z(end)];
    t = cat(1,t,rf-r0);
end
t(1:1274112,:) = [];





% plot3(r0(1),r0(2),r0(3),'or','MarkerFaceColor','r');
% plot3(rf(1),rf(2),rf(3),'ok','MarkerFaceColor','k');
% Line of Sight between initial and final state
% xx=linspace(r0(1),rf(1),10);
% yy=linspace(r0(2),rf(2),10);
% zz=linspace(r0(3),rf(3),10);
% plot3(xx,yy,zz,'k--','LineWidth',2);
% grid on;
% Distance
% d=sqrt(sum((rf-r0).^2));
% Information=strcat('Three dimensional Brownian Motion, d=',num2str(d),' units');
% title(Information ,'FontWeight','bold');
% view(-109,58);








%%%%% 
% Chen Chaotic System
%%%%%
function [T, Y]=FOChen(parameters, orders, TSim, Y0)
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
c=parameters(3); d=parameters(4);
% binomial coefficients calculation:
cp1=1; cp2=1; cp3=1;
for j=1:n
    c1(j)=(1-(1+q1)/j)*cp1;
    c2(j)=(1-(1+q2)/j)*cp2;
    c3(j)=(1-(1+q3)/j)*cp3;
    cp1=c1(j); cp2=c2(j); cp3=c3(j);
end
% initial conditions setting:
x(1)=Y0(1); y(1)=Y0(2); z(1)=Y0(3);
% calculation of phase portraits /numerical solution/:
for i=2:n
    x(i)=(a*(y(i-1)-x(i-1)))*h^q1 - memo(x, c1, i);
    y(i)=(-d*x(i)-x(i)*z(i-1)+c*y(i-1))*h^q2 - memo(y, c2, i);
    z(i)=(x(i)*y(i)-b*z(i-1))*h^q3 - memo(z, c3, i);
end
for j=1:n
    Y(j,1)=x(j);
    Y(j,2)=y(j);
    Y(j,3)=z(j);
end
T=h:h:TSim;
end
