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
h=figure;
view(3);
set(gca,'GridLineStyle','--')
hold on;
plot3(x,y,z,'LineWidth',1.5)
axis([min(x) max(x) min(y) max(y) min(z) max(z)]);
xlabel('x');
ylabel('y');
zlabel('z');
% initial radius
r0=[x(1) y(1) z(1)];
% final radius
rf=[x(end) y(end) z(end)];
plot3(r0(1),r0(2),r0(3),'or','MarkerFaceColor','r');
plot3(rf(1),rf(2),rf(3),'ok','MarkerFaceColor','k');
% Line of Sight between initial and final state
xx=linspace(r0(1),rf(1),10);
yy=linspace(r0(2),rf(2),10);
zz=linspace(r0(3),rf(3),10);
plot3(xx,yy,zz,'k--','LineWidth',2);
grid on;
% Distance
d=sqrt(sum((rf-r0).^2));
Information=strcat('Three dimensional Brownian Motion, d=',num2str(d),' units');
title(Information ,'FontWeight','bold');
view(-109,58);