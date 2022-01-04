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



