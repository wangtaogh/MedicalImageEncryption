
I = magic(15); % replace with image2
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
% for i = 2:1536256
    for i = 2:256

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

T1 = t(1:end,1);
T2 = t(1:end,2);
T3 = t(1:end,3);
% T1(1:1274112,:) = [];
% T2(1:1274112,:) = [];
% T3(1:1274112,:) = [];

T1(1:31) = [];
T2(1:31) = [];
T3(1:31) = [];

U1 = T1(1:end);
U2 = T2(1:end);
U3 = T3(1:end);

W1 = reshape(U1,15,15);
W2 = reshape(U2,15,15);
W3 = reshape(U3,15,15);

X1 = round(abs(W1))*10^14;
X2 = round(abs(W2))*10^14;
X3 = round(abs(W3))*10^14;

Y1 = mod(X1,256);
Y2 = mod(X2,256);
Y3 = mod(X3,256);

Z1 = Y1*I;
Z2 = Y2*I;
Z3 = Y3*I;

