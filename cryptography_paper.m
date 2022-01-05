image1 = imread("8.gif");
imsize = size(image1);
image2 = imresize(image1, [512 512]);
image2 = double(image2);
rows = [1:512];
columns = [1:512];

pixelIntensity = reshape(image2(:,:,1),512,512);

% For the Henon map to be chaotic, values of a and b : 
a = 1.4;
b = 0.3;

% Secret key for Henon Chaotic Map 
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
for n = 1:64
    b = 1+(n-1)*8;
    for n = 1:64                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
        a = 1+(n-1)*8;
        for i = b:b+7
            for j=a:a+7
                pixelIntensity(i,j) = pixelIntensity(indX(i),indY(j));
                image2(i,j,:) = image2(indX(i),indY(j),:);
            end
        end
    end
end
for m = 1:64
    b = 1+(m-1)*8;
    for n = 1:32
        a = 1+(n-1)*8;
        q1 = 1;
        q2 = 1;
        for i = a:a+7
            for j=b:b+7
            z = image2(i,j,:);
            image2(i,j,:) = image2(256+i-1,j,:);
            image2(256+i-1,j,:) = z;
            end
        end
    end
end
for m = 1:64
    b = 1+(m-1)*8;
    for n = 1:32
        a = 1+(n-1)*8;
        for i = b:b+7
            for j=a:a+7
            z = image2(i,j,:);
            image2(i,j,:) = image2(i,256+j-1,:);
            image2(i,256+j-1,:) = z;
            end
        end
    end
end
image3 = image2;


%% This simulation illustrates a fast implementation of three dimensional
% Brownian motion, the output is the Euclidean distance between initial
% and final positions. (rf-r0)

N=1000;
x=cumsum(randn(1,N));
y=cumsum(randn(1,N));
z=cumsum(randn(1,N));
%% Define starting point of the brownian particle - secret key  
x(1)=0;
y(1)=0;
z(1)=0;
r0=[x(1) y(1) z(1)];
% Final radius
rf=[x(end) y(end) z(end)];
t = rf-r0;
for i = 2:1536256
    disp(i);
    x=cumsum(randn(1,N));
    y=cumsum(randn(1,N));
    z=cumsum(randn(1,N));
    %% Define starting point of the brownian particle - secret key 
    x(1)=0;
    y(1)=0;
    z(1)=0;
    %% Storing r0 and rf values of previous iteration in m and n
    m = r0;
    n = rf;
    r0=[x(1) y(1) z(1)];
    % Final radius
    rf=[x(end) y(end) z(end)];
    t = cat(1,t,rf-r0);
end

T1 = t(1:end,1);
T2 = t(1:end,2);
T3 = t(1:end,3);

T1(1:1274112,:) = [];
T2(1:1274112,:) = [];
T3(1:1274112,:) = [];

U1 = T1(1:end);
U2 = T2(1:end);
U3 = T3(1:end);

W1 = reshape(U1,512,512);
W2 = reshape(U2,512,512);
W3 = reshape(U3,512,512);
%% Section 2
X1 = round(abs(W1));
X2 = round(abs(W2));
X3 = round(abs(W3));

Y1 = mod(X1,256);
Y2 = mod(X2,256);
Y3 = mod(X3,256);

Z1 = Y1*image3;
Z2 = Y2*image3;
Z3 = Y3*image3;

%%%%% 
% Chen Chaotic System
%%%%%
%
% Numerical Solution of the Fractional-Order Chen's System
%   D^q1 x(t) = a(y(t)-x(t))
%   D^q2 y(t) = dx(t) - x(t)z(t) + cy(t)
%   D^q3 z(t) = x(t)y(t) - bz(t)
% function [T, Y] = FOChen(parameters, orders, TSim, Y0)
% Input:    parameters - model parameters [a, b, c, d]
%           orders - derivatives orders [q1, q2, q3]
%           TSim - simulation time (0 - TSim) in sec
%           Y0 - initial conditions [Y0(1), Y0(2), Y0(3)]
% Output:   T - simulation time (0 : Tstep : TSim)
%           Y - solution of the system (x=Y(1), y=Y(2), z=Y(3))

% time step:
h=0.005; 
% orders of derivatives, respectively:
% q1=orders(1); q2=orders(2); q3=orders(3);
% constants of Chen's system:
% a=parameters(1); b=parameters(2); 
% c=parameters(3); d=parameters(4);

% Initialising constants for now :
a = 1;
b = 2;
c = 3;
d = a-c;
q1 = 0.5;
q2 = 0.3;
q3 = 0.4;

% Binomial coefficients calculation:
cp1=1; cp2=1; cp3=1;
for j=1:512
    c1(j)=(1-(1+q1)/j)*cp1;
    c2(j)=(1-(1+q2)/j)*cp2;
    c3(j)=(1-(1+q3)/j)*cp3;
    cp1=c1(j); cp2=c2(j); cp3=c3(j);
end
% initial conditions setting:
% x(1)=Y0(1); y(1)=Y0(2); z(1)=Y0(3);
x(1)=1; y(1)=3; z(1)=4; 
% calculation of phase portraits /numerical solution/:
for i=2:512
    x(i)=(a*(y(i-1)-x(i-1)))*h^q1 - memo(x, c1, i);
    y(i)=(-d*x(i)-x(i)*z(i-1)+c*y(i-1))*h^q2 - memo(y, c2, i);
    z(i)=(x(i)*y(i)-b*z(i-1))*h^q3 - memo(z, c3, i);
end
for j=1:512
    K(j,1)=x(j);
    K(j,2)=y(j);
    K(j,3)=z(j);
end
T = round(abs(K));

Z1_encrypt = [];
Z2_encrypt = [];
Z3_encrypt = [];
for i = 1:512
    for j = 1:512
        disp(i)
        Z1_encrypt = [Z1_encrypt bitor(Z1(j,i),T(j,1),'uint64')];
        Z2_encrypt = [Z2_encrypt bitor(Z2(j,i),T(j,2),'uint64')];
        Z3_encrypt = [Z3_encrypt bitor(Z3(j,i),T(j,3),'uint64')];
    end
end

Z1_encrypt = transpose(mod(Z1_encrypt,256));
Z2_encrypt = transpose(mod(Z2_encrypt,256));
Z3_encrypt = transpose(mod(Z3_encrypt,256));

% Choose random direction, say Z2 

Z2_reshaped = reshape(Z2_encrypt,512,512);
image3(1:512,1:512,:) = Z2_reshaped;
imshow(image3)
