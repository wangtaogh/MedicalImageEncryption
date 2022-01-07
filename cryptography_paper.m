image1 = imread("covid_pneumonia.jpeg");
imsize = size(image1);
image2 = imresize(image1, [512 512]);
rows = [1:512];
columns = [1:512];

pixelIntensity = reshape(image2(:,:,1),512,512);
image3 = image2;

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
y(513:end) = [];
[outX,indX] = sort(x);
[outY,indY] = sort(y);
for n = 1:64
    b = 1+(n-1)*8;
    for n = 1:64
        a = 1+(n-1)*8;
        for i = b:b+7
            for j=a:a+7
                pixelIntensity(i,j) = pixelIntensity(indX(i),indY(j));
                image2(i,j,:) = image3(indX(i),indY(j),:);
            end
        end
    end
end
% newIndX = empty(1,64);
% newIndY = empty(1,64);
newIndX = [];
newIndY = [];

corner_coordinates = [1:8:512];
for i = 1:1:512
    if ismember(indX(i), corner_coordinates) == 1
        newIndX = [newIndX indX(i)];
    end
    if ismember(indY(i), corner_coordinates) == 1
        newIndY = [newIndY indY(i)];
    end
end

image3 = image2;

for i = 1:1:64
    for j = 1:1:64
        for a = 1:1:7
            for b = 1:1:7
                image2(1+(j-1)*8+a,1+(i-1)*8+b,:) = image3(newIndX(i)+a,newIndY(j)+b,:);
            end
        end
    end
end

image3 = image2;
imshow(image3);

%% This simulation illustrates a fast implementation of three dimensional
% Brownian motion, the output is the Euclidean distance between initial
% and final positions. (rf-r0)

N=1000;
x_brown=cumsum(randn(1,N));
y_brown=cumsum(randn(1,N));
z_brown=cumsum(randn(1,N));
% Define starting point of the brownian particle - secret key
x_brown(1)=0;
y_brown(1)=0;
z_brown(1)=0;
r0=[x_brown(1) y_brown(1) z_brown(1)];
% Final radius
rf=[x_brown(end) y_brown(end) z_brown(end)];
t = rf-r0;
for i = 2:5:1536256
    disp(i);
    x_brown=cumsum(randn(1,N));
    y_brown=cumsum(randn(1,N));
    z_brown=cumsum(randn(1,N));
    % Define starting point of the brownian particle - secret key
    x_brown(1)=0;
    y_brown(1)=0;
    z_brown(1)=0;
    % Storing r0 and rf values of previous iteration in m and n
    m = r0;
    n = rf;
    r0=[x_brown(1) y_brown(1) z_brown(1)];
    % Final radius
    rf=[x_brown(end) y_brown(end) z_brown(end)];
    t = cat(1,t,rf-r0);
end
%%
T1 = t(1:end,1);
T2 = t(1:end,2);
T3 = t(1:end,3);

T1(1:end-262144,:) = [];
T2(1:end-262144,:) = [];
T3(1:end-262144,:) = [];


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

image3 = double(image3);
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
q1 = 0.997;
q2 = 0.1666;
q3 = 0.4765;

% Binomial coefficients calculation:
cp1=1; cp2=1; cp3=1;
for j=1:1:512
    c1(j)=(1-(1+q1)/j)*cp1;
    c2(j)=(1-(1+q2)/j)*cp2;
    c3(j)=(1-(1+q3)/j)*cp3;
    cp1=c1(j); cp2=c2(j); cp3=c3(j);
end
% initial conditions setting:
% x_chen(1)=Y0(1); y_chen(1)=Y0(2); z_chen(1)=Y0(3);
x_chen(1)=1; y_chen(1)=3; z_chen(1)=4;
% calculation of phase portraits /numerical solution/:
for i=2:1:512
    x_chen(i)=(a*(y_chen(i-1)-x_chen(i-1)))*h^q1 - memo(x_chen, c1, i);
    y_chen(i)=(-d*x_chen(i)-x_chen(i)*z_chen(i-1)+c*y_chen(i-1))*h^q2 - memo(y_chen, c2, i);
    z_chen(i)=(x_chen(i)*y_chen(i)-b*z_chen(i-1))*h^q3 - memo(z_chen, c3, i);
end
j=1:1:512;
K(j,1)=x_chen(j);
K(j,2)=y_chen(j);
K(j,3)=z_chen(j);

T = round(abs(K));

Z1_encrypt = [];
Z2_encrypt = [];
Z3_encrypt = [];
for i = 1:512
    for j = 1:1:512
        disp(i)
        Z1_encrypt = [Z1_encrypt bitor(Z1(j,i),T(j,1),'uint64')];
        Z2_encrypt = [Z2_encrypt bitor(Z2(j,i),T(j,2),'uint64')];
        Z3_encrypt = [Z3_encrypt bitor(Z3(j,i),T(j,3),'uint64')];
    end
end
% [i,j] = meshgrid(1:512,1:512);
% disp(i)
% Z1_encrypt = bitor(Z1(j,i),T(j,1),'uint64');
% Z2_encrypt = bitor(Z2(j,i),T(j,2),'uint64');
% Z3_encrypt = bitor(Z3(j,i),T(j,3),'uint64');
%%
Z1_encrypt = transpose(mod(Z1_encrypt,256));
Z2_encrypt = transpose(mod(Z2_encrypt,256));
Z3_encrypt = transpose(mod(Z3_encrypt,256));

% Choose random direction, say Z3

Z3_reshaped = reshape(Z3_encrypt,512,512);
image3(1:512,1:512,1) = Z3_reshaped;

converted_image3 = uint8(image3);
imshow(converted_image3)
% figure
% imhist(converted_image3)
J = histeq(converted_image3);
figure
subplot(1,2,1)
imshow(J)
subplot(1,2,2)
imhist(J,64)