%% Start simulation
clc;
time_start = tic;

%% Read image and resize it to 512x512
image1 = imread("Images/covid_pneumonia.jpeg");
image2 = imresize(image1, [512 512]);
% Creating a copy of the resized image
image3 = image2;

%% Henon Chaotic Map for Intra Block shuffling
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
    for m = 1:64
        a = 1+(m-1)*8;
        for i = b:b+7
            for j=a:a+7
                image2(i,j,:) = image3(indX(i),indY(j),:);
            end
        end
    end
end

%% Henon Chaotic Map for Inter Block shuffling
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

%% Brownian Motion Implementation
% This simulation illustrates a fast implementation of three dimensional
% Brownian motion, the output is the Euclidean distance between initial
% and final positions. (rf-r0)
N=1000;
t = nan( 1536256, 3 );
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
t(1,:) = rf-r0;
for i = 2:1:1536256
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
    t(i,:) = rf-r0;
end

%% Modulus,absolute functions and multiply operations
T1 = t(1:end,1);
T2 = t(1:end,2);
T3 = t(1:end,3);

T1(1:end-262144,:) = [];
T2(1:end-262144,:) = [];
T3(1:end-262144,:) = [];

U1 = T1(1:end);
U2 = T2(1:end);
U3 = T3(1:end);

V1 = U1 * 10^14;
V2 = U2 * 10^14;
V3 = U3 * 10^14;

W1 = reshape(V1,512,512);
W2 = reshape(V2,512,512);
W3 = reshape(V3,512,512);

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

Z1 = mod(Z1,256);
Z2 = mod(Z1,256);
Z3 = mod(Z1,256);

%% Chen Chaotic System
[simulation_time, K] = FOChen([35 3 28], [0.98 0.98 0.98], 20, [8 2 1]);

K(513:end,:) = [];
T = round(abs(K));

Z1_encrypt = bitor(Z1,T(:,1),'uint64');
Z2_encrypt = bitor(Z1,T(:,2),'uint64');
Z3_encrypt = bitor(Z1,T(:,3),'uint64');

%% Choose random direction, say Z3

image3(1:512,1:512,1) = Z3_encrypt;

converted_image = uint8(image3);
% imshow(converted_image3)
% figure
% imhist(converted_image3)
J = histeq(converted_image);
figure
subplot(1,2,1)
imshow(J)
subplot(1,2,2)
imhist(J,64)

%% End simulation
time_end = toc(time_start);
X = sprintf('Total time : %d',time_end);
disp(X);