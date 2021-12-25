image1 = imread("8.gif");
imsize = size(image1);
%imshow(image1);
image2 = imresize(image1, [512 512]);
imshow(image2);
%imshow(image2(1,3:end))
rows = [1:512];
columns = [1:512];
pixel_values = impixel(image2,rows,columns)
% For a grayscale image, values of R,G,B will be the same
pixelIntensity = pixel_values(:,1);

% secret key for Henon Chaotic Map 
x0 = input("Secret Key for the Henon Chaotic Map : ");
x2 = 1-a*x^2 