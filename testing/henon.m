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
% for n = 1:64
%     a = 1+(n-1)*8;
%     for m = 1:64                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
%         b = 1+(m-1)*8;
%         for i = a:a+7
%             for j=b:b+7
%                 image2(i,j,:) = image2(i+(indX(a)-1),j+(indY(b)-1),:);
%             end
%         end
%     end
% end
% Inter block shuffling 
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
imshow(image2)
