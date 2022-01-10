clc;
clear;
figure(2)

for i=1:10
    z_p=1;
    I = imread(['elements\lb_',num2str(z_p,'%2d'),'_',num2str(i,'%2d'),'.png']); 
    I=im2gray(I);
%     imshow(I)
%     title('Original Image')
    [m,n]=size(I);
    mask = false(size(I)); 
    mask(1,1) = true;
    mask(1,n) = true;
    mask(m,1) = true;
    mask(m,n) = true;
    W = graydiffweight(I, mask, 'GrayDifferenceCutoff', 25);
    thresh = 0.01;
    [BW, D] = imsegfmm(W, mask, thresh);
    subplot(2,5,i)
    imshow(BW)
%     title('Segmented Image')
    clear all;
end

figure(3)
for i=1:10
    z_p=1;
    I = imread(['elements\lb_',num2str(z_p,'%2d'),'_',num2str(i,'%2d'),'.png']); 
    subplot(2,5,i)
    I=im2gray(I);
    imshow(I)
end