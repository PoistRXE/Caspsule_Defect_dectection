function BW=background_Area(I)
    I=im2gray(I);
%     imshow(I)
%     title('Original Image')
    [m,n]=size(I);
    mask = false(size(I)); 
    mask(1,:) = true;
    mask(1,n) = true;
    mask(m,:) = true;
    mask(m,n) = true;
    W = graydiffweight(I, mask, 'GrayDifferenceCutoff', 25);
    thresh = 0.01;
    [BW, ~] = imsegfmm(W, mask, thresh);
end
