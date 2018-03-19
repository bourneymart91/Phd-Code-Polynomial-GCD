addpath(genpath(pwd));



%Import image
origimage = imread('Images/ct','jpg');
origimage = imresize(origimage,[50 50])

%Reduce image to 2-D
origimage = origimage(:,:,1);

%Plot image
figure, imagesc(origimage)
axis square
colormap gray
title('Original Image')
set(gca, 'XTick', [], 'YTick', [])






%Blur Kernel
ksize = 20;
kernel = zeros(ksize);

%Gaussian Blur
s = 5;
m = ksize/2;
[X, Y] = meshgrid(1:ksize);
kernel = (1/(2*pi*s^2))*exp(-((X-m).^2 + (Y-m).^2)/(2*s^2));

%Display Kernel
figure, imagesc(kernel)
axis square
title('Blur Kernel')
colormap gray







%Embed kernel in image that is size of original image
[h, w] = size(origimage);
kernelimage = zeros(h,w);
kernelimage(1:ksize, 1:ksize) = kernel;

%Perform 2D FFTs
fftimage = fft2(double(origimage));
fftkernel = fft2(kernelimage);

%Set all zero values to minimum value
fftkernel(find(fftkernel == 0)) = 1e-6;

%Multiply FFTs
fftblurimage = fftimage.*fftkernel;

%Perform Inverse 2D FFT
blurimage = ifft2(fftblurimage);

%Display Blurred Image
figure, imagesc(blurimage)
axis square
title('Blurred Image')
colormap gray
set(gca, 'XTick', [], 'YTick', [])
 
 







%Pad image
origimagepad = padimage(origimage, ksize);

%Embed kernel in image that is size of original image + padding
[h1, w1] = size(origimagepad);
kernelimagepad = zeros(h1,w1);

kernelimagepad(1:ksize, 1:ksize) = kernel;

%Perform 2D FFTs
fftimagepad = fft2(origimagepad);
fftkernelpad = fft2(kernelimagepad);

fftkernelpad(find(fftkernelpad == 0)) = 1e-6;

%Multiply FFTs
fftblurimagepad = fftimagepad.*fftkernelpad;

%Perform Reverse 2D FFT
blurimagepad = ifft2(fftblurimagepad);

%Remove Padding
blurimageunpad = blurimagepad(ksize+1:ksize+h,ksize+1:ksize+w);

%Display Blurred Image
figure, imagesc(blurimageunpad)
axis square
title('Blurred Image - with Padding')
colormap gray
set(gca, 'XTick', [], 'YTick', [])




% Get two rows of pixel values
[nRows, nCols] = size(blurimagepad);

vec = randi(nRows,2,1);
row_index_1 = vec(1);
row_index_2 = vec(2);

fx = blurimagepad(row_index_1,:)';
gx = blurimagepad(row_index_2,:)';

vec = randi(nCols,2,1);
col_index_1 = vec(1);
col_index_2 = vec(2);

fy = blurimagepad(:, col_index_1);
gy = blurimagepad(:, col_index_1);

m = GetDegree(fx);
n = GetDegree(gx);

limits_t = [0,min(m,n)];
rank_range = [0,0];

global SETTINGS
SETTINGS.MEAN_METHOD = 'Geometric Mean Matlab Method';
SETTINGS.BOOL_ALPHA_THETA = true;
SETTINGS.EX_NUM = '1';
SETTINGS.METRIC = 'Minimum Singular Values';
SETTINGS.PLOT_GRAPHS = true;

dx = o_gcd_mymethod_Univariate_2Polys(fx, gx, limits_t, rank_range);

display(dx)
