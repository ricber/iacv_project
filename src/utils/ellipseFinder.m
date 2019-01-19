function ellipse = ellipseFinder(im_gray, MIN, MAX, cannyLow, cannyHigh)
%ellipseFinder This function detects an ellipse given an image
%   This function uses the function ellipseDetection, that is, a script
%   written by Martin Simonovsky based on the paper "A New Efficient Ellipse Detection Method" 
%   (Yonghong Xie Qiang , Qiang Ji / 2002). The function take as input a
%   gray scale image and opens the crop tool to reduce the computation
%   load. It also takes as input the minimum and maximum length of the major axis
%   to speed up the detection of the ellipse and the parametrs for Canny edge 
%   detection. The user have to manually contour 
%   the ellipse with a rectangle as smaller as possible. The function returns
%   the position of the ellipse center
%   in the original image, the length of the semi-major and semi-minor axes and the
%   rotation angle. 

debug_print = false;
debug_im = false;

% crop the image because the algorithm to find ellipses is memory intensive
figure;
[cropped_image, rect] = imcrop(im_gray); %rect is the size and position of the cropped rectangle 
col_crop = rect(1,1); % x of cropped rectangle
row_crop = rect(1,2); % y of cropped rectangle

if debug_print
   fprintf("Cropped image. x: %d; y: %d; width: %d; height: %d", col_crop, row_crop, rect(1,3), rect(1,4)); 
end

% imtool(im_gray); % to measure the length of the ellipse major axis in pixels
 
% find edges of the cropped image with canny 
edges = edge(cropped_image, 'canny', [cannyLow cannyHigh]);
imshow(edges);

% override some default parameters
params.minMajorAxis = MIN;
params.maxMajorAxis = MAX;
params.numBest = 1; % return the best one
params.randomize = 0; % turn off randomization to obtain better results

% note that the edge (or gradient) image is used
bestFits = ellipseDetection(edges, params);

x_c_ellipse = bestFits(1,1) + col_crop;
y_c_ellipse = bestFits(1,2) + row_crop;
semiMajorAxis = bestFits(1,3);
semiMinorAxis = bestFits(1,4);
theta = bestFits(1,5);
ellipse = [x_c_ellipse, y_c_ellipse, semiMajorAxis, semiMinorAxis, theta];

if debug_print
    fprintf('Output %d best fits.\n', size(bestFits,1));
end

if debug_im
   ellipsePlot(bestFits(:,3),bestFits(:,4),bestFits(:,5)*pi/180,bestFits(:,1),bestFits(:,2),'r');
   figure, imshow(cropped_image);
   % ellipse drawing implementation: http://www.mathworks.com/matlabcentral/fileexchange/289 
   ellipsePlot(bestFits(:,3),bestFits(:,4),bestFits(:,5)*pi/180,bestFits(:,1),bestFits(:,2),'r');
end

