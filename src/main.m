%% Initializations

close all
clear
clc

image = 1; % choose the input image, 1 or 2 (THE PROGRAM CURRENTLY WORKS ONLY WITH IMAGE1)
debug_im = false; % debug mode image: shows debug images
debug_print = false; % debug mode print: shows debug strings
interactive = false;    % set this variable to true if you want to manually 
                        % select points, otherwise the hardcoded ones will
                        % be chosen


%% Load the image

fprintf("Loading the image...\n");

if image == 1
    im_rgb = imread('Image1.jpg');
else
    im_rgb = imread('Image2.jpg');
end

im_gray = rgb2gray(im_rgb); % a lot of algorithms work on gray scale

% the program is currently not working
% properly with the nornalization enabled
normalization = false;
if normalization
    norm_matrix = diag([1/size(im_rgb,2), 1/size(im_rgb,1), 1]);
    
    % An alternative normalization. BEWARE: hardcoded components are NOT
    % normalized with the following matrix. 
    % norm_matrix = diag([1/max(size(im_rgb)), 1/max(size(im_rgb)), 1]); 
else
    norm_matrix = diag([1, 1, 1]); % normalization off
end

%% Keypoints detection & selection

fprintf("Finding keypoint features...\n");

[xi, yi] = keypointsSelection(interactive, image, im_gray, 'harris');


%% Find ellipses

fprintf("Finding ellipses...\n");

[C1, H1, C2, H2] = ellipsesSelection(interactive, image, im_gray);

H1 =  norm_matrix*H1; % normalized rototranslation
H2 =  norm_matrix*H2; % normalized rototranslation


%% Lines tangent to ellipses

fprintf("Finding tangency points of lines with the two ellipses...\n");

ellipsesPoints = ellipsesPointsSelection(interactive, normalization, C1, H1, C2, H2);

% tangency points of parallel lines
npt_tng_ul = ellipsesPoints(1,:)'; % upper (u) normalized (n) point (pt) tangent (tng) to the left (l) ellipse
npt_tng_ur = ellipsesPoints(2,:)';
npt_tng_dr = ellipsesPoints(3,:)';
npt_tng_dl = ellipsesPoints(4,:)';

pt_tng_ul = norm_matrix\npt_tng_ul; % upper (u) NOT-normalized point (pt) tangent (tng) to the left (l) ellipse
pt_tng_ur = norm_matrix\npt_tng_ur;
pt_tng_dr = norm_matrix\npt_tng_dr;
pt_tng_dl = norm_matrix\npt_tng_dl;

% tangency points of crossing lines
npt_cross_ul = ellipsesPoints(5,:)';
npt_cross_ur = ellipsesPoints(6,:)';
npt_cross_dr = ellipsesPoints(7,:)';
npt_cross_dl = ellipsesPoints(8,:)';

% tangency lines
nline_tng_u = cross(npt_tng_ul, npt_tng_ur); % upper (u) normalized (n) tangent (tng) line
nline_tng_d = cross(npt_tng_dl, npt_tng_dr);

%% Polar lines of the ellipses and vanishing points on that plane

fprintf("Finding polar lines...\n");

% vertical polar lines
nline_v_polar_l = cross(npt_tng_ul, npt_tng_dl);
nline_v_polar_r = cross(npt_tng_ur, npt_tng_dr);

% vanishing points
nh_inf = cross(nline_tng_u, nline_tng_d); % horizontal vanishing point
nv_inf = cross(nline_v_polar_l, nline_v_polar_r); % vertical vanishing point

% LEFT polar
syms x1 x2
x = [x1; x2; 1];
eqns = [x.'*(inv(H1).'*C1*inv(H1)).'*nv_inf == 0, x.'*(inv(H1).'*C1*inv(H1)).'*x == 0];
vars = [x1 x2];
[sol_x1, sol_x2] = solve(eqns, vars);

% I want to know what's the left and right point of the polar
if (double(sol_x1(1)) < double(sol_x1(2)))
    npt_h_polar_ll = [double(sol_x1(1)); double(sol_x2(1)); 1];
    npt_h_polar_lr = [double(sol_x1(2)); double(sol_x2(2)); 1];
else
    npt_h_polar_ll = [double(sol_x1(2)); double(sol_x2(2)); 1];
    npt_h_polar_lr = [double(sol_x1(1)); double(sol_x2(1)); 1];
end

nline_h_polar_l = cross(npt_h_polar_ll, npt_h_polar_lr);

% RIGHT polar
syms x1 x2
x = [x1; x2; 1];
eqns = [x.'*(inv(H2).'*C2*inv(H2)).'*nv_inf == 0, x.'*(inv(H2).'*C2*inv(H2)).'*x == 0];
vars = [x1 x2];
[sol_x1, sol_x2] = solve(eqns, vars);

% I want to know what's the left and right point of the polar
if (double(sol_x1(1)) < double(sol_x1(2)))
    npt_h_polar_rl = [double(sol_x1(1)); double(sol_x2(1)); 1];
    npt_h_polar_rr = [double(sol_x1(2)); double(sol_x2(2)); 1];
else
    npt_h_polar_rl = [double(sol_x1(2)); double(sol_x2(2)); 1];
    npt_h_polar_rr = [double(sol_x1(1)); double(sol_x2(1)); 1];
end

nline_h_polar_r = cross(npt_h_polar_rl, npt_h_polar_rr);

% images of circumferences centers
ncenter_l = cross(nline_h_polar_l, nline_v_polar_l);
ncenter_r = cross(nline_h_polar_r, nline_v_polar_r);

center_l = norm_matrix\ncenter_l;
center_r = norm_matrix\ncenter_r;

nline_centers = cross(ncenter_l, ncenter_r);

% Here I fit the horizontal vanishing points to three lines, that is, the
% two bitangent lines and the line passing through the centers, so to
% obtain a better accuracy 
fnh_inf = fitVp([nline_tng_u, nline_tng_d, nline_centers]); % fitted (f) normalized (n) horizontal (h) vanishing point

fprintf("Drawing lines...\n");

% lines drawings
figure(1), imshow(im_gray);

lines_wheels = [(pt_tng_ul(1:2)/pt_tng_ul(3))'   (pt_tng_ur(1:2)/pt_tng_ur(3))';
                (pt_tng_dl(1:2)/pt_tng_dl(3))'   (pt_tng_dr(1:2)/pt_tng_dr(3))';
                (center_l(1:2)/center_l(3))'     (center_r(1:2)/center_r(3))'  ];          
linePlot(lines_wheels);
hold on;

pt_h_polar_ll = norm_matrix\npt_h_polar_ll;
pt_h_polar_lr = norm_matrix\npt_h_polar_lr;
pt_h_polar_rl = norm_matrix\npt_h_polar_rl;
pt_h_polar_rr = norm_matrix\npt_h_polar_rr;

lines_polar = [ (pt_tng_ul(1:2)/pt_tng_ul(3))'          (pt_tng_dl(1:2)/pt_tng_dl(3))';
                (pt_tng_ur(1:2)/pt_tng_ur(3))'          (pt_tng_dr(1:2)/pt_tng_dr(3))'];
linePlot(lines_polar, 'red');        
lines_polar = [ (pt_h_polar_ll(1:2)/pt_h_polar_ll(3))'  (pt_h_polar_lr(1:2)/pt_h_polar_lr(3))';
                (pt_h_polar_rl(1:2)/pt_h_polar_rl(3))'  (pt_h_polar_rr(1:2)/pt_h_polar_rr(3))']; 
linePlot(lines_polar, 'blue');


if debug_im
    if image == 1
        figure(1)
        pt_cross_ul = norm_matrix\npt_cross_ul; % upper (u) NOT-normalized point (pt) tangent (tng) to the left (l) ellipse
        pt_cross_ur = norm_matrix\npt_cross_ur;
        pt_cross_dr = norm_matrix\npt_cross_dr;
        pt_cross_dl = norm_matrix\npt_cross_dl;
       
        plot(pt_tng_ul(1)/pt_tng_ul(3), pt_tng_ul(2)/pt_tng_ul(3),'x');
        plot(pt_tng_ur(1)/pt_tng_ur(3), pt_tng_ur(2)/pt_tng_ur(3),'o');
        plot(pt_tng_dr(1)/pt_tng_dr(3), pt_tng_dr(2)/pt_tng_dr(3),'d');
        plot(pt_tng_dl(1)/pt_tng_dl(3), pt_tng_dl(2)/pt_tng_dl(3),'s');
    
        
        plot(pt_cross_ul(1)/pt_cross_ul(3), pt_cross_ul(2)/pt_cross_ul(3),'v');
        plot(pt_cross_ur(1)/pt_cross_ur(3), pt_cross_ur(2)/pt_cross_ur(3),'p');
        plot(pt_cross_dr(1)/pt_cross_dr(3), pt_cross_dr(2)/pt_cross_dr(3),'>');
        plot(pt_cross_dl(1)/pt_cross_dl(3), pt_cross_dl(2)/pt_cross_dl(3),'h');
        
      
        h_inf = norm_matrix\nh_inf;
        v_inf = norm_matrix\nv_inf;
        plot(h_inf(1)/h_inf(3), h_inf(2)/h_inf(3),'x');
        plot(v_inf(1)/v_inf(3), v_inf(2)/v_inf(3),'x');
        
        
        plot(pt_h_polar_ll(1)/pt_h_polar_ll(3), pt_h_polar_ll(2)/pt_h_polar_ll(3),'<');
        plot(pt_h_polar_lr(1)/pt_h_polar_lr(3), pt_h_polar_lr(2)/pt_h_polar_lr(3),'^');
        
        plot(pt_h_polar_rl(1)/pt_h_polar_rl(3), pt_h_polar_rl(2)/pt_h_polar_rl(3),'+');
        plot(pt_h_polar_rr(1)/pt_h_polar_rr(3), pt_h_polar_rr(2)/pt_h_polar_rr(3),'*');
    end
end
        
        
%% Wheels diameter and wheel-to-wheel distance ratio

fprintf("Finding wheels diameter and wheel-to-wheel distance ratio...\n");

% here I find the length of the wheels rays
ray_l = norm(pt_h_polar_ll(1:2)/pt_h_polar_ll(3) - center_l(1:2)/center_l(3));
ray_r = norm(center_r(1:2)/center_r(3) - pt_h_polar_rr(1:2)/pt_h_polar_rr(3));
% this is the distance between the centers of the wheels 
wheelbase = norm(center_l(1:2)/center_l(3) - center_r(1:2)/center_r(3));

% cross ratio
CR = (ray_l / (ray_l + wheelbase)) / ((ray_r + wheelbase) / ray_r); 

fprintf("The cross-ratio is: %f\n", CR); 

% ratio = (wheel_ray)/(wheelbase)
syms ratio 
eqn = (1/ratio)^2*CR + (1/ratio)*2*CR + (CR-1) == 0;
[sol_ratio] = solve(eqn);

diam_wheelbase_ratio = double(sol_ratio(1))*2;

fprintf("The wheel diameter and wheelbase ratio is: %f\n", diam_wheelbase_ratio);

%% Camera calibration

fprintf("Calibrating camera...\n");

% these are four symmetric keypoints of the car rear part
% they are used to compute a vanishing point
npt_rear_ul = norm_matrix*[xi(1); yi(1); 1];
npt_rear_ur = norm_matrix*[xi(2); yi(2); 1];
npt_rear_dr = norm_matrix*[xi(3); yi(3); 1];
npt_rear_dl = norm_matrix*[xi(4); yi(4); 1];

nline_rear_high = cross(npt_rear_ul, npt_rear_ur); % higher horizontal rear line
nline_rear_low = cross(npt_rear_dl, npt_rear_dr); % lower horizontal rear line

% vanishing point of the horizontal rear direction
nh_rear_inf = cross(nline_rear_high, nline_rear_low);

% rear vanishing lines plot
figure(1);
lines_rear_h = [    xi(1) yi(1)  xi(2) yi(2);
                    xi(3) yi(3)  xi(4) yi(4)];
linePlot(lines_rear_h, 'magenta');

% image of the line at the infinity of the plane containing the wheels
% circumferences
nl_inf_im = cross(fnh_inf, nv_inf);

% images of the circular points
syms x1 x2
x = [x1; x2; 1];
% TRY: additional constraint to increase robustness (second ellipse): 
% x.'*(inv(norm_matrix\H2).'*C2*inv(norm_matrix\H2))*x == 0];
eqns = [    nl_inf_im.'*x == 0, ... 
            x.'*(inv(H1).'*C1*inv(H1))*x == 0];
vars = [x1 x2];
[sol_x1, sol_x2] = solve(eqns, vars);

I_im = [double(sol_x1(1)); double(sol_x2(1)); 1];
J_im = [double(sol_x1(2)); double(sol_x2(2)); 1];

% w is a representation of the matrix of the IAC (image of absolute conic)
% as homogeneous 6-vector
syms w1 w4 w5 w6
w_sym = [   w1  0   w4;
            0   1   w5;
            w4  w5  w6];

% fitting of the IAC with orthogonal vanishing points and image of the
% circular points
eqns = [nv_inf'*w_sym*fnh_inf == 0, ...
        nv_inf'*w_sym*nh_rear_inf == 0, ... 
        fnh_inf'*w_sym*nh_rear_inf == 0, ...
        real(I_im)'*w_sym*imag(I_im) == 0, ...
        real(I_im)'*w_sym*real(I_im) - imag(I_im)'*w_sym*imag(I_im) == 0];

[A,b] = equationsToMatrix(eqns);

% fit a linear model without intercept
lm = fitlm(double(A),double(b), 'y ~ x1 + x2 + x3 + x4 - 1');
% get the coefficients
W = lm.Coefficients.Estimate;

% image of absolute conic
w = double([W(1,1) 0 W(2,1); 0 1 W(3,1); W(2,1) W(3,1) W(4,1)]);
w = norm_matrix*w*norm_matrix; % reverting the normalization

a = sqrt(w(1, 1));
u = -w(1, 3)/w(1, 1);
v = -w(2, 3);
fy = sqrt(w(3, 3) - w(1, 1)*u^2 - v^2);
fx = fy/a;

K = [   fx,    0,      u;
        0,     fy,     v;
        0,     0,      1];
    
%% Reconstruction of the 3d position of symmetric features

fprintf("Reconstructing 3D position of symmetric features...\n");

% 1st STEP
% Find the rotation matrix from the camera reference to the world
% reference R_c_w. This matrix can put points of the world in the reference
% frame of the camera just multiplying it before the point
h_rear_inf = norm_matrix\nh_rear_inf;
i_pi = K\h_rear_inf; % the direction of the x axis is the one of the rear vanishing point
i_pi = i_pi/norm(i_pi);
fh_inf = norm_matrix\fnh_inf;
k_pi = K\fh_inf ; % the direction of the z axis is the one of the horizontal vanishing point
k_pi = k_pi/norm(k_pi);
j_pi = cross(k_pi, i_pi); % the direction of the y axis
j_pi = j_pi/norm(j_pi);
R_c_w = [i_pi, j_pi, k_pi];



% 2nd STEP
% Find the translation vector
% I want to find the middle point of the two symmetric features
pt_rear_dl = norm_matrix\npt_rear_dl;
pt_rear_dr = norm_matrix\npt_rear_dr;
bp_pt_rear_dl = (K*R_c_w) \ pt_rear_dl; % back projection (bp) of point pt_rear_dl
bp_pt_rear_dr = (K*R_c_w) \ pt_rear_dr;

k = 4; % fixed distance of symmetric features 

v1 = bp_pt_rear_dr/norm(bp_pt_rear_dr); % directions of the viewing ray from the point pt_rear_dr
v2 = bp_pt_rear_dl/norm(bp_pt_rear_dl);

% here I find the symmetric features
syms y z dist_l dist_r
eqns = [dist_l*v1(1,1) == k/2, ...
        dist_l*v1(2,1) == y, ...
        dist_l*v1(3,1) == z, ... 
        dist_r*v2(1,1) == k/2, ...
        dist_r*v2(2,1) == y, ...
        dist_r*v2(3,1) == z];

[A,b] = equationsToMatrix(eqns);

A = double(A);
b = double(b);
sols = (A' * A) \ (A' * b);

% 3D points reconstructed
pt_3d_rear_dl = [k/2; sols(1); sols(2)];
pt_3d_rear_dr = [-k/2; sols(1); sols(2)];
pt_3d_rear_mid= [0; sols(1); sols(2)];

% 3D points and reference frame plotting
global ColorOrder, ColorOrder='rgb';
starts = [pt_3d_rear_mid'; pt_3d_rear_mid'; pt_3d_rear_mid'; zeros(3,3)];
ends = [(i_pi+pt_3d_rear_mid)'; (j_pi+pt_3d_rear_mid)'; (k_pi+pt_3d_rear_mid)'; [1 0 0]; [0 1 0]; [0 0 1]];

figure;
arrow3(starts, ends, 'o');
axis equal
hold on;
grid on;
xlabel('X')
ylabel('Y')
plot3(pt_3d_rear_dr(1,1), pt_3d_rear_dr(2,1), pt_3d_rear_dr(3,1), '*');
plot3(pt_3d_rear_dl(1,1), pt_3d_rear_dl(2,1), pt_3d_rear_dl(3,1), '*');
plotCamera('Location', [0 0 0],'Orientation', eye(3),'Size',0.1);




