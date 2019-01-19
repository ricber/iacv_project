function [C1, H1, C2, H2] = ellipsesSelection(interactive, image_sel, im_gray)
    % this function computes the ellipses real-time if interaction is true
    % otherwise it returns hardcoded ellipses
    
    if interactive
     
        % selection of Image1 or Image2
        if image_sel == 1
            % input parameters: gray scale image, min length major axis, max length mjor axis, Canny low thresh, Canny high thresh
            ellipse1 = ellipseFinder(im_gray, 455, 460, 0.05, 0.3); % left ellipse
            ellipse2 = ellipseFinder(im_gray, 285, 290, 0.09, 0.1); % right ellipse (backup: 0.001 0.3)
        else
            ellipse1 = ellipseFinder(im_gray, 633, 635, 0.07, 0.071); 
            ellipse2 = ellipseFinder(im_gray, 390, 393, 0.1, 0.16); 
        end

        ellipses = [ellipse1; ellipse2];
        
        figure(2), imshow(im_gray);
        ellipsePlot(ellipses(:,3), ellipses(:,4), ellipses(:,5)*pi/180, ellipses(:,1), ellipses(:,2),'r');

        % conics (ellipses) matrices
        x_c_1 = ellipse1(1,1); % x position of the center
        y_c_1 = ellipse1(1,2); % y position of the center
        semi_maj_1 = ellipse1(1,3); % semi major axis
        semi_min_1 = ellipse1(1,4); % semi minor axis
        theta1 = ellipse1(1,5)*pi/180; % angle

        C1 = [  1/(semi_maj_1)^2        0                       0;
                0                       1/(semi_min_1)^2        0;
                0                       0                       -1
             ];

        % rototranslation matrix
        R1 = [  cos(theta1)     -sin(theta1)    0;
                sin(theta1)     cos(theta1)     0;
                0               0               1]; % rotation
        T1 = [  1       0       x_c_1;
                0       1       y_c_1;
                0       0       1    ]; % translation
        H1 = T1*R1; % rototranslation

        x_c_2 = ellipse2(1,1);
        y_c_2 = ellipse2(1,2);
        semi_maj_2 = ellipse2(1,3); % semi major axis
        semi_min_2 = ellipse2(1,4); % semi minor axis
        theta2 = ellipse2(1,5)*pi/180; % angle

        C2 = [  1/(semi_maj_2)^2        0                       0;
                0                       1/(semi_min_2)^2        0;
                0                       0                       -1];

        % rototranslation matrix
        R2 = [  cos(theta2)     -sin(theta2)    0;
                sin(theta2)     cos(theta2)     0;
                0               0               1];
        T2 = [  1       0       x_c_2;
                0       1       y_c_2;
                0       0       1    ];
        H2 = T2*R2;
        
    else
        if image_sel == 1
            % HARDCODED ELLIPSES
            C1 = [  1.90532447957627e-05    0                       0;
                    0                       7.97193877551020e-05    0;
                    0                       0                       -1];
            H1 = [  0.255352600253544       0.966847997124550   2299.01000000000;
                    -0.966847997124550      0.255352600253544   1690.01000000000;
                    0                       0                   1               ];
            C2 = [  4.75861889698207e-05    0                       0;
                    0                       0.000540832882639265    0;
                    0                       0                       -1];
            H2 = [  0.217295721641613       0.976105818728815   3229.01000000000;
                    -0.976105818728815      0.217295721641613   1108.01000000000;
                    0                       0                   1               ];
        else
            error("NO hardcoded ellipses for this Image");
        end
    end

end