function [xi, yi] = keypointsSelection(interactive, image_sel, im_gray, method)
% this function finds the keypoints based on the method provided as
% parameter and allows the user to select the keypoints he desires
% if 'interactive' is 'true', otherwise it returns hardcoded points

    if interactive
        % selection of keypoints with user interaction
        
        % we have to select a ROI (Region Of Interest) because
        % otherwise we don't find the features we desire (they are
        % unique only in a small area of the image)
        if strcmp(method, 'harris')
            % 'MinQuality': Minimum accepted quality of corners in [0,1],
            % 'FilterSize': Gaussian filter dimension in the range [3, min(size(I))]
            corners = detectHarrisFeatures(im_gray, 'MinQuality', 0.01, 'FilterSize', 9, 'ROI', [730 1250 1000 410]); 
            figure, imshow(im_gray); hold on;
            plot(corners.selectStrongest(300)); 
        elseif strcmp(method, 'surf')
            corners = detectSURFFeatures(im_gray, 'MetricThreshold', 50); 
            figure, imshow(im_gray); hold on;
            plot(corners.selectStrongest(1000)); 
        elseif strcmp(method, 'fast')
            corners = detectFASTFeatures(im_gray, 'MinQuality', 0.01, 'MinContrast', 0.2); 
            figure, imshow(im_gray); hold on;
            plot(corners.selectStrongest(1000)); 
        elseif strcmp(method, 'kaze')
            % 'region' (default) | 'sharpedge' | 'edge'
            corners = detectKAZEFeatures(im_gray, 'Threshold', 0.001, 'Diffusion', 'edge');
            figure, imshow(im_gray); hold on;
            plot(corners.selectStrongest(1000));  
        end
        
        % Use normal button clicks to add points. A shift-, right-, or double-click adds a final point and ends the selection. 
        % Pressing Return or Enter ends the selection without adding a final point. 
        % Pressing Backspace or Delete removes the previously selected point.

        [xi, yi] = getpts; % get the points selected by the user
        
        % here I filter out the points outside the image area
        indecesX = xi(:,1) >= 0;
        xi = xi(indecesX,1);
        yi = yi(indecesX,1);
        indecesY = yi(:,1) >= 0;
        xi = xi(indecesY,1);
        yi = yi(indecesY,1);
        if size(xi,1)~=size(yi,1) || size(xi,1)~=4 || size(yi,1)~=4
            error('You have to select exactly 4 points in the following order: up left, up right, down right, down left');
        end

        close; % close the figure
    else
        if image_sel == 1
            % HARDCODED POINTS

            % these points are the four vertices of the license plate taken into
            % account starting from the upper left one going clockwise        
            % LICENSE PLATE
            [xi] = [915.6; 1224; 1206; 901.9];
            [yi] = [1362; 1459; 1568; 1464];

            % rear parking sensors and rear bottom red lights
            % PARKING SENSORS AND LIGHTS
            %[xi] = [839.2; 1226; 1411; 788.6];
            %[yi] = [1287; 1414; 1620; 1398];
            
            % rear upper plate points and rear lower red lights
            % LICENSE PLATE AND LIGHTS
            %[xi] = [915.6; 1224; 1411; 788.6];
            %[yi] = [1362; 1459; 1620; 1398];
        else
             error("NO hardcoded keypoints for this Image");
        end
    end
end