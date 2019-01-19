function linePlot(lines, color, linewidth)
% Colors: 
% 'red' or 'r'          Red
% 'green' or 'g'        Green	
% 'blue' or 'b'         Blue	
% 'yellow' or 'y'       Yellow	
% 'magenta' or 'm'      Magenta	
% 'cyan' or 'c'         Cyan	
% 'white' or 'w'        White	
% 'black' or 'k'        Black	
% 'none'                No color

    switch nargin
        case 1
           color = 'cyan';
           linewidth = 1;
        case 2
           linewidth = 1;          
    end

    for i = 1:size(lines, 1)
        myline = lines(i,:); % single line in the form of a vector 1x4 [x1 y1 x2 y2]
        
        line([myline(1) myline(3)], [myline(2) myline(4)], 'LineWidth', linewidth, 'Color', color);
    end