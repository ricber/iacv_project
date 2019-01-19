function vp = fitVp(lines)
% This function fits a vanishing point using multiple line by using LSA.
% 'lines' is a matrix 3*n containing all lines having a common direction
% Given the lines the vanishing points obeys the law: l.'*vp = 0
% considering the last element of vp equal to 1 => ax + by = -c
% where the line is [a b c]'

X = []; % matrix of nx2 elements (n is 'lines' size 2) 
Y = []; % matrix of nx1 elements

for i = 1:size(lines,2)
    % first get the line
    l = lines(:,i);
    
    % put the first two components inside x
    X(i, :) = [l(1), l(2)];
    Y(i, :) = -l(3);
    
end
W = (X.'*X)\(X.'*Y);
vp = [W(1,1), W(2,1), 1].';

