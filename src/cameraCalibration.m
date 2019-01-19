function K = cameraCalibration()

    
    [A,b] = equationsToMatrix(eqns);

    %{
    %A = double(A);
    %b = double(b);
    %sols = (A' * A) \ (A' * b);

    sols = double(A\b);

    w = [sols(1)    0           sols(2);
         0          1           sols(3);
         sols(2)    sols(3)     sols(4)];

    w = norm_matrix.'*w*norm_matrix;
    %}


    % fit a linear model without intercept
    lm = fitlm(double(A),double(b), 'y ~ x1 + x2 + x3 + x4 - 1');
    % get the coefficients
    W = lm.Coefficients.Estimate;

    % image of absolute conic
    w = double([W(1,1) 0 W(2,1); 0 1 W(3,1); W(2,1) W(3,1) W(4,1)]);
    w = norm_matrix*w*norm_matrix;

    a = sqrt(w(1, 1));
    u = -w(1, 3)/w(1, 1);
    v = -w(2, 3);
    fy = sqrt(w(3, 3) - w(1, 1)*u^2 - v^2);
    fx = fy/a;

    K = [fx,    0,      u;
         0,     fy,     v;
         0,     0,      1];
 
end