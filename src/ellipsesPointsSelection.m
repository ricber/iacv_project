function points = ellipsesPointsSelection(interactive, normalization, C1, H1, C2, H2)
    % this function takes into account two ellipses matrices and the
    % homographies to apply to them (rotation + translation) and returns 
    % the four tangency points of the two bitangent lines to the two ellipses
    % If 'interactive' is 'true' is computes the points real-time otherwise
    % it returns hardcoded points. 
    if interactive
        syms l1 l2
        l = [l1; l2; 1];
        eqns = [l.'*inv(inv(H1).'*C1*inv(H1))*l == 0, l.'*inv(inv(H2).'*C2*inv(H2))*l == 0];
        vars = [l1 l2];
        [sol_l1, sol_l2] = vpasolve(eqns, vars);

        L = [double(sol_l1.'); double(sol_l2.'); ones(1,4)];

        X = (inv(H1).'*C1*inv(H1))\L;
        Y = (inv(H2).'*C2*inv(H2))\L;
        
        pt_tng_ul = X(:,1)'; % upper (u) point (pt) tangent (tng) to the left (l) ellipse
        pt_tng_ur = Y(:,1)';
        pt_tng_dr = Y(:,2)';
        pt_tng_dl = X(:,2)';
        
        pt_cross_ul = X(:,4)'; % upper (u) point (pt) tangent (tng) to the left (l) ellipse
        pt_cross_dl = X(:,3)';
        pt_cross_dr = Y(:,4)';
        pt_cross_ur = Y(:,3)';
        
        points = [  pt_tng_ul;
                    pt_tng_ur;
                    pt_tng_dr;
                    pt_tng_dl;
                    pt_cross_ul;
                    pt_cross_ur;
                    pt_cross_dr;
                    pt_cross_dl];
    else
        % HARDCODED POINTS
        if normalization
            npt_tng_ul = [0.042910984364247;0.036621782724790;0.076778900152799]'; % upper (u) normalized (n) point (pt) tangent (tng) to the left (l) ellipse
            npt_tng_ur = [0.037183328035400;0.014780427238900;0.047241964279286]';
            npt_tng_dr = [-0.027520227328717;-0.014224691982798;-0.035365244782472]';
            npt_tng_dl = [-0.032857796645869;-0.036043613147457;-0.058878163395084]';

            npt_cross_ul = [0.058889183822435;0.049549978082241;0.104451342153868]';
            npt_cross_ur = [0.022895034435491;0.009190638771406;0.029168733294491]';
            npt_cross_dr = [-0.051120590754893;-0.026578035148216;-0.065862863747153]';
            npt_cross_dl = [-0.028693972041899;-0.030949242282073;-0.050993426930933]';

            points = [  npt_tng_ul;
                        npt_tng_ur;
                        npt_tng_dr;
                        npt_tng_dl;
                        npt_cross_ul;
                        npt_cross_ur;
                        npt_cross_dr;
                        npt_cross_dl];   
        else
            pt_tng_ul = [1.771365434563379e+02;1.133810393164662e+02;0.076778900153114]'; % upper (u) normalized (n) point (pt) tangent (tng) to the left (l) ellipse
            pt_tng_ur = [1.534927781303090e+02;45.760202731689590;0.047241964279341]';
            pt_tng_dr = [-1.136034984128143e+02;-44.039646378691586;-0.035365244782432]';
            pt_tng_dl = [-1.356369845547452e+02;-1.115910263049433e+02;-0.058878163395341]';

            pt_cross_ul = [2.430945508196192e+02;1.534067321430403e+02;0.104451342154129]';
            pt_cross_ur = [94.510702149577720;28.454217636234740;0.029168733294451]';
            pt_cross_dr = [-2.110257986364128e+02;-82.285596818959960;-0.065862863747221]';
            pt_cross_dl = [-1.184487165896862e+02;-95.818854105815000;-0.050993426931248]';

            points = [  pt_tng_ul;
                        pt_tng_ur;
                        pt_tng_dr;
                        pt_tng_dl;
                        pt_cross_ul;
                        pt_cross_ur;
                        pt_cross_dr;
                        pt_cross_dl]; 
        end  
    end
end