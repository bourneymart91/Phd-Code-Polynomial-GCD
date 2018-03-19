function [fxy_o, gxy_o, dxy_o, uxy_o, vxy_o, t, rank_range] = ...
    o_gcd_mymethod_Bivariate_2Polys(fxy, gxy, m, n, limits_t, rank_range)
% Compute the gcd of two input polynomials f(x,y) and g(x,y), where f(x,y)
% and g(x,y) are bivariate polynomials in Bernstein form, with total degree
% m and n respectively.
%
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of bivariate polynomial f(x,y) in the Bernstein form.
%
% gxy : (Matrix) Coefficients of bivariate polynomial g(x,y) in the Bernstein form.
%  
% m : (Int) Total degree of f(x,y)
%
% n : (Int) Total degree of g(x,y)
%
% limits_t : (Int Int)
%
% rank_range : [Float Float]
%
% % Outputs
%
% fxy : (Matrix) Coefficients of the bivariate polynomial f(x,y) 
%
% gxy : (Matrix) Coefficients of bivariate polynomial g(x,y)
%
% dxy : (Matrix) Coefficients of the GCD of f(x,y) and g(x,y)
%
% uxy : (Matrix) Coefficients of the polynomial u(x,y) where f(x,y)/d(x,y) = u(x,y)
%
% vxy : (Matrix) Coefficients of the polynomial v(x,y) where g(x,y)/d(x,y) = v(x,y)
% 
% t : (Int) Total degree of the GCD d(x,y)


% Get the degree of the GCD using only the first sylvester subresultant
% matrix
%boolComputeDegreeSingularValuesS1 = true;
%if ( boolComputeDegreeSingularValuesS1 == true)
%    t = ComputeDegreeSingularValuesS1(fxy, gxy, m, n);
%end




% Compute the degree of the GCD
[t, GM_fx, GM_gx, alpha, th1, th2, rank_range] = ...
    GetGCDDegree_Bivariate_2Polys(fxy, gxy, m, n, limits_t, rank_range);


% % 
% Normalize fxy and gxy to obtain fxy_n and gxy_n
fxy_n = fxy ./ GM_fx;
gxy_n = gxy ./ GM_gx;

% Get low rank approximation of S_{t}(f,g) and use perturbed coefficients.
% Update f(w,w) and g(w,w).
[fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
    GetLowRankApproximation_2Polys(fxy_n, gxy_n, alpha, th1, th2, m, n, t);


% %
% Get the GCD polynomial d(x,y)
[fxy_lra, gxy_lra, uxy_lra, vxy_lra, dxy_lra, alpha_lra, th1_lra, th2_lra] = ...
    APF_2Polys(fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr, m, n, t);


% Get outputs
fxy_o = fxy_lra;
gxy_o = gxy_lra;
dxy_o = dxy_lra;
uxy_o = uxy_lra;
vxy_o = vxy_lra;


end



function [t] = ComputeDegreeSingularValuesS1(fxy, gxy, m, n)


    GM_fx = GetMean(fxy, m, n-1);
    GM_gx = GetMean(gxy, n, m-1);
    
    % Divide entries of f(x,y) and g(x,y) by geometric mean
    fxy_n = fxy ./ GM_fx;
    gxy_n = gxy ./ GM_gx;
    
    [alpha, th1, th2] = Preprocess(fxy_n, gxy_n, m, n, 1);
    
    % Get f(x,y) and g(x,y) with thetas to get f(w_{1},w_{2}) and g(w_{1},w_{2})
    fww = GetWithThetas(fxy_n, m, th1, th2);
    gww = GetWithThetas(gxy_n, n, th1, th2);
    
    
    S1 = BuildSylvesterMatrix_2Polys(fww, alpha.*gww, m, n, 1);
    
    vSingularValues = svd(S1);
    vMetric = log10(vSingularValues);
    t = GetGCDDegree_OneSubresultant(vMetric);
    
end