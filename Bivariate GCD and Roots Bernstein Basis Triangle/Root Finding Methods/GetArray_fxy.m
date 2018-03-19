function [arr_fxy, vDeg_t_arr_fxy] =  GetArray_fxy(fxy_matrix,M)
%
% % Inputs
%
% fxy_matrix : (Matrix) Coefficients of f(x,y)
%
% M : (Int) Degree of f(x,y)
%
% % Outputs
%
% arr_fxy : (Array of Matrices)



% Set the iteration number
ite = 1;

% Set the first entry of q to be the input polynomial f(x,y)
arr_fxy{1} = fxy_matrix;

% Get the total degree of f(x,y)
vDeg_t_arr_fxy(ite,1) = M(ite);

% Get the degree of f(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(arr_fxy{ite});
vDeg_x_arr_fxy(ite,1) = m1;
vDeg_y_arr_fxy(ite,1) = m2;

rank_range = [0 0];

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while vDeg_x_arr_fxy(ite,1) > 0
    
    
    if (vDeg_x_arr_fxy(ite,1) == 1)
        % The derivative with respect to x is a constant
        
        % The GCD is a constant
        %arr_fxy{ite+1,1} = Differentiate_wrt_x(arr_fxy{ite},vDeg_t_arr_fxy(ite));
        arr_fxy{ite+1,1} = Differentiate_wrt_x(arr_fxy{ite});
        
        % Deconvolve
        %         arr_uxy{ite+1,1} = Deconvolve_Bivariate(arr_fxy{ite}, arr_fxy{ite+1},...
        %             vDeg_t_arr_fxy(ite), vDeg_t_arr_fxy(ite)-1);
        arr_uxy{ite+1,1} = Deconvolve_Bivariate(arr_fxy{ite}, arr_fxy{ite+1});
        
        
        % Get total degree of d(x,y) and degree with respect to x and y
        vDeg_t_arr_fxy(ite+1,1) = vDeg_t_arr_fxy(ite) - 1;
        vDeg_x_arr_fxy(ite+1,1) = vDeg_x_arr_fxy(ite) - 1;
        vDeg_y_arr_fxy(ite+1,1) = vDeg_y_arr_fxy(ite);
        
        break;
    end
    
    
    LineBreakLarge()
    fprintf([mfilename ' : ' sprintf('GCD Calculation Loop iteration = %i \n', ite)]);
    fprintf([mfilename ' : ' sprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite-1,ite-1)]);
    
    fxy = arr_fxy{ite};
    
    % Get the derivative of f(x,y) with respect to x.
    %gxy = Differentiate_wrt_x(arr_fxy{ite},vDeg_t_arr_fxy(ite));
    gxy = Differentiate_wrt_x(arr_fxy{ite});
    hxy = Differentiate_wrt_y(arr_fxy{ite});
    
    % Get the total degree of f(x,y), g(x,y) and h(x,y)
    m =  vDeg_t_arr_fxy(ite);
    n = m - 1;
    o = m - 1;
    
    % Get the upper and lower limit of the degree of the GCD(f, f')
    if ite > 1
        lowerLimit_t = vDeg_t_arr_fxy(ite) - vNumberDistinctRoots(ite-1);
        upperLimit_t = m - 1;
    else
        lowerLimit_t = 1;
        upperLimit_t = m - 1;
    end
    
    LineBreakLarge();
    fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i}: %i \n', ite, lowerLimit_t)]);
    fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i}: %i \n', ite, upperLimit_t)]);
    LineBreakLarge();
    
    % GCD is only a scalar with respect to x so set equal to g(x,y).
    
    limits_t = [lowerLimit_t, upperLimit_t];
    
    
    
    [fxy_o, gxy_o, hxy_o, dxy_o, uxy_o, vxy_o, wxy_o, t, rank_range] = ...
        o_gcd_mymethod_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, limits_t, rank_range);
    
    %[arr_fxy{ite,1},~,arr_fxy{ite+1,1},arr_uxy{ite,1},arr_vxy{ite,1},t,t1,t2] = o_gcd_mymethod(arr_fxy{ite},gxy,m,n,);
    % arr_fxy{ite,1} = fxy_o;
    arr_fxy{ite+1,1} = dxy_o;
    arr_uxy{ite,1} = uxy_o;
    arr_vxy{ite,1} = vxy_o;
    
    
    % Get total degree of d(x,y) and degree with respect to x and y
    vDeg_x_arr_fxy(ite+1,1) = t;
    vDeg_y_arr_fxy(ite+1,1) = t;
    vDeg_t_arr_fxy(ite+1,1) = t;
    
    % Get number of distinct roots of f(ite)
    vNumberDistinctRoots(ite,1) = vDeg_t_arr_fxy(ite) - vDeg_t_arr_fxy(ite+1);
    
    LineBreakLarge()
    fprintf([mfilename ' : ' sprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,vDeg_t_arr_fxy(ite+1))]);
    fprintf([mfilename ' : ' sprintf('Number of distinct roots in f_{%i} : %i \n', ite, vNumberDistinctRoots(ite))]);
    fprintf([mfilename ' : ' sprintf('Degree of f_{%i} : %i \n',ite , vDeg_t_arr_fxy(ite+1))])
    LineBreakLarge()
    
    % Increment the iteration number
    ite = ite + 1;
    
end



end
