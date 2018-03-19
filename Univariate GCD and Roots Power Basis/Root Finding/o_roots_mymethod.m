function [root_multiplicity_matrix] = o_roots_mymethod(fx)
% Given the polynomial f(x) calculate its real roots by square free
% decomposition.
%
% Inputs.
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% Outputs.
%
% root_mult_array : (Matrix) Output two columns, the first containing the root, 
%   the second its corresponding multiplicity.

global SETTINGS


% Initialise an iteration counter
ite = 1;

% Initialise an array 'q' which stores the gcd outputs from each gcd
% calculation
arr_fx{1,1} = fx;

% let vGCD_Degree store the degrees corresponding to the array of
% GCDs stored in q.
vDegt_fx(1,1) = GetDegree(arr_fx{1});

% Let theta_vec store all theta values used in each iteration.
vTheta(1,1) = 1;

% Get the number of distinct roots of f_{1}. Since this is unknown at this
% time, set number of distinct roots to be m_{1} = deg(f_{1}).
vNumberDistinctRoots_fx(1,1) = GetDegree(arr_fx);

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.

rank_range = [-16 0];


while GetDegree(arr_fx{ite,1}) > 0
    
    % If degree of f(ite_num) is greater than one
    if vDegt_fx(ite,1) > 1
        
        
        fprintf(['\n' mfilename ' : ' sprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite-1,ite-1)]);
        
        
        % Perform GCD computation
        
        % Get upper and lower bounds of the next GCD computation
        % M_{i+1} > M_{i} - d_{i-1}
        if ite > 1
        
            lowerLimit_t = vDegt_fx(ite) - vNumberDistinctRoots_fx(ite-1);
            upperLimit_t = vDegt_fx(ite) - 1;


            
        else
            
            lowerLimit_t = 1;
            upperLimit_t = vDegt_fx(ite)-1;
            
        end
        

        fprintf([ mfilename ' : ' sprintf('Minimum degree of f_{%i}: %i \n', ite, lowerLimit_t)]);
        fprintf([ mfilename ' : ' sprintf('Maximum degree of f_{%i}: %i \n\n', ite, upperLimit_t)]);
        

        
        
        %        
        [arr_fx{ite,1},~,arr_fx{ite+1,1}, arr_ux{ite,1} ,arr_vx{ite,1}, ~, vTheta(ite+1,1),vDegt_fx(ite+1,1),~,~, rank_range] ...
            = o_gcd_mymethod_Univariate_2Polys( arr_fx{ite,1} , Differentiate(arr_fx{ite,1}) , [lowerLimit_t,upperLimit_t], rank_range);
        
        
        fprintf([ mfilename ' : ' sprintf('Computed degree of f_{%i}: %i \n', ite, vDegt_fx(ite+1) )]);
        
        % Get number of distinct roots of f(ite), given by
        vNumberDistinctRoots_fx(ite,1) = vDegt_fx(ite,1) - vDegt_fx(ite+1,1);

        fprintf([ mfilename ' : ' sprintf('Number of distinct roots in f_{%i} : %i \n',ite,vNumberDistinctRoots_fx(ite))]);
        
        % increment iteration number.
        ite = ite+1;
        
        LineBreakMedium();
        
    elseif vDegt_fx(ite) == 1
        % if m = 1, then n = 0, GCD has maximum degree 0.
        dx = 1;
        
        %theta_vec(ite_num+1) = 1;
        vDegt_fx(ite+1) = 0;
        
        arr_fx{ite+1,1} = 1;
        arr_ux{ite,1} = arr_fx{ite};
        ite = ite+1;
        
        break;
        
        
    end
end


% Get the degree structure of the polynomials h_{i}
vDeg_arr_hx = diff(vDegt_fx);

% Get the degree structure of the polynomials w_{i}
vDeg_arr_wx = diff([vDeg_arr_hx; 0]);

vMultiplicities = find(vDeg_arr_wx~=0);

% %
% %
% %
% %
% Deconvolve the first set of polynomials q_{i}(x), which are the outputs
% of the series of GCD computations, to obtain the set of polynomaials h_{x}

% We can either obtain h(x) from the series of deconvolutions on f(x){i} or
% we can use the already calculated values ux{i}
%method = 'By Deconvolution';

switch SETTINGS.ROOTS_HX_COMPUTATION_METHOD
    
    case 'From Deconvolutions'
        
        fprintf([mfilename ' : ' sprintf('Deconvolution Method : %s',SETTINGS.DECONVOLUTION_METHOD_FX_HX)]);
        arr_hx = Deconvolve_Set(arr_fx,SETTINGS.DECONVOLUTION_METHOD_FX_HX);
        
        
    case 'From GCD Computation'
        fprintf([mfilename ' : ' sprintf('Deconvolution Method : %s',SETTINGS.ROOTS_HX_COMPUTATION_METHOD)]);
        arr_hx = arr_ux;
        
    otherwise
        str = sprintf([mfilename ' : ' sprintf('ROOTS_UX is either From Deconvolutions or From ux') '\n']);
        error(str);
        
end


% Get the number of polynomials in h_{x}
[nPolys_arr_hx] = size(arr_hx,1);

% If only one entry in h_{x}
if nPolys_arr_hx == 1
    
    
    % if number of cols in h1 is only 1, then do not perform second
    % deconvolution, since only one entry in h1.
    % Note - this is a rare exception.
    vRoots = [];
    
    factor_x = arr_hx{1};
    
    % Normalise the polynomial coefficients by the leading coefficient x^m
    factor_x = factor_x./factor_x(end);
    
    % Get the root from the factor ( x - r);
    rt = - factor_x(1);
    
    % get the roots with respect to y, and their multiplicities all set
    % to one.
    arr_wx = rt;
    
    % add the roots to the array of roots
    vRoots = [vRoots ; arr_wx];
    
else
    % perform deconvolutions
    
    % Deconvolve the second set of polynomials
    arr_wx = Deconvolve_Set(arr_hx, SETTINGS.DECONVOLUTION_METHOD_HX_WX);
    
    % w1 yields the simple, double, triple roots of input polynomial f.
    % w1{i} yields the roots of multiplicity i.
    
    % set the w1{max} = h1{max}
    arr_wx{ite-1,1} = arr_hx{ite-1};
    
    % get number of entries in w1
    [nPolys_arr_wx] = size(arr_wx,1);
    
    % initialise an empty set
    vRoots = [];
    
    % for each multiplicity in w1.
    for i = 1:1:nPolys_arr_wx
        
        % Get the polynomial w_{i}
        poly_wi = arr_wx{i};
        
        % Get the degree of w_{i}
        vDeg_arr_wx(i,1) = GetDegree(arr_wx{i});
        
        
        % If w_{i} is of degree one, then is of the form (ax+b)
        % and has only one root. Add it to the list of roots.
        if vDeg_arr_wx(i,1) == 1;
            
            
            % Normalise the coefficients of w_{i} by dividing by the
            % leading coefficient. Since coefficients are orders in
            % asending powers, LC is the second coefficient.
            poly_wi = poly_wi./poly_wi(2);
            
            % Get the root
            rt = - poly_wi(1);
            
            % Add the root to the [root, mult] matrix
            vRoots = [vRoots ; rt];
            
        elseif (vDeg_arr_wx > 1)
            % The given multiplicity contains more than one root, such that
            % number of coefficients in greater than 2, use MATLAB roots
            % function to find roots.
            
            % Get the roots in terms of w_{i} using Matlab Roots function.
            rts = roots(flipud(poly_wi));
            
            % Add the computed roots to the array of roots.
            vRoots = [vRoots ; rts];
        end
    end
end

% % Get multiplicities of the roots
% Obtaining multiplicities of the calculated roots



% create a matrix where the first column contains the multiplicities, and
% the second column contains the number of roots of that multiplicity
nPolys_wi = size(vDeg_arr_wx,1);
mat = [(1:1:nPolys_wi)' vDeg_arr_wx];

count = 1;
root_multiplicity_matrix = [];
for i = 1:1:size(mat,1)
    % Get the number of roots of multiplicity i
    nRoots_of_Multi = mat(i,2);
    for j = 1:1:nRoots_of_Multi
        root_multiplicity_matrix = [root_multiplicity_matrix ; vRoots(count,1) i];
        count= count +1;
    end
end

% % Print the calculated roots and the corresponding multiplicities.


end