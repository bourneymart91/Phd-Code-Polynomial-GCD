function [root_multiplicity_matrix, arr_fx, arr_hx, arr_wx] = ...
    o_roots_mymethod(fx)
% Calculate the roots of the input polynomial f(x) in Bernstein form.
%
% Inputs.
%
% fx:   Vector of polynomial coefficients in Bernstein basis, where the
%       ith entry of fx, fx(i) is the coefficient a_{i}B^{m}_{i}, and m is the
%       degree of fx.
%
% Outputs.
%
% root_mult_array : The Calculated roots of the polynomial f(x) and their
% corresponding multiplicities.

global SETTINGS

% Get the array of polynomials f_{i}(x) by series of GCD computations
arr_fx = Get_arr_fx(fx);

% Get number of polynomials f_{i}(x)
nPolys_fx = length(arr_fx);

arr_fx_n = cell(nPolys_fx,1);

% Normalised polynomials f_{i}(x)
for i = 1:1:length(arr_fx)
    
    % Get ith polynomial
    fx = arr_fx{i};
    
    % Normalise
    fx = fx ./ fx(1);
    
    % Save normalised polynomial to array
    arr_fx_n{i} = fx;
    
end


% Deconvolve the set of polynomials f_{i}(x) to get the set of polynomials
% h_{i}(x)
arr_hx = Deconvolve_Set(arr_fx_n, SETTINGS.DECONVOLUTION_METHOD_HX);


% Get array w_{i}(x), then the set of roots and corresponding multplicities
[arr_wx, root_multiplicity_matrix] = GetArray_wx(arr_hx);



end


function arr_fx = Get_arr_fx(fx)
%
% % Inputs
%
% fx : (Vector) Vector of coefficients of polynomial f_{0}(x)
%
% % Outputs
%
% arr_fx : (Array of Vectors) Vectors containing coefficients of the set of
% polynomials f_{i}(x)
 
% Initialise an iteration counter for looping through GCD computations.
ite = 1;


% Initialise an array 'q' which stores the gcd outputs from each gcd
% calculation
arr_fx{1,1} = fx;

% let degree_vector store the degrees corresponding to the array of
% GCDs stored in q.
vDegree_arr_fx(1) = GetDegree(arr_fx{1});

% Get the number of distinct roots of f_{1}. Since this is unknown at this
% time, set number of distinct roots to be m_{1} = deg(f_{1}).
vNumberDistinctRoots(1) = GetDegree(arr_fx{1});

% Initialise rank_range metric used in computing the degree of the \GCD.
% This is iteratively updated where the lower limits is the numeric value
% \rho_{t} and the upper limit is the numeric value \rho_{t+1}
rank_range = [0 0];


% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.

while vDegree_arr_fx(ite) > 0
    
    
    % If degree of f_{i} is greater than one
    if vDegree_arr_fx(ite) > 1
        
        
        fprintf([mfilename ' : ' sprintf('Compute f_{%i} the GCD of f_{%i} and derivative f_{%i}\n\n',ite, ite - 1, ite - 1)]);
        
        
        if ite > 1
            [lowerLimit_t, upperLimit_t] =...
                GetLimits(vDegree_arr_fx(ite), vNumberDistinctRoots(ite-1));
        else
            [lowerLimit_t, upperLimit_t] =...
                GetDefaultLimits(vDegree_arr_fx(ite));
        end
        
        
        
        
        % Get GCD
        
        % With replacement of f_{i}(x)
        [arr_fx{ite}, ~, arr_fx{ite+1,1}, ~,~,~,~, vDegree_arr_fx(ite+1), rank_range] = ...
            o_gcd_2Polys_mymethod( arr_fx{ite}, Bernstein_Differentiate(arr_fx{ite}), [lowerLimit_t, upperLimit_t], rank_range);
        
        % Without replacement of f_{i}(x)
        %[~, ~, arr_fx{ite+1,1}, ~,~,~,~, vDegree_arr_fx(ite+1), rank_range] = ...
        %    o_gcd_2Polys_mymethod( arr_fx{ite}, Bernstein_Differentiate(arr_fx{ite}), [lowerLimit_t, upperLimit_t], rank_range);
        
        
        
        % Get number of distinct roots of f(ite)
        vNumberDistinctRoots(ite) = ...
            vDegree_arr_fx(ite) - vDegree_arr_fx(ite+1);
        
        fprintf([mfilename ' : ' sprintf('Degree of f_{%i} : %i \n',ite, vDegree_arr_fx(ite+1))]);
        
        % increment iteration number.
        ite = ite + 1;
        
        
    elseif vDegree_arr_fx(ite) == 1
        
        % if m=1, then n = 0, GCD has maximum degree 0.
        fprintf([mfilename ' : ' 'Only one subresultant exists \n'])
        dx = 1;
        
        
        vDegree_arr_fx(ite+1) = length(dx)-1;
        arr_fx{ite+1} = dx;
        arr_ux{ite} = arr_fx{ite};
        ite = ite+1;
        
        break;
        
        
    end
end

end


function [lowerLimit_t, upperLimit_t] = GetLimits(m, d)
%
% % Inputs
%
% m : (Int)
%
% d : (Int)

% Get upper and lower bounds of the GCD Computation.
% M_{i+1} > M_{i} - d_{i-1}
try
    
    
    
    lowerLimit_t = m - d;
    if lowerLimit_t < 0
        lowerLimit_t = 0;
    end
    upperLimit_t = m - 1;
    
    fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i}: %i \n', ite, lowerLimit_t)]);
    fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i}: %i \n\n', ite, upperLimit_t)]);
    
catch
    GetDefaultLimits(m)
    
end

end

function [lowerLimit_t, upperLimit_t] = GetDefaultLimits(m)

lowerLimit_t = 1;
upperLimit_t = m-1;

end


function [arr_wx, root_mult_matrix] = GetArray_wx(arr_hx)
%
% % Inputs
%
% arr_hx : (Array of Vectors)
%
% % Outputs
%
% arr_wx : (Array of Vectors)
%
% root_mult_matrix

% Get number of polynomials in the array h_{i}(x)
[nPolynomialsArray_hx] = size(arr_hx, 1);

global SETTINGS


% Initialise an empty root array
root_mult_matrix = [];

if nPolynomialsArray_hx == 1
    
    % if number of cols in h1 is only 1, then do not perform second set of
    % deconvolutions, since only one entry in h1.
    
    % Get the polynomial a(w)
    wx = arr_hx{1};
    
    % Get degree of w(x)
    n = GetDegree(wx);
    
    % Normalise the coefficients of a(w)
    wx = wx ./ wx(1);
    
    % Get roots of factor w(x)
    vRoots = GetRoots(wx);
    vMultiplicity = 1*ones(n,1);
    
    % Add the roots to the array of roots
    root_mult_matrix = [root_mult_matrix ; [vRoots vMultiplicity]];
    
    
else
    
    
    % Deconvolve the second set of polynomials h_{i}(x) to obtain the set of
    % polynomials w_{i}(x)
    arr_wx = Deconvolve_Set(arr_hx, SETTINGS.DECONVOLUTION_METHOD_WX);
    
    % set the w1{max} = h1{max}
    arr_wx{nPolynomialsArray_hx, 1} = arr_hx{nPolynomialsArray_hx, 1};
    
    
    % Get number of entries in w1
    nPolynomialsArray_wx = nPolynomialsArray_hx;
    
    % For each polynomial w_{i}(x)
    for i = 1 : 1 : nPolynomialsArray_wx
        
        
        % Get the polynomial a(w), whose simple roots have multiplicty
        % i in the polynomial f which we started with.
        wx = arr_wx{i};
        
        % Normalise the polynomial coefficients
        wx = wx./wx(1);
        
        n = GetDegree(wx);
        
        % if the polynomial of said multiplicity is of length 2, degree 1, then
        % only one root exists for this multiplicity. Add it to the list wp1.
        if (n == 1)
            
            % Convert to power form, so that coefficients are in terms of y^{i}
            % rather than (1-y)^{m-i}y^{i}.
            wx_pwr = [wx(1,:) ; wx(2,:) - wx(1,:)];
            
            % Obtain the root in terms of y, and set multiplicity to one.
            root_wx = -wx_pwr(1,:) ./ wx_pwr(2,:);
            root_multiplicity = i;
            
            root_mult_pair = [root_wx root_multiplicity];
            
            % Add the root to the [root, mult] matrix
            root_mult_matrix = [root_mult_matrix ; root_mult_pair];
            
        elseif (n > 1)
            % The given multiplicity contains more than one root, such that
            % number of coefficients in greater than 2, use MATLAB roots
            % function to find roots.
            % display('Multiplicity contains more than one root. ')
            
            vRoots = GetRoots(wx);
            vMultiplicity = i*ones(n, 1);
            
            % add the roots to the array of roots
            root_mult_matrix = [root_mult_matrix ; [vRoots vMultiplicity]];
            
        end
    end
end

end


function roots_wrt_x = GetRoots(wx)



% get aw including its binomial coefficients.
ax_bi = GetWithBinomials(wx);

% get the roots in terms of z^{i} where z^{i} =(\frac{y}{(1-y)})^{i}.
% Note we must flip the coefficients to conform with input format of
% MATLAB roots function.
rt_wrt_z = roots(flipud(ax_bi));

% Get the roots in terms of y.
roots_wrt_x = rt_wrt_z ./ (1 + rt_wrt_z); %Edit 27/07

% Get the roots with respect to y, and their multiplicities all set
% to one.
roots_wrt_x = roots_wrt_x ;



end