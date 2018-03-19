function [root_mult_matrix] = o_roots_Musser(fx)
% Given the polynomial f(x) compute its roots by Square free factorization.
% This algorithm is referred to as Musser's Algorithm in 
% http://www.cs.berkeley.edu/~fateman/282/F%20Wright%20notes/week6.pdf 
%
%
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x) as a column vector where first
%   coefficient has lowest power [a_{0} ; ... ; a_{m}]^{T}
%
% % Outputs 
% 
% root_mult_array : (Matrix) A matrix of roots and multiplicities, where the first
%   column contains all of the roots of f(x) and the second contains the
%   corresponding multiplicities.

global SETTINGS


% Set Iteration number
ite = 1;

% 
f{1} = fx;

% Get the derivative of f(x)
gx{1} = Differentiate(f{1});

% Get the degree of f(x)
m = GetDegree(f{1});

n = GetDegree(gx{1});

% Set upper and lower limits of GCD(f,g) # Since number of distinct roots
% is unknown, upper and lower limits are unknown.
lowerLimit_t = 1;
upperLimit_t = min(m,n);
limits_t = [lowerLimit_t, upperLimit_t];

rank_range = [0 -16];

% Perform GCD computation.
[fx_n,gx_n,dx, ux_o, vx_o, alpha, theta, t , lambda,mu, rank_range] ...
    = o_gcd_mymethod_Univariate_2Polys(f{1}, gx{ite}, limits_t, rank_range);


LineBreakMedium();

arr_gx{ite} = dx;

arr_hx{ite} = Deconvolve(f{ite},arr_gx{ite});

while (GetDegree(arr_hx{ite}) > 0 )
    
    % Get the degree of polynomial f(x)
    m = GetDegree(arr_gx{ite});
    
    % Get the degree of polynomial g(x)
    n = GetDegree(arr_hx{ite});
    
    % Set Limits
    lowerLimit_t = 1;
    upperLimit_t = min(m,n);
    limits_t = [lowerLimit_t, upperLimit_t];
    
    if (GetDegree(arr_hx{ite}) ==0 || GetDegree(arr_gx{ite}) == 0)
        arr_hx{ite+1} = 1;
        arr_vx{ite} = arr_hx{ite};
        arr_ux{ite} = arr_gx{ite};
    else
        
          
          [fx_n, gx_n, arr_hx{ite+1}, arr_ux{ite}, arr_vx{ite}, ~, ~, ~ , ~,~, rank_range] ...
              = o_gcd_mymethod_Univariate_2Polys(arr_gx{ite}, arr_hx{ite}, limits_t, rank_range);
        
    end
    
    % The polynomial g can be obtained in two ways, as u(x) from the GCD
    % triple (d(x),u(x),v(x)) or by deconvolution.
   
    
    switch SETTINGS.ROOTS_HX_COMPUTATION_METHOD_IN_MUSSER_ALGORITHM

        case 'From GCD Computation'
            arr_gx{ite+1} = arr_ux{ite};
            
        case 'From Deconvolution'
            arr_gx{ite+1} = Deconvolve(arr_gx{ite}, arr_hx{ite+1});
            
        otherwise
            error('err')
    end
    
    ite = ite+1;
    
    LineBreakMedium();
    
end

root_mult_matrix = [];

for i = 1:1:length(arr_vx)
    
    try
        %fprintf('Roots of multiplicity %i \n',i)
        
        factor = arr_vx{i};
        m = GetDegree(factor);
        if m > 1
            % Get roots by matlab method
            root = roots(flipud(factor));
            
            % Get number of roots
            nRoots = length(root);
            
            % Get the new roots and mutliplicities
            new_roots_mults = [root i.*ones(nRoots,1)];
            
            % Add the new roots to the array of roots
            root_mult_matrix = [root_mult_matrix ; new_roots_mults];
            
        else
            % Divide by x coefficient
            factor = factor./factor(2);
            % Get root
            root = -1.*factor(1);
            
            new_root_mult = [root i];
            
            % Add the new root to the array of roots
            root_mult_matrix = [root_mult_matrix ; new_root_mult];
            
            
        end

        
    catch
    end
end

end