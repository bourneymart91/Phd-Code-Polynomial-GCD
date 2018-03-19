function [root_mult_array] = o_roots_Musser(fx)
% Given the polynomial f(x) compute its roots by Square free factorization.
% This algorithm is referred to as Musser's Algorithm in
% http://www.cs.berkeley.edu/~fateman/282/F%20Wright%20notes/week6.pdf
%
%
% % Inputs
%
% fx : Coefficients of polynomial f(x) as a column vector where first
% coefficient has lowest power [a_{0} ; ... ; a_{m}]^{T}
%
% % Outputs
%
% root_mult_array : A matrix of roots and multiplicities, where the first
% column contains all of the roots of f(x) and the second contains the
% corresponding multiplicities.

global SETTINGS

% Set Iteration number
ite = 1;

%
arr_fx{1} = fx;

% Get the derivative of f(x)
arr_gx{1} = Bernstein_Differentiate(arr_fx{1});

% Get the degree of f(x)
m = GetDegree(arr_fx{1});

% Get the degree of g(x)
n = GetDegree(arr_gx{1});

% Set upper and lower limits of GCD(f,g) # Since number of distinct roots
% is unknown, upper and lower limits are unknown.
lowerLimit_t = 1;
upperLimit_t = min(m,n);
limits_t = [lowerLimit_t, upperLimit_t];

rank_range = [-10,0];

% Perform GCD computation.
[fx_o,gx_o,dx_o, ux_o, vx_o, alpha, theta, t] ...
    = o_gcd_2Polys_mymethod(arr_fx{1}, arr_gx{ite}, limits_t, rank_range);

dx = dx_o;
ux = ux_o;

LineBreakMedium();
arr_gx{ite,1} = dx;

arr_hx{ite,1} = Deconvolve(arr_fx{1},arr_gx{ite});

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
        arr_hx{ite+1,1} = 1;
    else
        
        %[fx_n,gx_n,dx, ux, ~, ~, ~, ~ ] ...
        %    = o_gcd_mymethod(g{ite},h{ite},deg_limits);
        
        [fx_o, gx_o, dx_o, ux_o, vx_o, alpha_o, theta_o, t, rank_range] ...
            = o_gcd_2Polys_mymethod(arr_gx{ite}, arr_hx{ite}, limits_t, rank_range);
        
        dx = dx_o;
        ux = ux_o;
        
        arr_hx{ite+1,1} = dx;
    end
    
    % The polynomial g can be obtained in two ways, as u(x) from the GCD
    % triple (d(x),u(x),v(x)) or by deconvolution.
    
    
    switch SETTINGS.GET_HX_METHOD
        case 'From ux'
            arr_gx{ite+1,1} = ux;
            
        case 'From Deconvolution'
            arr_gx{ite+1,1} = Deconvolve(arr_gx{ite},arr_hx{ite+1});
            
        otherwise
            error('err')
    end
    
    %fprintf([mfilename ' : ' sprintf('g_{%i} degree : %i \n',ite+1,GetDegree(g{ite+1})) ]);
    %fprintf([mfilename ' : ' sprintf('h_{%i} degree : %i \n',ite+1,GetDegree(h{ite+1})) ]);
    
    %arr_mx{ite,1} = Deconvolve(arr_hx{ite},arr_hx{ite+1});
    ite = ite+1;
    
    LineBreakMedium();
    
end

w_batch = Deconvolve_Set(arr_hx, 'Batch');
arr_wx = w_batch;


%
root_mult_array = [];

for i = 1:1:length(arr_wx)
    
    try
        %fprintf('Roots of multiplicity %i \n',i)
        
        factor = arr_wx{i};
        m = GetDegree(factor);
        if m == 0
            % do nothing
        elseif m > 1
            % Get roots by matlab method
            root = roots(flipud(factor));
            
            % Get number of roots
            nRoots = length(root);
            
            % Get the new roots and mutliplicities
            new_roots_mults = [root i.*ones(nRoots,1)];
            
            % Add the new roots to the array of roots
            root_mult_array = [root_mult_array ; new_roots_mults];
            
        else
            % Divide by x coefficient
            factor = factor./factor(2);
            
            % Convert to power form, so that coefficients are in terms of y^{i}
            % rather than (1-y)^{m-i}y^{i}.
            a_pwr = [factor(1,:) ; factor(2,:)-factor(1,:)];
            
            % Obtain the root in terms of y, and set multiplicity to one.
            root = -a_pwr(1,:)./a_pwr(2,:);
            
            new_root_mult = [root i];
            
            % Add the new root to the array of roots
            root_mult_array = [root_mult_array ; new_root_mult];
            
            
        end
        
        
    catch
    end
end
