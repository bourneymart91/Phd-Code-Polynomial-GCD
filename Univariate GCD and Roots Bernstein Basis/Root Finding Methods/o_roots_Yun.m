function [root_mult_array] = o_roots_Yun(fx)
%
% % Inputs
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% % Outputs
%
% root_mult_array : (Matrix)

fx_dash = Bernstein_Differentiate(fx);

m =  GetDegree(fx);
n = m - 1;

% Set upper and lower limits of GCD(f,g) # Since number of distinct roots
% is unknown, upper and lower limits are unknown.
lowerLimit_t = 1;
upperLimit_t = min(m,n);
limits_t = [lowerLimit_t, upperLimit_t];


% Perform GCD computation.
[fx_o,gx_o,dx, ux_o, vx_o, alpha, theta, t] ...
    = o_gcd_2Polys_mymethod(fx, fx_dash, limits_t);



i = 0;

arr_wx{i+1} = dx;
arr_ux{i+1} = Deconvolve(fx, arr_wx{i+1});
arr_vx{i+1} = Deconvolve(fx_dash, arr_wx{i+1});
arr_tx{i+1} = arr_vx{i+1} - Bernstein_Differentiate(arr_ux{i+1});

i = 1;



while (GetDegree(arr_tx{i}) >= 1)
    
    m = GetDegree(arr_ux{i});
    n = GetDegree(arr_tx{i});
    limits_t = [1 min(m,n)];
    
    [fx_o,gx_o,dx, ux_o, vx_o, alpha, theta, t] ...
    = o_gcd_2Polys_mymethod(arr_ux{i}, arr_tx{i}, limits_t);
    
    
    arr_wx{i+1,1} = dx;
    
    arr_ux{i+1,1} = Deconvolve(arr_ux{i}, arr_wx{i+1});
    arr_vx{i+1,1} = Deconvolve(arr_vx{i}, arr_wx{i+1});
    
    arr_tx{i+1,1} = arr_vx{i+1} - Bernstein_Differentiate(arr_ux{i+1});
    
    i = i+1;
    
end












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




end



