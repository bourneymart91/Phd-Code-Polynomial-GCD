function [Sk] = BuildSubresultant_3Polys(fx, gx, hx, k)
% BuildSubresultant_3Polys(fx, gx, hx, k)
%
% This function builds the k-th Sylvester subresultant matrix S_{k}(f,g),
% in the Bernstein Basis.
%
%
% % Inputs
%
% fx : (Vector) The coefficients of the polynomial f(x) 
%
% gx : (Vector) The coefficients of the polynomial g(x) 
%
% hx : (Vector) The coefficients of the polynomial h(x) 
%
% k:  (Int) Index of subresultant S_{k} to be built
%
% % Outputs
%
% Sk : (Matrix) The k-th Sylvester subresultant matrix S_{k}(f,g,h)

if (nargin ~= 4)
    error('Not enough input arguments')
end
% Global Variables.
global SETTINGS





switch SETTINGS.SYLVESTER_EQUATIONS
    
    case '2'
        Sk = BuildSubresultant_3Polys_2Eqns(fx, gx, hx, k);
        
    case '3'
        Sk = BuildSubresultant_3Polys_3Eqns(fx, gx, hx, k);
        
    otherwise
        error('err')
end


end



function [Sk] = BuildSubresultant_3Polys_2Eqns(fx, gx, hx, k)


% Get the degree of the polynomial f(x)
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);

global SETTINGS

switch SETTINGS.SYLVESTER_MATRIX_VARIANT
    
    case 'T'
        
        T = BuildT_3Polys_2Eqns(fx, gx, hx, k);
        Sk = T;
        
    case 'DT'
        
        D = BuildD_3Polys_2Eqns(m, n - k, o - k);
        T = BuildT_3Polys_2Eqns(fx, gx, hx, k);
        Sk = D*T;
        
    case 'TQ'
        
        
        T = BuildT_3Polys_2Eqns(fx, gx, hx, k);
        Q = BuildQ_3Polys(m, n, o, k);
        
        

        Sk = T*Q;
        
    case 'DTQ'
        
        D = BuildD_3Polys_2Eqns(m, n - k, o - k);
        T = BuildT_3Polys_2Eqns(fx, gx, hx, k);
        Q = BuildQ_3Polys(m, n, o, k);
        
        
        Sk = D*T*Q;
        
        
        
        
        
    otherwise
        error('SETTINGS.SYLVESTER_MATRIX_VARIANT is either standard or rearranged')
end


% % Split the matrix into R1 and R2 row-partitions
% R1 = Sk(1 : m + n - k + 1, :);
% R2 = Sk(m + n - k + 2 : end, :);
% 
% c1 = R1(:, 1:n - k + 1);
% c2 = R1(:, n - k + 2 : n + o - (2*k) + 2);
% c3 = R1(:, m + n - (2 *k) + 3 : end);
% 
% c4 = R2(:, 1 : n - k + 1);
% c5 = R2(:, n - k + 2 : n + o - (2*k) + 2);
% c6 = R2(:, m + n - (2*k) + 3 : end);
% 
% v_c1 = c1(c1~=0);
% v_c3 = c3(c3~=0);
% v_c5 = c5(c5~=0);
% v_c6 = c6(c6~=0);
% 
% ratio_c1 = GetMaxMinRatio(v_c1);
% ratio_c3 = GetMaxMinRatio(v_c3);
% ratio_c5 = GetMaxMinRatio(v_c5);
% ratio_c6 = GetMaxMinRatio(v_c6);
% 
% % Get all nonzero Entries in R1
% vR1 = R1(R1~=0);
% vR2 = R2(R2~=0);
% 
% maxR1 = max(max(abs(vR1)));
% minR1 = min(min(abs(vR1)));
% 
% maxR2 = max(max(abs(vR2)));
% minR2 = min(min(abs(vR2)));
% 
% ratio_R1 = maxR1 ./ minR1; 
% ratio_R2 = maxR2 ./ minR2;
% 
% LineBreakMedium()
% fprintf('Ratio R1 : %e \n', ratio_R1);
% fprintf('Ratio R2 : %e \n', ratio_R2);
% fprintf('Ratio R1/R2 : %e \n', ratio_R1 ./ ratio_R2);
% LineBreakMedium()
% fprintf('Ratio c1 : %e \n', ratio_c1);
% fprintf('Ratio c3 : %e \n', ratio_c3);
% fprintf('Ratio c5 : %e \n', ratio_c5);
% fprintf('Ratio c6 : %e \n', ratio_c6);
% LineBreakMedium()
% 
%     function ratio =  GetMaxMinRatio(v_c1)
% 
%         ratio = max(abs(v_c1)) ./ min(abs(v_c1));
%     
%     end


end




function [Sk] = BuildSubresultant_3Polys_3Eqns(fx, gx, hx, k)
%
% % Inputs
% 
% fx : (Vector) The coefficients of the polynomial f(x)
%
% gx : (Vector) The coefficients of the polynomial g(x)
%
% hx : (Vector) The coefficients of the polynomial h(x)
%
% k : (Integer) Index of subresultant matrix S_{k}(f,g,h) which is
% constructed by this function.
%
% % Outputs
% 
% Sk : (Matrix) The k-th subresultant matrix S_{k}(f,g,h)



% Get the degree of the polynomials f(x), g(x) and h(x)
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);

global SETTINGS

switch SETTINGS.SYLVESTER_MATRIX_VARIANT
    
    case 'T'
        
        
        T = BuildT_3Polys_3Eqns(fx, gx, hx, k);
        
        Sk = T;
        
    case 'DT'
        
        
        D = BuildD_3Polys_3Eqns(m, n, o, k);
        T = BuildT_3Polys_3Eqns(fx, gx, hx, k);
        
        
        Sk = D*T;
        
    case 'TQ'
        
        T = BuildT_3Polys_3Eqns(fx, gx, hx, k);
        Q = BuildQ_3Polys(m, n, o, k);
        
        
        Sk = T*Q;
        
    case 'DTQ'
         
        
        D = BuildD_3Polys_3Eqns(m, n, o, k);
        T = BuildT_3Polys_3Eqns(fx, gx, hx, k);
        Q = BuildQ_3Polys(m, n, o, k);
        
        
        Sk = D*T*Q;
        
        
        
        
        
    otherwise
        error('SETTINGS.SYLVESTER_MATRIX_VARIANT is either standard or rearranged')
end

end



