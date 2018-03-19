function lambda = GetArithmeticMean(fxy, m, n_k)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% m : (Int) Degree of f(x,y)
%
% n_k : (Int)
%
% % Outputs
%
% lambda : (Float) Arithmetic Mean


lambda = ... 
        sum(sum(fxy))...
        * (nchoosek(m+n_k+2,2) ./ ( nchoosek(m+2,2)^2 * nchoosek(n_k+2,2)  ));


    
%Tf = BuildDT1Q1(fxy,m,n_k);
%vec = Tf(Tf~=0);
%lambda = mean(vec);


end

