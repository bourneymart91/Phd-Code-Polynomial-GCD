function gx = Bernstein_Integrate(fx)
% Integrate a polynomial f(x) given in the Bernstein Basis.
%
% Inputs.
%
% fx : Coefficients of input polynomial f(x)
%
%
% Outputs.
%
% gx : Coefficients of integral of f.


% Get Degree of input polynomial.
n = size(fx,1) - 1;

% Initialise vector of integrated polynomial coefficients.
gx = zeros(1,n+1);

% Loop through each coefficient

for k = 0:1:n+1
    i = k+1;
    if k > 0
        I = (1/(n+1)) * ...
            sum(fx(1:(i-1)));
    else
        I = 0;
    end
    gx(i) = I;
end


end