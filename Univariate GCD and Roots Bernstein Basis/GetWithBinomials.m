function fx_bi = GetWithBinomials(fx)
% Given the polynomial coefficients f(x) in the Bernstein basis, multiply
% the coefficients a_{i} by nchoosek(m,i), to obtain coefficients in the
% scaled Bernstein form.
%
%
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% % Outputs
%
% fx_bi : (Vector) Coefficients of the polynomial f(x) in scaled Bernstein
% form.


% Get the degree of polynomial f(x).
m = GetDegree(fx);

% Get the binomial coefficients corresponding to f(x)
bi_m = GetBinomials(m);

% Multiply coefficients of f(x) by the binomial coefficients
fx_bi = fx.*bi_m;


end