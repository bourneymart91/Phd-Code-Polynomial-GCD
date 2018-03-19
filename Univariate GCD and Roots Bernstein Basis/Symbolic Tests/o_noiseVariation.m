function [] = o_noiseVariation(ex)
%% noise variation: A File to experiment with noise in both Bernstein and Power Bais.
% Given an example number, take the first polynomial from the example file,
% and experiment with various noise levels. Graph the noisy polynomials in
% power and bernstien basis.

% Set lower bound of interval.
a = 0;

% Set upper bound of interval
b = 1.2;

% Set increment for evaluation
inc = 0.01;

% Get the set of polynomial roots.
f_roots = Private_Examples(ex);

% Vary noise on coefficients of polynomial f, in Bernstein basis.
o_noiseVariation_Bernstein(f_roots,a,b,inc);

% Vary noise on coefficients of polynomial f, in Power basis.
o_noiseVariation_Power(f_roots,a,b,inc);

end


function [] = o_noiseVariation_Bernstein(f_root_mult_arr,a,b,inc)
%% Trial varying levels of noise on the coefficients of f in the Bernstein Basis.
% Graph the noisy polynomial, which are perturbed forms of f.

% Get polynomial f in scaled Bernstein Basis.
f_exact_bi = BuildPolyFromRoots(f_root_mult_arr);

% Get degree of polynomial f
m = GetDegree(f_exact_bi);

% Calculate binomial coefficients corresponding to polynomial f.
Bi_m = GetBinomials(m);

% Get polynomial f in standard Bernstein Basis.
f_exact = GetWithoutBinomials(f_exact_bi);

%% Add varying levels of noise to f.
fx10 = AddVariableNoiseToPoly(f_exact,1e-10,1e-10);
fx9  = AddVariableNoiseToPoly(f_exact,1e-9,1e-9);
fx8  = AddVariableNoiseToPoly(f_exact,1e-8,1e-8);
fx2  = AddVariableNoiseToPoly(f_exact,1e-2,1e-2);
fx1  = AddVariableNoiseToPoly(f_exact,1e-1,1e-1);
fx   = f_exact;


% Evaluate the noisy functions at increments over interval [a,b]
f10  = EvaluateFunction_BernsteinBasis(a,b,inc,fx10);
f9 = EvaluateFunction_BernsteinBasis(a,b,inc,fx9);
f8 = EvaluateFunction_BernsteinBasis(a,b,inc,fx8);
f2 = EvaluateFunction_BernsteinBasis(a,b,inc,fx2);
f1 = EvaluateFunction_BernsteinBasis(a,b,inc,fx1);
f  = EvaluateFunction_BernsteinBasis(a,b,inc,fx);


x = a:inc:b;

% Graph the noisy functions in bernstein basis
figure(100)
grid on
hold on
title('Bernstein Basis')
plot(x,f10,'-','DisplayName','f(x) 1e-10');
plot(x,f9,'-','DisplayName','f(x) 1e-9');
plot(x,f8,'-','DisplayName','f(x) 1e-8');
plot(x,f2,'-','DisplayName','f(x) 1e-2');
plot(x,f1,'-','DisplayName','f(x) 1e-1');
plot(x,f,'-r','DisplayName','f(x)');
hold off
legend(gca,'show')
hold off

end


function [] = o_noiseVariation_Power(f_roots,a,b,inc)
%% Trial varying levels of noise on the coefficients of f in the Power Basis.
% Graph the noisy polynomial, which are perturbed forms of f.

% Flatten the roots. Multplicities are discarded. eg [r1 4] becomes [r1,r1,r1,r1]
my_roots = zeros(1,sum(f_roots(:,2)));

% initialise counter
k = 0;

% for each unique root
for i = 1:1:size(f_roots,1)
    % for each multiplicity of the given root.
    for j = 1:1:f_roots(i,2)
        k = k+1;
        % place the root into the vector of roots.
        my_roots(k) = f_roots(i,1);
    end
end

% convert my_roots to polynomial coefficients.
f_exact = poly(my_roots);

% Get length of polynomial.
m = GetDegree(f_exact);

% Add varying levels of noise to the polynomial coefficients.
fx10 = VariableNoise(f_exact,1e-10,1e-10);
fx9  = VariableNoise(f_exact,1e-9,1e-9);
fx8  = VariableNoise(f_exact,1e-8,1e-8);
fx2  = VariableNoise(f_exact,1e-2,1e-2);
fx1  = VariableNoise(f_exact,1e-1,1e-1);
fx = f_exact;


% Evaluate the polynomials at a set of points over interval [a,b]
f10 = EvaluateFunction_PowerBasis(a,b,inc,fx10);
f9  = EvaluateFunction_PowerBasis(a,b,inc,fx9);
f8  = EvaluateFunction_PowerBasis(a,b,inc,fx8);
f2  = EvaluateFunction_PowerBasis(a,b,inc,fx2);
f1  = EvaluateFunction_PowerBasis(a,b,inc,fx1);
f   = EvaluateFunction_PowerBasis(a,b,inc,fx);

% Graph the noisy functions in power basis
figure(101)
grid on
hold on
title('power basis')
x = a:inc:b;
plot(x,f10,'-','DisplayName','f(x) 1e-10');
plot(x,f9,'-','DisplayName','f(x) 1e-9');
plot(x,f8,'-','DisplayName','f(x) 1e-8');
plot(x,f2,'-','DisplayName','f(x) 1e-2');
plot(x,f1,'-','DisplayName','f(x) 1e-1');
plot(x,f,'-r','DisplayName','f(x) exact');

hold off
legend(gca,'show')
hold off

end


function a = Private_Examples(ex_num)
%% Set of private examples for the function NoiseVariation().

switch ex_num
    
    case 0
        a = [
            0.1 1
            0.3 1
            0.5 1
            0.7 1
        ];
    case 1
        a = [
         0.7    2   
        ];
    
    case 21
        a = [
            0.9   5
            ];
    case 22 
        a = [
            0.90    40
            0.56    4
            0.69    2
            ];
end


end