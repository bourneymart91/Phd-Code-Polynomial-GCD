function o1_trials()
% Function purely for producing report outputs, wherein the matlab roots
    %function is trialled many times

% Set number of trials
Num_Trials = 1000;

% Set upper and lower noise bounds
emin = 1e-007;
emax = 1e-009;

% Initialise a new figure
figure()
hold on


% Root Methods
% 0 - Matlab Method.
% 1 - MultRoots Method.
RootMethod = 1;


for i=1:1:Num_Trials
    
    % Choose root solving method
    switch RootMethod
        case 0
            % Get roots by matlab method
            roots_1 = o_roots_matlab(1,emin,emax);
            scatter(real(roots_1),imag(roots_1),'DisplayName','Matlab Roots');            
        case 1
            % Get roots by multroot method.
            roots_2 = o_roots_multroot(1,emin,emax);
            scatter(real(roots_2(:,1)),imag(roots_2(:,1)),'DisplayName','Multroots Roots');
    end
    
    % Note using this produces thousands of graphs, graphing should be
    % off before using this call.
    % roots_3 = o_roots_mymethod(1,emin,emax,1,1,1,1,1,0,1);
    % scatter(real(roots_3(:,1)),imag(roots_3(:,1)),'DisplayName','MyRoots Roots');
    
    
end
hold off

end



function [roots_calc] = o_roots_matlab(ex_num,emin,emax)
% Get roots of Bernstein Basis Polynomial by matlab function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ex_num : example number
% emin   : minimum noise/signal \epsilon_{i}
% emax   : maximum noise/signal \epsilon_{i}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get exact roots of polynomial f
f_roots_exact = Private_Examples(ex_num);

% Get coefficients of exact polynomial f in scaled bernstein basis.
f_exact_bi = BuildPolyFromRoots(f_roots_exact);

% Get degree of polynomial f.
m = length(f_exact_bi) - 1;

% Produce a vector of binomial coefficients corresponding to the
% coefficients in polynomial f.
Bi_m = zeros(1,m+1);
for i = 0:1:m
    Bi_m(i+1) = nchoosek(m,i);
end

% Divide the polynomial f in scaled bernstein basis, to obtain in regular
% bernstein basis.
f_exact = f_exact_bi./Bi_m;

% Add noise to the coefficients of polynomial f in the regular bernstein
% basis.
fx = AddVariableNoiseToPoly(f_exact,emin,emax);

% Get roots wrt to y using MATLAB method.
roots_calc = roots(fx.*Bi_m);

% convert roots to bernstein basis
roots_calc = 1- roots_calc./(1+roots_calc);


end

function [roots_calc] = o_roots_multroot(ex_num,emin,emax)
% Get roots of Bernstein Basis Polynomial by mutlroot function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ex_num : example number

% emin   : minimum noise/signal \epsilon_{i}

% emax   : maximum noise/signal \epsilon_{i}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath 'multroot/multroot'

% Get exact roots of polynomial f
f_root_mult_arr_exact = Private_Examples(ex_num);

% Get coefficients of exact polynomial f in scaled bernstein basis.
f_exact_bi = BuildPolyFromRoots(f_root_mult_arr_exact);

% Get degree of polynomial f.
m = length(f_exact_bi) - 1;

% Produce a vector of binomial coefficients corresponding to the
% coefficients in polynomial f.
Bi_m = zeros(1,m+1);
for i = 0:1:m
    Bi_m(i+1) = nchoosek(m,i);
end

% Divide the polynomial f in scaled bernstein basis, to obtain in regular
% bernstein basis.
f_exact = f_exact_bi./Bi_m;

% Add noise to the coefficients of polynomial f in the regular bernstein
% basis.
fx = AddNoiseToPoly(f_exact,emin);

% Get roots wrt to y
roots_calc = multroot(fx.*Bi_m);

% convert roots wrt t, Bernstein basis
roots_calc = [1- roots_calc(:,1)./(1+roots_calc(:,1)) roots_calc(:,2)];


end




function a = Private_Examples(ex_num)

switch ex_num
    case 1 
        a = [
            0.1     10
            -0.1    10
        ];
    
end


end


