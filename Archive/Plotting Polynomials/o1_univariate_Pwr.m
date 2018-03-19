function [] = o1_univariate_Pwr(ex_num)
% Given a set of polynomial roots. Plot the polynomial. in univariate power
% basis.

example_style = input('Enter type of problem: \n roots (r) or coefficients (c) ','s')

switch example_style
    case 'r'
        
        example_num = input('Enter Example Number : ','s')
        [roots,~] = Examples_Univariate_Roots(example_num)
        
        % Obtain the coefficients of the power basis polynomial
        f = BuildPoly_Pwr(roots);
        
    case 'c'
        % Obtain the coefficients of the power basis polynomial
        example_num = input('Enter Example Number : ','s')

        
        f = [2 20 10 7 36 5 9];

        
end


% Establish an interval and step size.

lwr_bound = -10;
upr_bound =10;
intvl = 0.001;

x_vec = lwr_bound:intvl:upr_bound;
% Initialise the y component for polynomial f
y_vec_f = zeros(size(x_vec,1));
% Initialise the y component for polynomial g
y_vec_g = zeros(size(x_vec,1));


% Get values f(x) for a series of x values over the defined interval

% For each entry in the x vector
for i = 1:1:length(x_vec)
    x_val = x_vec(i);
    y_vec_f(i) = Evaluate_PowerPoly_Univariate(x_val,f);
   
end


%% Plot graph of the polynomial
figure(1)
p1 = plot(x_vec,(y_vec_f),'-r','LineWidth',2,'DisplayName','f(y)');
hold on
axis([lwr_bound,upr_bound,-100,+100])
legend([p1])



end



