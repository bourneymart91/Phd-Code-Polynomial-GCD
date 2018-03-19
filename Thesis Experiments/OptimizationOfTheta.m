function [] = OptimizationOfTheta()
% This experiment is designed to show how the inclusion of theta scales the
% coefficients of the polynomials.

fx = ...
    [
    143.545; 
    163.154; 
    202.5789654;
    10.123748;
    42.546987;
    10.246579;
    2500.064
    5454;
    154.3547;
    5676.457
    ];

hx = ...
    [
        15646;
        5647;
        12347;
        6872;
        1657;
        1657;
        13657;
        6657;
        14;
        124.54;
        154.45
    ];


fx = abs(fx);

m = GetDegree(fx);
nplusm = GetDegree(hx);
n = nplusm - m;
f = [1 -1 0];

PartOne = ...
    [
    ones(m+1,1) , zeros(m+1,1), -1.*(0:1:m)'
    ];

PartTwo = ...
    [
    zeros(m+1,1), -1.*ones(m+1,1), (0:1:m)'
    ];

A = ...
    [
    PartOne;
    PartTwo
    ];

b = [log10(abs(fx)); -log10(abs(fx))];

x = linprog(f,-A,-b);
theta = 10.^(x(3));

fw = GetWithThetas(fx,theta);

figure('name','Plotting coefficients log10')
hold on
plot((1:1:m+1),log10(fw),'DisplayName','f(w)')
plot((1:1:m+1),log10(fx),'DisplayName','f(x)')
legend(gca,'show');
hold off

figure('name','Plotting coefficients')
hold on
plot((1:1:m+1),(fw),'DisplayName','f(w)')
plot((1:1:m+1),(fx),'DisplayName','f(x)')
legend(gca,'show');
hold off


% One 
hw = GetWithThetas(hx,theta);
Tfw = BuildT1(fw,n);
gw = SolveAx_b(Tfw,hw)

residual = norm(Tfw*gw - hw)

% Two
gx_from_gw = GetWithoutThetas(gw,theta);
Tfx = BuildT1(fx,n);
residual2 = norm(Tfx*gx_from_gw - hx)

% Three
gx = SolveAx_b(Tfx,hx);
residual3 = norm(Tfx*gx - hx)

cond(Tfw)
cond(Tfx)

end
