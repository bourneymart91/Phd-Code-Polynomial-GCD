function [fww] = GetWithoutThetas(fxy,m,th1,th2)

% Multiply the rows of fxy by th1

pre_thetas = diag(1./(th1.^(0:1:m)));
post_thetas = diag(1./(th2.^(0:1:m)));

fww = pre_thetas * fxy * post_thetas;



end