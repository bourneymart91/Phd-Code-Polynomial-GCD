function [] = o_gcd_experiment_method(fx,gx)



[lambda, mu, alpha, theta] = Preprocess(fx,gx);

fx_n = fx./lambda;
gx_n = gx./mu;

fw = GetWithThetas(fx_n,theta);
gw = GetWithThetas(gx_n,theta);


% Build the first subresultants
S1 = BuildT(fw,alpha.*gw,1);

% Get rank of S1
svd(S1);


% Get optimal column for removal


% Get u(x) and v(x)


end