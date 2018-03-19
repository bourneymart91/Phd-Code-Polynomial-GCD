function [] = PrintAllExamples_GCD()

% Example Numbers
ex_num_arr = {'1','2','3','4','5','6','7','8','9','10','11'};

for i = 1:1:length(ex_num_arr)
    
    ex_num = ex_num_arr{i};
    
    [f_roots_mult_arr,g_roots_mult_arr,d_roots_mult_arr,...
        u_roots_mult_arr,v_roots_mult_arr] = Bivariate_GCD_Examples(ex_num);
    
    
    symbolic_f = GetSymbolicPolyFromSymbolicRoots(f_roots_mult_arr);
    symbolic_g = GetSymbolicPolyFromSymbolicRoots(g_roots_mult_arr);
    symbolic_d = GetSymbolicPolyFromSymbolicRoots(d_roots_mult_arr);
    symbolic_u = GetSymbolicPolyFromSymbolicRoots(u_roots_mult_arr);
    symbolic_v = GetSymbolicPolyFromSymbolicRoots(v_roots_mult_arr);
    
    
    [m,m1,m2] = GetDegree_SymbolicPoly(symbolic_f);
    [n,n1,n2] = GetDegree_SymbolicPoly(symbolic_g);
    [t,t1,t2] = GetDegree_SymbolicPoly(symbolic_d);
    [m_t,m1_t1,m2_t2] = GetDegree_SymbolicPoly(symbolic_u);
    [n_t,n1_t1,n2_t2] = GetDegree_SymbolicPoly(symbolic_v);
    
    fprintf('Example : %s \n',ex_num);
    
    fprintf('\\begin{tabular}{|l|l|l|} \n');
	fprintf('\\hline \n');
	fprintf('	$m = %s $ \n',int2str(m));
	fprintf('& 	$m_{1} = %s$ \n', int2str(m1));
	fprintf('&	$m_{2} = %s $ \n', int2str(m2))
	fprintf('\\\\ \n')
	fprintf('\\hline \n')
		fprintf('$n = %s $ \n', int2str(n))
	fprintf('& 	$n_{1} = %s $ \n', int2str(n1))
	fprintf('& 	$n_{2} = %s $ \n', int2str(n2))
	fprintf('\\\\ \n')
	fprintf('\\hline \n')
		fprintf('$ t = %s $ \n', int2str(t))
	fprintf('&	$t_{1} = %s $ \n', int2str(t1))
	fprintf('& 	$t_{2} = %s $ \n', int2str(t2))
	fprintf('\\\\ \n')
    fprintf('\\hline \n')
    fprintf('\\end{tabular} \n')


    fprintf('\\begin{align*} \n')
    fprintf('f &= %s \n',(latex(symbolic_f)));
    fprintf('\\\\ \n')
    fprintf('g &= %s \n',(latex(symbolic_g)));
    fprintf('\\\\ \n')
    fprintf('d &= %s \n',(latex(symbolic_d)));
    fprintf('\\\\ \n')
    fprintf('u &= %s \n',(latex(symbolic_u)));
    fprintf('\\\\ \n')
    fprintf('v &= %s \n',(latex(symbolic_v)));
    fprintf('\\\\ \n')
    fprintf('\\end{align*} \n')
end

end

function [m,m1,m2] = GetDegree_SymbolicPoly(symbolic_poly)

syms x y;

m = double(feval(symengine, 'degree', symbolic_poly));
m1 = double(feval(symengine, 'degree', symbolic_poly,x));
m2 = double(feval(symengine, 'degree', symbolic_poly,y));

end