function [] = PrintAllExamples_Roots()

% Example Numbers
ex_num_arr = {'1','2','3'};

for i = 1:1:length(ex_num_arr)
    
    ex_num = ex_num_arr{i};
    
    [f_roots_mult_arr] = Bivariate_Roots_Examples(ex_num);
    
    
    symbolic_f = GetSymbolicPolyFromSymbolicRoots(f_roots_mult_arr);
    
    
    
    [m,m1,m2] = GetDegree_SymbolicPoly(symbolic_f);
   
    
    fprintf('Example : %s \n',ex_num);
    
    fprintf('\\begin{center} \n');
    fprintf('\\begin{tabular}{|l|l|l|} \n');
	fprintf('\\hline \n');
	fprintf('	$m = %s $ \n',int2str(m));
	fprintf('& 	$m_{1} = %s$ \n', int2str(m1));
	fprintf('&	$m_{2} = %s $ \n', int2str(m2))
	fprintf('\\\\ \n')
    fprintf('\\hline \n')
    fprintf('\\end{tabular} \n')
    fprintf('\\end{center} \n');

    fprintf('\\begin{align*} \n')
    fprintf('f &= %s \n',(latex(symbolic_f)));
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