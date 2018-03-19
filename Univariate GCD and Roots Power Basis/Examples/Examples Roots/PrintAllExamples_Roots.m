function [] = PrintAllExamples_Roots()

% Example Numbers
ex_num_arr = {'1','2','3','4','5','6','7','8','9','10','11'};
addpath('../Examples')
for i = 1:1:length(ex_num_arr)
    
    ex_num = ex_num_arr{i};
    
    [f_roots_mult_arr] = Univariate_Roots_Examples(ex_num);
    
    
    symbolic_f = GetSymbolicPolyFromSymbolicRoots(f_roots_mult_arr);
    
    
    [m] = GetDegree_SymbolicPoly(symbolic_f);
    
    fprintf('Example : %s \n',ex_num);
    fprintf('\\\\ \n');
    
    fprintf('\\begin{center} \n');
    fprintf('\\begin{tabular}{|l|} \n');
	fprintf('\\hline \n');
	fprintf('	$m = %s $ \n',int2str(m));
	fprintf('\\\\ \n')
	fprintf('\\hline \n')
    fprintf('\\end{tabular} \n')
    fprintf('\\end{center} \n')

    fprintf('\\begin{align*} \n')
    fprintf('f &= %s \n',(latex(symbolic_f)));
    fprintf('\\\\ \n')
    fprintf('\\end{align*} \n')
end

end

function [m] = GetDegree_SymbolicPoly(symbolic_poly)

m = double(feval(symengine, 'degree', symbolic_poly));


end