function [] = plotRowDiag_degreeTotal(arr_R1_diag, limits_k, limits_t)
%
% % Inputs
%
% arr_R1_diag : (Array of Vectors)
%
% limits_k : [Int Int]
%
% limits_t : [Int Int]

myLowerLimit = limits_k(1);
myUpperLimit = limits_k(2);

nSubresultants = myUpperLimit - myLowerLimit + 1;

figure_name = sprintf([mfilename ' : ' 'Plot Row Diagonals of S']);
figure('name',figure_name)

hold on

for i = 1 : 1 : nSubresultants
   
    k = myLowerLimit + (i-1);
    
    vR1_diags = arr_R1_diag{i};
    
    vec_k = k .* ones(length(vR1_diags),1);
    
    plot(vec_k,log10(vR1_diags),'*');
    
end

vline(limits_t, {'-r','-r'});

hold off

end