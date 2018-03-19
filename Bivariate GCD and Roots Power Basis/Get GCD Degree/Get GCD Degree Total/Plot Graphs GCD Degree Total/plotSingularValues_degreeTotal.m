function plotSingularValues_degreeTotal(arr_SingularValues, limits_k, limits_t)
% 
% % Inputs
%
% arr_SingularValues : (Array of Vectors) Each vector contains singular
% values of S_{k}(f,g)
%
% limits_k : [(Int) (Int)]
%
% limits_t : [(Int) (Int)]

%
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);


nSubresultants = upperLimit_k - lowerLimit_k + 1;

figure_name = sprintf([mfilename ' : ' 'Plotting Singular Values']);
figure('name', figure_name)
hold on


for i = 1 : 1 : nSubresultants
   
    k = lowerLimit_k + (i-1);
    
    % Get vector of singular values of S_{k}
    vSingularValues = arr_SingularValues{i};
    
    v_ks = k.*ones(length(vSingularValues),1);
    
    % plot singular values
    plot(v_ks, log10(vSingularValues),'*')
    
end

xlim(limits_k);
vline(limits_t, {'-r','-r'});


hold off


end