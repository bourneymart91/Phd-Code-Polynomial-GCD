function [] = plotRowNormsR1_degreeRelative(arr_R1_RowNorms, limits_k1, limits_k2, limits_t1, limits_t2)
%
% % Inputs
%
% arr_R1_RowNorms : (Array)
%
% myLimits_t1 : (Int Int)
%
% myLimits_t2 : (Int Int)
%
% limits_t1 : (Int Int)
%
% limits_t2 : (Int Int)


lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);
lowerLimit_k2 = limits_k2(1);
upperLimit_k2 = limits_k2(2);

% Get number of subresultant matrices
nSubresultants_k1 = upperLimit_k1 - lowerLimit_k1 + 1;
nSubresultants_k2 = upperLimit_k2 - lowerLimit_k2 + 1;

figure_name = sprintf([mfilename ' : ' 'Plot Row Norms of R1 from QR decomposition']);
figure('name',figure_name)
hold on

for i1 = 1:1:nSubresultants_k1
    for i2 = 1:1:nSubresultants_k2
   
        % Get k1 and k2
        k1 = lowerLimit_k1 + (i1-1);
        k2 = lowerLimit_k2 + (i2-1);
        
        % Get vector of row norms of R_{1} from S_{k1,k2}
        temp_vec = arr_R1_RowNorms{i1,i2};
        
        v_k1 = k1.* ones(length(temp_vec),1);
        v_k2 = k2.* ones(length(temp_vec),1);
        
        plot3(v_k1, v_k2, log10(temp_vec));
        
    end
end

hold off

hold off





end