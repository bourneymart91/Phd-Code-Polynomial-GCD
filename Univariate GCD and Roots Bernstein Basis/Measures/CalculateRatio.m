
function [value] = CalculateRatio(Ak)

% This function calculates the residual r of an approximate linear
% algebraic equation whose coefficient matrix is ds and right hand 
% side vector is line. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    

         [~,R] = qr(Ak);
         R = abs(R);
         [R1_rows,~] = size(diag(R));
         R1 = R(1:R1_rows,1:R1_rows);
         
         
         R_sum_rows = norm(R1,2);
         
         R_diag = abs(diag(R));
         
%          N = zeros(size(R1));
%          
%          N = sqrt(sum(abs(R1).^2,2))
%          N = sort(N,'descend')
         
         
         value = 1/min(R_sum_rows);
         
%         value1 = 1/min(R_diag);
%         value2 = 1/min(R_sum_rows);
%          
   
%         figure(999) 
%             x = 1:size(R1,1);
%             y = N(x);
%             plot(x,y)
end

