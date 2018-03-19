function fxy = GetAsMatrix_Version2(v_fxy, m1, m2)
% Given a vector of coefficients of the polynomial f(x,y), return the
% matrix of coefficients.
%
% % Inputs
%
% v_fxy : (Vector) Coefficients of the polynomial f(x,y) as a vector
%
% m1 : (Int) Degree of f(x,y) with respect to x
%
% m2 : (Int) Degree of f(x,y) with respect to y
%
% % Outputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)

% Initialise matrix
fxy = zeros(m1+1,m2+1);

for j = 1:1:(m2+1)

    % Get start and end point (in the vector vfxy) of the coefficients for 
    % the jth column
    start_index = (j-1) * (m1+1) +1;
    end_index = (j) * (m1+1);
    
    % Get coefficients for the jth column
    temp_vec = v_fxy(start_index: end_index);
    
    % Insert column of coefficients 
    fxy(:,j) = temp_vec;
    
end

end