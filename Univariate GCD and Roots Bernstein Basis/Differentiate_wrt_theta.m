function Partial_fw_wrt_theta = Differentiate_wrt_theta(fw,theta)
% Differentiate f(\omega,\theta) with respect to theta
%
% % Inputs
%
% fw : (Vector) Coefficients of the polynomial f(\omega)
%
% theta : (Float) Optimal value of \theta
%
%
% % Outputs
%
% Partial_fw_wrt_theta : (Vector) Coefficients of the partial derivative of
% f(\omega) with respect to \theta.


bool_diff_method = 'power basis';
%bool_diff_method = 'Bernstein basis';

switch bool_diff_method
    case 'power basis'
        
        Partial_fw_wrt_theta = Differentiate_wrt_theta_powerbasis(fw,theta);
        
    case 'Bernstein basis'
        
        Partial_fw_wrt_theta = Differentiate_wrt_theta_Bernsteinbasis(fw,theta);
        
    otherwise
        
        error('err')
        
end


end