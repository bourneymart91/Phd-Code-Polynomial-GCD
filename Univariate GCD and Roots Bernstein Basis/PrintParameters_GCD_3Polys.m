
% Print the parameters.
LineBreakLarge()
fprintf('INPUT VARIABLES\n');
fprintf('\t EXAMPLE NUMBER : %s \n', ex_num);
fprintf('\t EXAMPLE NUMBER VARIANT: %s \n', ex_num_variant);
fprintf('\t EMIN : %s \n' , num2str(SETTINGS.EMIN));
fprintf('\t EMAX : %s \n' , num2str(SETTINGS.EMAX));
fprintf('\t MEAN METHOD : %s \n', SETTINGS.MEAN_METHOD);
fprintf('\t ALPHA_THETA : %s \n', num2str(SETTINGS.BOOL_ALPHA_THETA));
fprintf('\t LOW RANK APPROX METHOD : %s \n', SETTINGS.LOW_RANK_APPROXIMATION_METHOD);
fprintf('\t APF METHOD : %s \n ', SETTINGS.APF_METHOD);
fprintf('\t LOG: %s \n', SETTINGS.BOOL_LOG);
fprintf('\t SYLVESTER MATRIX VARIANT : %s \n', SETTINGS.SYLVESTER_MATRIX_VARIANT);
fprintf('\t SYLVESTER n EQUATIONS: %s \n', SETTINGS.SYLVESTER_EQUATIONS);
fprintf('\t SCALING METHOD: %s \n', SETTINGS.SCALING_METHOD);
LineBreakLarge()