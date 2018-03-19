function [] = PrintoutRoots(type,root_mult_arr)

fprintf('\nROOTS CALCULATED %s \n',type);
fprintf('\t Root \t \t \t \t\t \t \t \t \t \t \t \t Multiplicity \n')
fprintf('%22.15f + %22.15f i  \t \t %3g \n', [real(root_mult_arr(:,1)),imag(root_mult_arr(:,1)),...
    root_mult_arr(:,2)]');
fprintf('\n');