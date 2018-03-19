function [] = PrintRoots(roots_calc,type)

fprintf('\nROOTS CALCULATED %s \n',type);
fprintf('\t Root \t \t \t \t\t \t \t \t \t \t \t \t Multiplicity \n')
fprintf('%22.15f + %22.15f i  \t \t %3g \n',[real(roots_calc(:,1)),imag(roots_calc(:,1)),...
    roots_calc(:,2)]');
fprintf('\n');