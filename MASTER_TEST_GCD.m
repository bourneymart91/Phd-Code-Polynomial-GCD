
home_folder = 'C:\Users\Martin\Dropbox (Personal)\PhD\GitHub\PhD-Code';
univar_pb = 'C:\Users\Martin\Dropbox (Personal)\PhD\GitHub\PhD-Code\Univariate GCD and Roots Power Basis';
univar_bb = 'C:\Users\Martin\Dropbox (Personal)\PhD\GitHub\PhD-Code\Univariate GCD and Roots Bernstein Basis';
bivar_pb = 'C:\Users\Martin\Dropbox (Personal)\PhD\GitHub\PhD-Code\Bivariate GCD and Roots Power Basis';
bivar_bb_rect = 'C:\Users\Martin\Dropbox (Personal)\PhD\GitHub\PhD-Code\Bivariate GCD and Roots Bernstein Basis';
bivar_bb_tri =  'C:\Users\Martin\Dropbox (Personal)\PhD\GitHub\PhD-Code\Bivariate GCD and Roots Bernstein Basis Triangle';

% Task 1.
cd(univar_pb)
close all; clc; 
o_gcd_Univariate_2Polys('1',1e-10,1e-12,'Geometric Mean Matlab Method', 'y','Standard STLN','None')

% Task 2.
cd(univar_bb)
close all; clc; 
o_gcd_Univariate_2Polys('1',1e-10,1e-12,'Geometric Mean Matlab Method', 'y','Standard STLN','Standard APF Nonlinear','DTQ')

% Task 3.
cd(bivar_pb)
close all; clc; 
o_gcd_Bivariate_2Polys('1',1e-10,1e-12,'Geometric Mean Matlab Method', 'y','Standard STLN','None','Relative')

% Task 4.
cd(bivar_bb_rect)
close all; clc; 
o_gcd_Bivariate_2Polys('1',1e-10,1e-12,'Geometric Mean Matlab Method','y','Standard STLN','None','DTQ')

% Task 5
cd(bivar_bb_tri)
close all; clc; 
 o_gcd_2Polys('1',1e-10,1e-12,'Geometric Mean Matlab Method','y','Standard STLN','None','DTQ')