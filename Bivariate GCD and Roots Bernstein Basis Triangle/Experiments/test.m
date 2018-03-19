close all;
clc; 

ex_num = '18';
el = 1e-8;
eu = 1e-10;

mean_method = 'Geometric Mean Matlab Method';
boolAlphaTheta = true;

arrFormat = {'T', 'DT', 'TQ', 'DTQ', 'DTQ Denominator Removed'};
nMethods = length(arrFormat);

for i = 1:1:nMethods
   
    subresultant_format = arrFormat{i};
    o_gcd_Bivariate_2Polys(ex_num, el, eu, mean_method, boolAlphaTheta, 'None', 'None', subresultant_format)


end