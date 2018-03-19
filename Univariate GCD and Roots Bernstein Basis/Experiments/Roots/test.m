
ex_num = '4';

fx = Examples_Roots_FromCoefficients(ex_num);

gx = Bernstein_Differentiate(fx);

m = GetDegree(fx);
n = GetDegree(gx);

k = 3;


for k = 1:1:min(m,n)
vLambda(k) = GetGeometricMeanMatlabMethod(fx, n - k);
vMu(k) = GetGeometricMeanMatlabMethod(gx, m - k);

end


bool_log = true

if bool_log
   vLambda = log10(vLambda)
   vMu = log10(vMu)
end

figure()
hold on
plot(vLambda)
plot(vMu)
%plot(vLambda./vMu,'-s')
hold off
