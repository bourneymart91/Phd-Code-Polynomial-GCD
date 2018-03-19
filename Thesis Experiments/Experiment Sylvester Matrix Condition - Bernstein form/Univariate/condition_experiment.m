function [] = condition_experiment(m, n)
%
% % Inputs
%
% m : (Int) 
%
% n : (Int)
%






fx = ones(m+1, 1);
gx = ones(n+1, 1);
arrSylvesterBuildMethod = {'DTQ', 'TQ', 'DT'};
nMethods = length(arrSylvesterBuildMethod);

arr_condition = cell(nMethods + 1,  1);

for i = 1 : 1 : nMethods
   
    method_name = arrSylvesterBuildMethod{i};
    arr_condition{i} = GetConditionMatrix(fx, gx, method_name);
    
end




figure()
hold on

for i = 1: 1: nMethods

   
    method_name = arrSylvesterBuildMethod{i};
    
    %
    vCondition = arr_condition{i};
    plot(log10(vCondition),'-s','DisplayName', method_name)
    
end
legend(gca,'show');
hold off

end

function vCondition = GetConditionMatrix(fx, gx, method_name)
%
% Get the condition of each of the Sylvester subresultant matrices

m = GetDegree(fx);
n = GetDegree(gx);


vCondition = zeros(min(m,n) + 1, 1);

for k = 1 : 1 : min(m,n)

    Sk = BuildSubresultant(fx, gx, k, method_name);
    vCondition(k+1) =  cond(Sk);

end


end





function Sk = BuildSubresultant(fx,gx,k,subresultant_type)

m = GetDegree(fx);
n = GetDegree(gx);



switch subresultant_type 
    case 'DTQ'
        
        DT1Q1 = BuildDT1Q1(fx,n-k);
        DT2Q2 = BuildDT1Q1(gx,m-k);
        Sk = [DT1Q1 DT2Q2];
        
    case 'TQ'
        
        T1 = BuildT1(fx,n-k);
        T2 = BuildT1(gx,m-k);
        Q1 = BuildQ1(n-k);
        Q2 = BuildQ1(m-k);
        Sk = [T1*Q1 T2*Q2];
        
    case 'DT'
        
        D = BuildD(m,n-k);
        T1 = BuildT1(fx,n-k);
        T2 = BuildT1(gx,m-k);
        Sk = D*[T1 T2];
        
    case 'T'
        T1 = BuildT1(fx,n-k);
        T2 = BuildT1(gx,m-k);
        Sk = [T1 T2];
end

end


function DT1Q1 = BuildDT1Q1(fx,n_k)

m = GetDegree(fx);
D = BuildD(m,n_k);
T1 = BuildT1(fx,n_k);
Q1 = BuildQ1(n_k);

DT1Q1 = D*T1*Q1;

end


function D = BuildD(m,n_k)

D = diag(1./GetBinomials(m+n_k));

end


function T1 = BuildT1(f,n_k)

m = GetDegree(f);

f_bi = f .* GetBinomials(m);


T1 = zeros(m+n_k+1,n_k+1);


for j = 1:1:n_k+1
    T1(j:j+m,j) = f_bi;
end

end



function Q1 = BuildQ1(n_k)
% Build the matrix Q_{1}    

Q1 = diag(GetBinomials(n_k));
end


function [binoms] = GetBinomials(m)
binoms = zeros(m+1,1);
for i = 0:1:m
    binoms(i+1) = nchoosek(m,i);
end
end

function [m] = GetDegree(fx)

m = length(fx) - 1;

end

