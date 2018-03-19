function [c] = Bernstein_DegreeElevate(fx,r)
% Function performs degree elevation on polynomial f(x) r times.
%
% Inputs.
%
% fx : The coefficients of polynomial f(x).
%
% r  : Number of degree elevations such that the output polynomial is of
%      degree m + r.

% Get the degree of polynomial f(x)
m = GetDegree(fx);

c = zeros(m+r+1,1);

% Each coefficient c(k+1)
for k = 0:1:m+r
   % Initialise temporary sum.
   temp_sum = 0;
   for j = max(0,k-r):1:min(m,k)
       temp_sum = temp_sum + (fx(j+1) * nchoosek(m,j) * nchoosek(r,k-j)) ./ nchoosek(m+r,k) ;
   end
   c(k+1) = temp_sum;
end


end