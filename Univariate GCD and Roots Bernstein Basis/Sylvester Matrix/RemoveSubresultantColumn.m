function [Ak , ck] = RemoveSubresultantColumn(Sk, i)    
%% Remove the column from the subresultant at index i
%   Sk  =   Subresultant S_{k}
%   i   =   index of column to be removed
    
    ck = Sk(:,i);
    Sk(:,i) = [];
    Ak = Sk;

end

 

