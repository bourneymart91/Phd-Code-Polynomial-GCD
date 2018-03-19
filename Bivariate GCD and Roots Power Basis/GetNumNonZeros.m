function count = GetNumNonZeros(m1,m2,m)

count = 0;
for i = 0:1:m1
    for j = 0:1:m2
        if i + j <= m
            count = count + 1;
        end
    end
end

% try
%     count2 = (m1+1) * (m2+1) - nchoosek(m1+m2-m+1,2);
% catch
%     count2 = (m1+1) * (m2+1);
% end
% 
% if count ~= count2
%    error('err'); 
% end


end