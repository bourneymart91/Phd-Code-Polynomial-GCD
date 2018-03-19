function [] = PrintOutRoots(wx,x)
%Print out the set of w_{i} in term of x

fprintf('----------------------------------------------------------------\n')
fprintf('Printing the roots with respect to %s \n',x)
[~,c] = size(wx);

for i = 1:1:c
    fprintf('The factors of multiplicity %i: \n',i)
    wxi = cell2mat(wx(i));
    wxi./wxi(1,end)
end

end