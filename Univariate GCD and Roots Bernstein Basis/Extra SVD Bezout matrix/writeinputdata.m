
function [] = writeinputdata(d)

% This function writes the matrix d that defines a polynomial to the
% screen.  

disp (' ');
disp ('    exact root      multiplicity');
disp ('--------------------------------');
for i=1:size(d,1)
    fprintf('% 16.8e %10.0f\n',d(i,:));
end
disp (' ');