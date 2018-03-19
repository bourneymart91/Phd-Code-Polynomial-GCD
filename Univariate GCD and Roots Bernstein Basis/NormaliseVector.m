function [vNormalised] = NormaliseVector(vector)


vNormalised = vector ./ norm(vector);


end