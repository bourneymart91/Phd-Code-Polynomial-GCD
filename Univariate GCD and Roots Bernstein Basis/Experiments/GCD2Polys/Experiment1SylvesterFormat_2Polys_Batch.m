function [] = Experiment1SylvesterFormat_2Polys_Batch()
% Perform a batch of examples

% Initialise array of example numbers
arrExamples = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14'};

% Get number of examples
nExamples = length(arrExamples);

% For each example
for i = 1 : 1 : nExamples
    
    ex_num = arrExamples{i};
    
    Experiment1SylvesterFormat_2Polys(ex_num)
    
end

end
