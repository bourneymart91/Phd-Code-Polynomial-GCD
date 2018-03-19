function [] = Batch_Experiment1_SylvesterFormats_2Polys(bool_preproc)

arrExampleNumber = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14'};


for i = 1: 1 : length(arrExampleNumber)
    
    ex_num = arrExampleNumber{i};
    
    Experiment1SylvesterFormat_2Polys(ex_num, bool_preproc)
    
end

end