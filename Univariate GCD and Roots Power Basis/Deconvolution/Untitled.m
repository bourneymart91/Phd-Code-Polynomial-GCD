tmp_data = importdata('Test_Deconvolution.txt',' ')

mat = tmp_data.data;
error_separate = mat(:,2);
error_batch = mat(:,3);
error_batch_STLN = mat(:,4);
error_batch_constrained = mat(:,5);
error_batch_constrained_STLN = mat(:,6);

figure()
hold on
plot(log10(error_separate),'-s','DisplayName','Separate')
plot(log10(error_batch),'DisplayName','Batch')
plot(log10(error_batch_STLN),'DisplayName','Batch STLN')
plot(log10(error_batch_constrained),'DisplayName','Batch Constrained')
plot(log10(error_batch_constrained_STLN),'-s','DisplayName','Batch Constrained STLN')
legend(gca,'show')
hold off