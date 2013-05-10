% Compare original and resulting results
% 
% Input: 
% reconstructed_weights - vector returned by the reconstruction algorithm
% correct_weights - true frequencies 
% comparison_output_file - where to save results 
% 
function CompareSolutionToTrueMixture(reconstructed_weights, correct_weights, comparison_output_file)

plot(reconstructed_weights, correct_weight,'.')
xlabel('found weight');ylabel('correct weight')
title('based on iterations - last iteration of 358  bacteria')
xlabel('found freq')
ylabel('correct freq')

print('-dpdf', comparison_output_file); 
hgsave(comparison_output_file); 




