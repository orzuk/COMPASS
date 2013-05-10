function changeNumberOfProcessors(dirName,numProcessors)
% changeNumberOfProcessors('canonicalFull',50)
dr = dir(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/',dirName,'/structFor_*_Correction_1.mat']);
%keyboard
for i=1:length(dr)
  clear auxData
  load(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/',dirName,'/',dr(i).name],'auxData');
  auxData.numProcessors = numProcessors;
  save(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/',dirName,'/',dr(i).name],'auxData','-append')
end
