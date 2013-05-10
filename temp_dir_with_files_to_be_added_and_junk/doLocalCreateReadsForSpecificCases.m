function doLocalCreateReadsForSpecificCases(auxData,userDir)
disp('does nr>10^5 and read length =400')
%%%%%%%%%%%%%%%%%%%%%555
% create the reads for Nreads>10^5 and readLength>=100
a = find(auxData.nr>10^5 & auxData.nr<Inf);
k = 1;
for i=1:length(a)
  dr_tmp = dir([userDir,'/CS/tmp/',auxData.outFileDirName,'/structFor_',auxData.outFileDirName,'*_readlen_400_*_Nreads_',num2str(auxData.nr(a(i))),'*_createReads*']);
  for j=1:length(dr_tmp)
    dr(k).name = dr_tmp(j).name;
    k = k+1;
  end
end


for i=1:length(dr)
  distributeBAC_generalOrFourth([userDir,'/CS/tmp/',auxData.outFileDirName,'/',dr(i).name]);
end
