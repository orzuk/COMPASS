function freqset=compareResults(dist,resFile,correctFile,basicSeqKey,basicSeqNameDir)
  
load(resFile)
load(correctFile,'correctWeight','ind_bac_in_mix')

disp(['using the last iteration which consists of ',num2str(length(store_kp{end})),' bacteria. If that is too small - use preceding iterations. see compareResults.m'])
disp('thresholds the found frequency at 10^-3. see compareResults.m')

solution = abs(found{end});
foundBAC = find(abs(found{end})>10^-3)';

b = [ind_bac_in_mix;foundBAC];
b = unique(b);
  
readLength = 50;
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFiles2(readLength,b,basicSeqNameDir,basicSeqKey); 

[junk,i_f1] = ismember(foundBAC,b);
[junk,i_c1] = ismember(ind_bac_in_mix,b);

%keyboard
[freqset]=Compare2(normalizedBac,i_c1,correctWeight(ind_bac_in_mix),i_f1,solution(foundBAC),dist);

