clear
readLength = 50;

mp = ['ACTG'];
for i=1:1000
  Sequence1{i} = mp(unidrnd(4,1,1500));
end

%%%%%%55
% create reads

ind_bac_in_mix = [1,2];
correctWeight = zeros(1,1000);
correctWeight(1:2) = [0.1,0.9];
[fracRelevantReadsForInfinity,junk,values,junk1]=prepareReadsForRandomData(readLength,Sequence1,ind_bac_in_mix,correctWeight);

% test a group that does not include 2
[junk1,normalizedBac,uniqueReads,uniqueReads_length]=prepareReadsForRandomData(readLength,Sequence1([1,3:end]));

[junk,i1,i2] = intersect(values,uniqueReads,'rows');
y = fracRelevantReadsForInfinity(i1);

[currResNoNorm,firstFlag]=runOneGroupOf1000post(1,normalizedBac,y,0);


plot(currResNoNorm,'.')

% namely - we get the correctweight of the single bacteria which appears in the currect 1000 

