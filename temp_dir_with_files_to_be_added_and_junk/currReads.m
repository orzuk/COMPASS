function [fracRelevantReads,sumRelevantReads]=currReads(uniqueReads,uniqueReads_length,values,inifiniteNumberOfReadsFlag,dataIn)

% dataIn - is diffenent for the two cases - infinite and finite

if ~exist('inifiniteNumberOfReadsFlag')
  inifiniteNumberOfReadsFlag = 0;
  if exist('dataIn')
    disp('problem. currReads.m')
    pause
  end
  dataIn = struct;
  dataIn.normalizeByTotalNumberOfReadsFlag = 0;
end

%keyboard
okFlag = 0;
pt = size(uniqueReads,1)-1;
while okFlag==0
  try
    part = 1:pt:size(uniqueReads,1);
    part(end) = size(uniqueReads,1)+1;
    
    i1 = [];
    i2 = [];
    for i=1:length(part)-1
      i
      tmpInd = [part(i):part(i+1)-1]';
      [junk,i11,i22] = intersect(values,uniqueReads(tmpInd,:),'rows');
      i1 = [i1;i11];
      i2 = [i2;tmpInd(i22)];
    end
    okFlag = 1;
  catch me 
    pt = round(pt/2)
    disp('reduced size in runOneGroupOf1000')
    okFlag = 0;
  end
end




%keyboard
if inifiniteNumberOfReadsFlag==0

  %keyboard
  if exist('fracRelevantReadsIn')
    disp('problem - fracRelevantReads!!! but it is finite?!. currReads.m')
    pause
  end
  numRelevantReads = zeros(size(values,1),1);
  clear values
  for i=1:length(i1)
    numRelevantReads(i1(i)) = uniqueReads_length(i2(i));
  end
  
  if dataIn.normalizeByTotalNumberOfReadsFlag==0
    disp('normalize by the LOCAL number of reads')
    fracRelevantReads = numRelevantReads/sum(numRelevantReads);
  else
    disp('normalize by the TOTAL number of reads')
    fracRelevantReads = numRelevantReads/sum(uniqueReads_length);
  end
  sumRelevantReads = sum(numRelevantReads);


else % infinite
  if isfield(dataIn,'normalizeByTotalNumberOfReadsFlag')
    disp('problem normalizeByTotalNumberOfReadsFlag. currReads.m');
    pause
  end
  
  fracRelevantReads = zeros(size(values,1),1);
  for i=1:length(i1)
    fracRelevantReads(i1(i)) = dataIn.fracRelevantReadsForInfinity(i2(i));
  end
  sumRelevantReads = sum(fracRelevantReads);
  disp('notice that for infinite number of reads - fracRelevantReads is normalized to the TOTAL number of reads and not to the relevant number of reads!!!. currReads.m')
end
