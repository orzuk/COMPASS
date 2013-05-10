function [X,sumRelevantReads]=divAndCon(uniqueReads,uniqueReads_length,curr_kp,setParameters)
%keyboard

readLength = setParameters.readLength;
currNumProcessors = setParameters.numProcessors;


n = length(curr_kp);

if setParameters.keepOriginalOrderFlag==0
  currInd = randperm(n);
  disp('using permuted order')
else
  currInd = 1:n;
  disp('using the original order')
end

part = 1:setParameters.groupSize:n;
part(end) = n+1;


tmpInd = cell(length(part)-1,1);
for i=1:length(part)-1
  tmpInd{i} = part(i):part(i+1)-1;
  tmpInd{i} = curr_kp(currInd(tmpInd{i}));
end

firstFlag = 1;

disp('do also parallel')
pause

if currNumProcessors>1 
    %keyboard
    
    % using matlab's parallel tool
    w = ['matlabpool open local ',num2str(currNumProcessors),';parfor i=1:length(tmpInd),i,[xx,ss]=solveForGroupOrFourth(i,i,setParameters,[],[],tmpInd,basicSeqNameDir,basicSeqKey,readLength,uniqueReads,uniqueReads_length);x{i} = xx{i};sumRelevantReads(i) = ss(i);;end;matlabpool close'];
    eval(w);
    
  %keyboard
  else % single
    %keyboard 
    
    for i=1:length(tmpInd)
      [res.x,res.sumRelevantReads]=solveForGroupOrFourth(i,i,setParameters,[],[],tmpInd,basicSeqNameDir,basicSeqKey,readLength,uniqueReads,uniqueReads_length);
      x{i} = res.x{i};
      sumRelevantReads(i) = res.sumRelevantReads(i);
    end
end


%keyboard


if size(x{1},2)==1 % standard
  X = zeros(1,n);
  for i=1:length(part)-1
    tmpInd = part(i):part(i+1)-1;
    tmpInd = currInd(tmpInd);
    if ~isempty(x{i})
      X(tmpInd) = x{i};
    end
  end
else % x_min x_max
  %keyboard
  if length(part)~=2 % should be a single one
    disp('problem. iterateParallelDistributedSequenceFilesOr.m')
    keyobard
  end
  %keyboard
  X = zeros(size(x{i}));
  for i=1:length(part)-1
    tmpInd = part(i):part(i+1)-1;
    tmpInd = currInd(tmpInd);
    if ~isempty(x{i})
      X(tmpInd,:) = x{i};
    end
  end
  
end




