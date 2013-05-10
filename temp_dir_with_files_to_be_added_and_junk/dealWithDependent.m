function [x_extra,ind_extra]=dealWithDependent(dependentInd,auxData,basicAllocationFileName,basicSeqNameDir,basicSeqKey,readLength,uniqueReads,uniqueReads_length)

if exist('basicAllocationFileName') && ~isempty(basicAllocationFileName) % deployed
  load(basicAllocationFileName)
else
  disp('using user input in dealWithDependent')
end

dep = [];
for i=1:length(dependentInd)
  dep = [dep,dependentInd{i}];
end

clear currDep
iterDepCounter = 1;
currDep{iterDepCounter} = dep;
formerLength = inf;
%keyboard
while length(currDep{iterDepCounter})>auxData.minimalNumberOfDependent & length(currDep{iterDepCounter})<formerLength
  
  formerLength = length(currDep{iterDepCounter});
  clear res_extra
  
  [res_extra.x,res_extra.sumRelevantReads,res_extra.independentInd,res_extra.dependentInd]=solveForGroup4(1,1,auxData,[],currDep(iterDepCounter),basicSeqNameDir,basicSeqKey,readLength,uniqueReads,uniqueReads_length);
  
  tmp_x_extra{iterDepCounter} = res_extra.x{1};
  tmp_ind_extra{iterDepCounter} = res_extra.independentInd{1};
  
  iterDepCounter = iterDepCounter+1;
  currDep{iterDepCounter} = res_extra.dependentInd{1};
  
  clear res_extra
  disp('the extra is considered a regular bucket - with the same threshold! it it ok?')
  disp(['size of group: ',num2str(length(currDep{iterDepCounter}))])
end
%keyboard

if length(currDep{iterDepCounter})==formerLength
  currDep(iterDepCounter) = [];
  tmp_x_extra(iterDepCounter-1) = [];
  tmp_ind_extra(iterDepCounter-1) = [];
  disp(['iterations stopped - did not find a small group than: ',num2str(length(currDep{iterDepCounter}))])
end

x_extra = [];
ind_extra = [];

if exist('tmp_x_extra')
  for i=1:length(tmp_x_extra)
    x_extra = [x_extra;tmp_x_extra{i}];
    ind_extra = [ind_extra,tmp_ind_extra{i}];
  end 
end

% add those that were left, and their size is smaller than minimalNumberOfDependent
disp(['adding those whose size is smaller than minimalNumberOfDependent: ',num2str(length(currDep{end}))])
x_extra = [x_extra;zeros(length(currDep{end}),1)];
ind_extra = [ind_extra,currDep{end}];

