function profileMixture(setParameters)
% randomize the seed - change 
rand('seed',sum(100*clock));


%keyboard
clear store* 
k = 1;
curr_kp = 1:setParameters.numBACtoConsider; % remaining bacteria
while length(curr_kp)>setParameters.smallestSetOfCollected 
  %keyboard
    
    if length(curr_kp)>setParameters.upperLimitForRepeat
      % regular
      disp('running one random partition')
      [currX,currSumRelevantReads] = divAndCon(uniqueReads,uniqueReads_length,curr_kp,setParameters);
    else
      [currX,currSumRelevantReads] = doRepeatsFourth(uniqueReads,uniqueReads_length,curr_kp,parallelType,basicSaveName,k,userDir,runFileName,setParameters);
      %what about the value of k?
    end
    

  
  end
  
  next_kp = curr_kp(find(currX>setParameters.thresholdForCollectingBAC));
  store_kp{k} = curr_kp;
  store_X{k} = currX;
  store_sumRelevantReads{k} = currSumRelevantReads;
  k = k +1;
  
  curr_kp = next_kp;
  
  if length(store_kp{k-1})==length(curr_kp)
    disp('did not reduce size, so breaks. distributeBAC_general3.m')
    save([setParameters.saveName,'_break'],'store_kp','store_X','store_sumRelevantReads')
    break
  end
end
if length(store_kp{k-1})~=length(curr_kp) % size was reduced
  store_kp{k} = curr_kp;
  store_X{k} = [];
end

%keyboard
disp('does it take the last one ok? distributeBAC_general3.m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve again
lengths_kp = zeros(1,length(store_kp));
for i=1:length(lengths_kp)
  lengths_kp(i) = length(store_kp{i});
end
%1



%keyboard
do = find(lengths_kp<setParameters.smallestSetOfCollected*2);

% change 12.12.11
if isempty(do) && lengths_kp(end)<3000
  do = length(lengths_kp);
end
% end change 12.12.11
%keyboard
% deal with the case that n-1 is larger than setParameters.smallestSetOfCollected*2 and the n'th is really small

if setParameters.realReadsFromFileFlag==0
  if isempty(find(lengths_kp(do)>100))
    do = do(1)-1;
    disp(['running with larger group than ',num2str(setParameters.smallestSetOfCollected*2),' of size: ',num2str(lengths_kp(do)),' distributeBAC_general3.m'])
  end
else % real reads
   %real reads - take all do
  
end

% solve - regular but the negative solution are the ones which are unique_inds
setParameters.solveOr = 1;
setParameters.do_lp_flag = 0;
setParameters.numProcessors = 1;
setParameters.moveDependentToNextStageFlag = 0;
%setParameters.keepOriginalOrderFlag==1;
for i=do
  curr_kp = store_kp{i};
  setParameters.groupSize = length(curr_kp)-1;
  [currX,currSumRelevantReads] = divAndCon(uniqueReads,uniqueReads_length,curr_kp,setParameters.basicSeqNameDir,setParameters.basicSeqKey,parallelType,[basicSaveName,'_final_',num2str(i)],userDir,[runFileName,'_final_',num2str(i)],setParameters);
  found{i} = zeros(1,setParameters.numBACtoConsider);
  found{i}(curr_kp) = currX;
  found_sumRelevantReads(i) = currSumRelevantReads;
end

disp('standard way of saving')
save(setParameters.saveName,'found','found_sumRelevantReads','store_kp','store_X','store_sumRelevantReads','correctWeight','ind_bac_in_mix')





