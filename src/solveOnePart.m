function [x,sumRelevantReads]=solveOnePart(allocateToOneFirst,allocateToOneLast,auxData,basicAllocationFileName,saveFileName,tmpInd,basicSeqNameDir,basicSeqKey,readLength,uniqueReads,uniqueReads_length)


if ~isfield(auxData,'moveDependentToNextStageFlag')
  auxData.moveDependentToNextStageFlag = 0;
end

if ~isfield(auxData,'solveOr')
  auxData.solveOr = 0;
end

if ~isfield(auxData,'do_lp_flag')
  auxData.do_lp_flag = 0;
end


for i=allocateToOneFirst:allocateToOneLast
  
  %keyboard
  [normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd{i},basicSeqNameDir,basicSeqKey);
  %keyboard
  
  dataIn = struct;
    
  disp('taking measurement as is')
  [fracRelevantReads,sumRelevantReads(i)] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);
  if ~isempty(find(fracRelevantReads))
    disp('found reads')
  else
    disp('non found')
  end
  
  if sumRelevantReads(i)>0 
    [x{i}]=runOneGroupOf1000ForCompilationFourth(normalizedBac,fracRelevantReads);    
  else
    x{i} = zeros(size(normalizedBac,2),1);
  end % solve standard way
  
  
  if auxData.moveDependentToNextStageFlag | auxData.solveOr
      %keyboard
      A = full(normalizedBac'*normalizedBac);
      
      
      %%%%%%%%%
      % move dependent to next stage
      if auxData.moveDependentToNextStageFlag
        run_lin_prog = 0;
        [unique_inds x_min x_max] = set_solution_bounds(A, x{i}, fracRelevantReads, 0, run_lin_prog);
        oth = 1:size(A,2);
        oth(unique_inds) = [];
        
        x{i}(oth) = 2; % the are moved to the next stage
        disp(['moved ',num2str(length(oth)),' to the next stage. solveForGroupOrFourth.m'])
        %disp('check this. does it move to the next stage. cn just calc the rank? solveForGroupOrFourth. ')
        %pause 
      end
      %%%%%%%%%%%%%%% end 
      
      if auxData.solveOr 
        %%%%%%%%%%%%%%
        % solve Or
        %keyboard
        if auxData.do_lp_flag
          disp('solve the lp. solveForGroupOrFourth.m')
          run_lin_prog = 1;
          [unique_inds x_min x_max] = set_solution_bounds(A, x{i}, fracRelevantReads, 0, run_lin_prog);
          xx = x{i};
          x{i} = [];
          x{i}(:,1) = xx;
          x{i}(:,2) = x_min;
          x{i}(:,3) = x_max;
        else
          disp('does not solve the lp. just provides the result. Those indices for which x is negative are the unique_inds');
          run_lin_prog = 0;
          [unique_inds x_min x_max] = set_solution_bounds(A, x{i}, fracRelevantReads, 0, run_lin_prog);
          x{i}(unique_inds) = -x{i}(unique_inds);
          
        end
        
     end
      %%%%%%%%%%
  end
end

