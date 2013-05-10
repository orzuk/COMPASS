function [uniqueReads,uniqueReads_length,auxData]=createReadsData(auxData,userDir)
%disp('does not correctMeasurementForReadErrors createReadsData.m')

if auxData.createReadsFlag 
  if auxData.createSpecificMixtureFlag==1 % for specific mixtures
    disp('running specific mixture. createReadsData.m')
    [red] = createReadsForSpecificMixture(auxData);
  else
    [junk,red] = newReadsBasedSeed2Or(auxData.reads_data_fileName,userDir,auxData.basicSeqNameDir,auxData.basicSeqKey,auxData.numBACtoConsider,auxData.addNoiseFlag,auxData.ErrorStruct);
  end
  
else % load the reads themselves
  load(auxData.reads_data_fileName);
end


[uniqueReads,uniqueReads_inds] = extract_sub_kmers(red, auxData.readLength*ones(size(red,1),1),auxData.readLength, 1,0);

%keyboard

% change 19.3.12
%if auxData.correctMeasurementForReadErrorsFlag
%  auxData.readsForCorrectReads = red;
%end
% end change 19.3.12

clear red

[junk_vals junk_inds uniqueReads_length] = get_duplicates(uniqueReads_inds(:,1));
clear junk_vals junk_inds uniqueReads_inds
        %keyboard


