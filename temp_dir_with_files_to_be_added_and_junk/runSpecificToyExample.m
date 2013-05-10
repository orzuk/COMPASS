function runSpecificToyExample(correctWeight,Nreads,basicSeqKey,basicSeqNameDir,addNoiseFlag,readLength,saveFile,forBlastFlag,ErrorStruct)

auxData = struct;
auxData.correctWeight = correctWeight;
auxData.Nreads = Nreads;
auxData.basicSeqKey = basicSeqKey;
auxData.basicSeqNameDir = basicSeqNameDir;
auxData.addNoiseFlag = addNoiseFlag;
auxData.readLength = readLength;
if auxData.addNoiseFlag==1
  auxData.ErrorStruct = ErrorStruct;
end
%keyboard
[red]=createReadsForSpecificMixture(auxData);


[uniqueReads,uniqueReads_inds] = extract_sub_kmers(red, auxData.readLength*ones(size(red,1),1),auxData.readLength, 1,0);
%keyboard
if isempty(find(uniqueReads(1,:)-zeros(1,size(uniqueReads,2),'uint64')))
  disp('no reads for some bacteria. try rerun');
  return
end

clear red

[junk_vals junk_inds uniqueReads_length] = get_duplicates(uniqueReads_inds(:,1));
clear junk_vals junk_inds uniqueReads_inds

if forBlastFlag
  %keyboard
  reads_uni = int2nt(unpack_seqs(uniqueReads,auxData.readLength,64));
  save(saveFile,'uniqueReads','uniqueReads_length','reads_uni');
else
  save(saveFile,'uniqueReads','uniqueReads_length');
end



