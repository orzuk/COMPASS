function [readsFileName,red]=newReadsBasedSeed2Or(reads_data_fileName,userDir,basicSeqNameDir,basicSeqKey,N,addNoiseFlag,ErrorStruct)
%keyboard
[Nbac_in_mixture,readlen,npower,bac_dist_flag,Nread]=findParametersReads(reads_data_fileName);

load(reads_data_fileName)

load([userDir,'/CS/BAC/cellsimbac2']) 
%keyboard
outfilename = reads_data_fileName;
%[readsFileName,red] = fun_prep_data_for_simulations_loadSpecificBAC2(basicSeqKey,basicSeqNameDir,indsimbac_vec,Nbac_in_mixture,Nread,readlen,npower,bac_dist_flag,outseed,outfilename,N,1);
%keyboard

[readsFileName,red] = fun_prep_data_for_simulations_loadSpecificBAC2OrFourth(basicSeqKey,basicSeqNameDir,indsimbac_vec,Nbac_in_mixture,Nread,readlen,npower,bac_dist_flag,outseed,outfilename,N,1,0,addNoiseFlag,ErrorStruct);
