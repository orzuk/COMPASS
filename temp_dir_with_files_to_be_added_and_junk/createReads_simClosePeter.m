function createReads_simClosePeter(userDir,mixtureFile,mixtureName,Nreads)
%keyboard
load(mixtureFile)





forBlastFlag = 0;
%userDir = getuserdir;
basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];

addNoiseFlag = 1;
readLength = 100;

ErrorStruct = []; ErrorStruct.error_model = 'exponential';
ErrorStruct.baseline_error = 0.005; 
ErrorStruct.final_error = 0.03; 
                    
p_one_nuc_error_mat = [ -1   0.3  0.22  0.18
                    0.5    -1  0.22   0.6
                    0.35  0.15    -1  0.22
                    0.15  0.55  0.56   -1]; % , 16, length(p));
ErrorStruct.p_one_nuc_error_mat = p_one_nuc_error_mat;

dataName = [mixtureName,'_reads_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];
saveFile = [userDir,'/CS/BAC/toyExample/',mixtureName,'/reads/',dataName];

runSpecificToyExample(correctWeight,Nreads,basicSeqKey,basicSeqNameDir,addNoiseFlag,readLength,saveFile,forBlastFlag,ErrorStruct);
