function oneNodeForL1OrFourth(fileNameL1data)

% rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));addpath('/homes/csfaculty/shental/CS/mFilesBAC/');cd /homes/csfaculty/shental/CS/BAC/cvx; mcc -m oneNodeForL1OrFourth.m -a addL1toL2_res.m -a prepareGroupOf1000DistributedSequenceFilesOrFourth.m -a currReadsFourth.m -a myCorrectReadsNew3_for_SL4  -a testL1_2.m -a ./builtins -a ./commands -a ./functions -a ./keywords -a ./lib -a ./sdpt3 -a ./sedumi -a ./sets -a ./structures

load(fileNameL1data)

addL1toL2_res(uniqueReads,uniqueReads_length,tmpInd,fileNameL1Res,auxData)
