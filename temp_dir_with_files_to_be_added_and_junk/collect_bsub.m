% for the master node:

% ['bsub -c 1440 -q ',master_queueName,' -M 33554432',...
%                         ' -o ',userDirForFile,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName,'.o ',...
%                         ' -e ',userDirForFile,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName,'.e ',...
%                         '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_generalOrSecond.sh  /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b ',auxDataFileName]


% in the iterations (which works on a set of 400 of 1000s)

% First try:

%   notice that it's done in batch - maybe this is the problem (iterateParallelDistributedSequenceFilesOr2.m)
% fprintf(fid,'bsub -W 1:30 -M 1000000  -q  %s -W 1:59 -J "%s[1-%d]" -o "%s_partDist_Z%s.o" ~/runInput_donotdelete %s_partDist_%s;\n',auxData.queueName,tmpRunFileName,length(partDist)-1,runFileName,'%I',runFileName,'\$LSB_JOBINDEX');

% Second try: look for those that failed and run thme again (submitProblematic.m)

% fprintf(fid,'for i in %s; do bsub -W 1:30 -M 1000000   -q %s -W 1:59 -J %s -o %s_partDist_Z%s.o ~/runInput_donotdelete %s_partDist_%s;done;',num2str(currPartToRun),auxData.queueName, [tmpRunFileName,'_',tmpName,'_$i'], ...
%           runFileName,'$i',runFileName,'$i');

  
