
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%
%%%%%%
%%%%%
% 25AprilAfternoon



clear
basicName = {'reads-1000-noise-0',...
'reads-1000-noise-100',...
'reads-1000-noise-10',...
'reads-1000-noise-1',...
'reads-1000-noise-S1',...
'reads-100-noise-0',...
'reads-100-noise-100',...
'reads-100-noise-10',...
'reads-100-noise-1',...
'reads-100-noise-S1',...
'reads-10-noise-0',...
'reads-10-noise-100',...
'reads-10-noise-10',...
'reads-10-noise-1',...
'reads-10-noise-S1'};
auxDataTest = struct;

auxDataTest.datasetName = 'primers750_primer_tail';
auxDataTest.numBACtoConsider = 212040;
auxDataTest.dataSetPrimers750Flag = 1;

auxDataTest.readLengthPre = 100;
auxDataTest.readLengthPost = 50;
auxDataTest.readLength = 50;

auxDataTest.currNumProcessors = 20;
auxDataTest.basicDirNameForFile = 'sim25AprilAfternoon_50';

clear fileName inputName outputName
for i=1:length(basicName)
  fileName{i} = ['myFormat_',basicName{i}];
  inputName{i} = basicName{i};
  outputName{i} = ['test_',auxDataTest.basicDirNameForFile,'_',basicName{i}];
end

for i=1:length(basicName)
  createUniqueReads_divideReads(['~/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,'/',inputName{i}],['~/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,'/',fileName{i}],auxDataTest.readLengthPre,auxDataTest.readLengthPost)
end


% copy to br
PWD = pwd;
cd(['~/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile])
clear w
w = 'scp '
for i=1:length(basicName)
  w = [w,' ',fileName{i},'.mat '];
end
w = [w,' orzuk@tin.broadinstitute.org:/seq/orzuk2/compressed_sensing/metagenomics/next_gen','/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,'/'];

cd(PWD);


% prepare the files

for i=1:length(outputName)
  testSimFromReadsLikeRealData_generic(outputName{i},['simFollowingSolexa/',auxDataTest.basicDirNameForFile,'/',fileName{i}],auxDataTest)
end

%%%%%5
% run the files

for i=1:length(outputName)
  unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName{i},'_',outputName{i},'/run_sol_',outputName{i},'_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength)])
end

% collect results
for i=1:length(outputName)
  unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName{i},'_',outputName{i},'/sol_',outputName{i},'_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength),'/sol_',outputName{i},'_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength),'.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,'/results/'])
end

% copy to sol
%unix(['cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,';tar cf results_', ...
%      auxDataTest.basicDirNameForFile,'.tar ','/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,'/results/'])

unix(['cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,';tar cf results_', ...
      auxDataTest.basicDirNameForFile,'.tar results/'])

% copy locally
unix(['scp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,'/results_', ...
      auxDataTest.basicDirNameForFile,'.tar ',...
      'shental@sol.cslab.openu.ac.il:~/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile])
% run locally

for i=1:length(outputName)
  extractResSim_generic(outputName{i},inputName{i},fileName{i},auxDataTest)
  drawnow
  %pause
end

for i=1:length(outputName)
  plotResSim_generic(outputName{i},inputName{i},fileName{i},auxDataTest)
  drawnow
  pause
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% divide S1


