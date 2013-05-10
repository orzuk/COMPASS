function jobList=runProblematic3(listOfProblematic,runFileName,tmpRunFileName,auxData)

delThese = [];
runThese = [];
for j=1:length(listOfProblematic)
  delThese = [delThese,' ',runFileName,'_partDist_Z',num2str(listOfProblematic(j)),'.o'];
  runThese = [runThese,' ',num2str(listOfProblematic(j))];
end


%unix(['rm -f ',delThese,';sleep 10']);
unix(['rm -f ',delThese]);


tmpName = num2str(runThese);
tmpName(find(tmpName==' ')) = '_';

for jj=1:length(listOfProblematic)
  submitProblematic(runFileName,tmpName,tmpRunFileName,listOfProblematic(jj),auxData);
end

% add them to jobList
jobList = [];
for jj=1:length(listOfProblematic)
  jobList{end+1} = [tmpRunFileName,'_',tmpName,'_',num2str(listOfProblematic(jj))];
end

