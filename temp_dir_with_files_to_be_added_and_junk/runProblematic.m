function jobList=runProblematic(listOfProblematic,runFileName,tmpRunFileName,auxData)

delThese = [];
runThese = [];
for j=1:length(listOfProblematic)
  delThese = [delThese,' ',runFileName,'_partDist_Z',num2str(listOfProblematic(j)),'.o'];
  runThese = [runThese,' ',num2str(listOfProblematic(j))];
end

%runThese(1) = [];
unix(['rm -f ',delThese,';sleep 10']);

% submit them again
tmpName = num2str(runThese);
tmpName(find(tmpName==' ')) = '_';
fid = fopen([runFileName,tmpName,'_shell.sh'],'w');
while fid<0
  pause(5)
  fid = fopen([runFileName,tmpName,'_shell.sh'],'w');
end       
%keyboard
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_${RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'for i in %s; do bsub   -q %s -J %s -o %s_partDist_Z%s.o ~/runInput_donotdelete %s_partDist_%s;sleep 20;done;',runThese,auxData.queueName, [tmpRunFileName,'_',tmpName,'_$i'], ...
        runFileName,'$i',runFileName,'$i');
%fprintf(fid,'bsub   -q  %s -J "%s[%s]" -o "%s_partDist_Z%s.o" ~/runInput_donotdelete %s_partDist_%s;\n',auxData.queueName,tmpRunFileName,num2str(runThese),runFileName,'%I',runFileName,'\$LSB_JOBINDEX');
fprintf(fid,'rm -rf $TMPDIR \n');
fclose(fid)
pause(10)
unix(['chmod 0700 ',runFileName,tmpName,'_shell.sh']);
disp('does it run?')
pause(5)

unix([runFileName,tmpName,'_shell.sh'])
pause(60)

% add them to jobList
jobList = [];
for jj=1:length(listOfProblematic)
  jobList{end+1} = [tmpRunFileName,'_',tmpName,'_',num2str(listOfProblematic(jj))];
end



