function rerun(totalListOfNames,totalListOfRuns,basicDir,basicPrefix,basicRerunName,userDir)
%keyboard
unix(['cd ',userDir,basicDir,'; grep -l "Error:Initialize component instance failed" ',basicPrefix,'*.e | cut -f1 -d"." > ',basicPrefix,'_listOfErrors']);


fid = fopen([userDir,basicDir,basicPrefix,'_listOfErrors'],'r');
a = fscanf(fid,'%s\n');
fclose(fid);

fid = fopen([userDir,basicDir,basicRerunName],'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'\n');


b = findstr(a,[basicPrefix,'_']);
b = [b,length(a)+1];
for i=1:length(b)-1
  for j=1:length(totalListOfNames)
    if strcmp(a(b(i):b(i+1)-1),totalListOfNames{j})
      totalListOfRuns{j}
      fprintf(fid,'sleep 30;\n')
      fprintf(fid,'%s\n',totalListOfRuns{j});
    end
  end
end
fclose(fid)


