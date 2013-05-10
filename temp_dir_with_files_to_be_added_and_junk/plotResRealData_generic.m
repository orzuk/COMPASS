function plotResRealData_generic(resultsDir,outputName,auxDataTest)
%keyboard


load([resultsDir,'/',auxDataTest.basicDirNameForFile,'/res_',outputName])
%keyboard
numBACtoConsider = auxDataTest.numBACtoConsider;

figure(4)
plot(f_l1,f_l2,'.')
xlabel('l1');ylabel('l2')

name = outputName;
name(find(name=='_')) = '-';
title(name)
