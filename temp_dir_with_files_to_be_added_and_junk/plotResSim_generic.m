function plotResSim_generic(outputName,inputName,fileName,auxDataTest)
%keyboard




load(['~/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,'/simRes/res_',inputName])
load(['~/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,'/',inputName],'freq','bact')

numBACtoConsider = auxDataTest.numBACtoConsider;
c = zeros(numBACtoConsider,1);
c(bact) = freq;

subplot(2,2,1)
plot(c,f_l1,'b.')
hold on
plot(c,f_l2,'kx')
legend('l1','l2')
title(inputName)
subplot(2,2,2)
plot(c,f_l1,'kx')
line([0,max(c)],[0 max(c)])
xlabel('correct');ylabel('l1')
subplot(2,2,3)
plot(c,f_l2,'kx')
line([0,max(c)],[0 max(c)])
xlabel('correct');ylabel('l2')

subplot(2,2,4)
plot(f_l1,f_l2,'.')
line([0,max(c)],[0 max(c)])
xlabel('l1');ylabel('l2')
title(fileName)

