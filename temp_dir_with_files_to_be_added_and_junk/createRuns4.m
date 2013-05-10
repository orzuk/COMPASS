%%%%%%%%%%%
% create
clear


cr('s4',[0],[10],[0.5],100,[10^4 10^5 10^6 Inf])



w = [];
k = 1;
for i=[1:10:391]
  w = [w,'./listOfR_s4_',num2str(i),';sleep 30;'];
  if mod(k,5)==0
    w = [w,';sleep 1h; '];
  end
  k = k+1;
end



%rerun s1
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
basicPrefix = 's4';
basicDir = '/CS/BAC/dataForSim/set/sec/';
basicRerunName = [basicPrefix,'_rerun1'];
load([userDir,basicDir,basicPrefix,'_list'])
rerun(totalListOfNames,totalListOfRuns,basicDir,basicPrefix,basicRerunName,userDir);


% del *.e files and run
unix(['cd ',userDir,basicDir,';rm ',basicPrefix,'*.e'])
unix(['cd ',userDir,basicDir,';chmod 0700 ',userDir,basicDir,basicRerunName])


% check the reads are the same
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
basicSeqNameDir = [userDir,'/CS/BAC/dat450000/'];
basicSeqKey =  [userDir,'/CS/BAC/dat450000/key450000'];

cr('s4',[0],[10],[0.5],100,[10^4 10^5 10^6 Inf])

loadRes('s4',0,10,0.5,100,[Inf])
loadRes('s2',1,10,0,50,[10^5])
