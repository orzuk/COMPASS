clear 
addpath(genpath('~/CS/BAC/cvx'))
addpath('~/CS/mFilesBAC')
userDir = getuserdir;
basicSeqNameDir = [userDir,'/CS/BAC/dat450000/'];
basicSeqKey =  [userDir,'/CS/BAC/dat450000/key450000'];



file = '~/CS/BAC/dataForSim/simSec_Nbacmix_100_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat';
     
load ~/CS/BAC/testSimSec
load(file) 

plotRes(file,'~/CS/BAC/testSimSec')

n = length(found{end});
[junk,red] = newReadsBasedSeed2(file,userDir,basicSeqNameDir,basicSeqKey,n);
[uniqueReads,uniqueReads_inds] = myUniqueUINT8(red);
clear red
uniqueReads_length = zeros(size(uniqueReads,1),1);
for i=1:size(uniqueReads_length,1)
  uniqueReads_length(i) = length(uniqueReads_inds{i});
end
clear uniqueReads_inds

n = length(found{end});
readLength = 50;
tmpInd = store_kp{end};

[normalizedBac_red values_red] = prepareGroupOf1000DistributedSequenceFiles2(readLength,tmpInd,basicSeqNameDir,basicSeqKey); 

[fracRelevantReads_red,sumRelevantReads_red] = currReads(uniqueReads,uniqueReads_length,values_red);
[x_red]=runOneGroupOf1000post(1,normalizedBac_red,fracRelevantReads_red,1);

f_red = zeros(1,n);
f_red(tmpInd) = x_red;
plot(f_red,correctWeight,'.')

% compare the costs
sumsqr(normalizedBac_red*correctWeight(tmpInd)-fracRelevantReads_red)
sumsqr(normalizedBac_red*f_red(tmpInd)'-fracRelevantReads_red)


Yc= normalizedBac_red*correctWeight(tmpInd);
plot(Yc,fracRelevantReads_red,'.')


z = find(Yc>6*10^-4 & fracRelevantReads_red<0.2*10^-3);
fracRelevantReads2 = fracRelevantReads_red;
fracRelevantReads2(z) = Yc(z);

[x_red2]=runOneGroupOf1000post(1,normalizedBac_red,fracRelevantReads2,1);

f_red = zeros(1,n);
f_red(tmpInd) = x_red2;
plot(f_red,correctWeight,'.')


spec = char(values_red(z(1),:));
q = tmpInd(find(normalizedBac_red(z(1),:)));


load ~/CS/BAC/bacteria_s16_data_uni

Sequence_uni{q}

cell2mat(strfind(Sequence_uni(q),spec))./len_uni(q)


for i=1000
  a(i,:) = randi([1,10],1,1000);
end


Header_uni{q}



% give lower weight based on the number bac it appears in 
A = normalizedBac_red;
A(find(A)) = 1;
weight = (size(A,2)-sum(A,2))./size(A,2);
z = weight*ones(1,size(normalizedBac_red,2));


a = sum(normalizedBac_red,2);
weight = 1-a./sum(a);
delta = max(weight)-min(weight);
weight = (weight-min(weight))./delta;
%weight = weight.^100;
weight= (weight)/sum(weight);
%weight = sqrt(weight);
%weight = weight.^1.45;
z = weight*ones(1,size(normalizedBac_red,2));

[x_red2] = runOneGroupOf1000post(1,normalizedBac_red.*z,fracRelevantReads_red.*weight,1);

f_red2 = zeros(1,n);
f_red2(tmpInd(b)) = x_red2;
plot(f_red2,correctWeight,'.')

plot(normalizedBac_red.*z*x_red2,fracRelevantReads_red.*weight,'.')

% by what is 51 different from it's close ones.



sumsqr(normalizedBac_red.*z*x_red2-fracRelevantReads_red.*(weight./sum(weight)))

sumsqr(normalizedBac_red.*z*(correctWeight(tmpInd(b))/sum(correctWeight(tmpInd(b))))-fracRelevantReads_red*(weight./sum(weight)))


%%%%%%%%%%%
n = length(found{end});
readLength = 50;
tmpInd = store_kp{end};
x = tmpFun(readLength,tmpInd,basicSeqNameDir,basicSeqKey) 
f = zeros(1,n);
f(tmpInd) = x;

xx{1} = x;
currInd = tmpInd;
for i=2:10
  currInd = currInd(find(xx{i-1}>10^-3));
  
  xx{i}=tmpFun(readLength,currInd,basicSeqNameDir,basicSeqKey,uniqueReads,uniqueReads_length) ;
  f(i,:) = zeros(1,n);
  f(i,currInd) = xx{i};
  plot(f(i,:),correctWeight,'.');
  length(currInd)
  pause
end
