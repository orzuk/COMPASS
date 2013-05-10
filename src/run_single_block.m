%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Simple Reconstruction for one block
userDir = '../../matlab/libs/bcs_nextgen/data/';


num = 100; % Number of bacteria in the mixture
list = randperm(400000);
ind_bac_in_mix = list(1:num);
other = list(num+1:num+900); % indices of bacteria not in mixture

tmpInd = [ind_bac_in_mix,other];

correctWeight = zeros(1,[length(ind_bac_in_mix)+length(other)]);
correctWeight(ind_bac_in_mix) = 1/num*ones(1,num);
correctWeight = correctWeight';
basicSeqNameDir = fullfile(userDir, 'packed64');
basicSeqKey= fullfile(userDir,'keyNoNonACGT.mat');

readLength = 50;

[uniqueReads,uniqueReads_length,auxData.fracRelevantReadsForInfinity] ... % simulate reads (in the inifinte limit)
    = createReadsForInfiniteNumberOrFourth(ind_bac_in_mix,correctWeight,readLength,basicSeqNameDir,basicSeqKey);

[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOr(readLength,tmpInd,basicSeqNameDir,basicSeqKey); % build mixture matrix

dataIn = struct;
dataIn.fracRelevantReadsForInfinity = auxData.fracRelevantReadsForInfinity;
auxData.inifiniteNumberOfReadsFlag

[fracRelevantReads,sumRelevantReads(i)] = currReads(uniqueReads,uniqueReads_length,values,1,dataIn);

numVariables = size(normalizedBac,2);
cvx_begin
cvx_quiet(true)
variable x(numVariables)
minimize( norm(normalizedBac*x-fracRelevantReads) );
subject to
x >= 0;
cvx_end

x = x./sum(x);

cvx_begin
cvx_quiet(true)
variable xx(numVariables)
minimize( norm(normalizedBac*xx-fracRelevantReads) );
cvx_end

xx = xx./sum(xx);


% L2 over infinite reads is correct. no need for x>0
% larger samples - how do invert the matrix? how do you
how do we store the large matrix
what happens for finite read

