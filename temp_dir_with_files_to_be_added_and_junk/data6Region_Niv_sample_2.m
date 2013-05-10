
% check rank in each region

clear
%%%%%%%%%%%55
bcs_data_path = '~/CS/BAC/'
load([bcs_data_path,'full16S/data_primersAmit_0MM_readLength_26'])


% load Niv reads - sample 3

clear fracRelevantReads sumRelevantReads reads F tmpMat
reads = cell(6,1);
fracRelevantReads = reads;
dataIn = struct;
for i=1:6
  tmp = load([bcs_data_path,'Niv/SR50/barcode2_region',num2str(i),'_unireads_freq_norm.mat']);
  reads{i}.uniqueReads = pack_seqs(tmp.reads_uni,64);
  reads{i}.uniqueReads_length = tmp.reads_uni_freq;
  
  [fracRelevantReads{i},sumRelevantReads{i}] = currReadsFourth(reads{i}.uniqueReads,reads{i}.uniqueReads_length,values{i},0,dataIn);

  F{i} = fracRelevantReads{i}./sum(fracRelevantReads{i});
  
  tmpMat{i} = sparse(size(fracRelevantReads{i},1),6);
  tmpMat{i}(:,i) = -F{i};
end

sumRelevantReads

%Ap=[M{1},tmpMat{1};M{2},tmpMat{2};M{3},tmpMat{3};M{4},tmpMat{4};M{5},tmpMat{5};M{6},tmpMat{6}];
%Ap=[M{2},tmpMat{2};M{3},tmpMat{3};M{4},tmpMat{4};M{5},tmpMat{5};M{6},tmpMat{6}];
%Ap=[M{1},tmpMat{1};M{2},tmpMat{2};M{3},tmpMat{3};M{5},tmpMat{5};M{6},tmpMat{6}];

%Ap=[M{3},tmpMat{3};M{4},tmpMat{4};M{5},tmpMat{5};M{6},tmpMat{6}];
%Ap=[M{4},tmpMat{4};M{5},tmpMat{5}];
%Ap=[M{4},tmpMat{4};M{5},tmpMat{5};M{6},tmpMat{6}];
%Ap=[M{3},tmpMat{3};M{5},tmpMat{5};M{6},tmpMat{6}];
%Ap=[M{5},tmpMat{5};M{2},tmpMat{2}];
%Ap=[M{6},tmpMat{6}];

%1
%pause
Ap=[M{1},tmpMat{1};M{2},tmpMat{2};M{3},tmpMat{3};M{4},tmpMat{4};M{5},tmpMat{5};M{6},tmpMat{6}];
Ap(:,1) = []; % do we need this?
[x] = solve_L2_3(Ap,zeros(size(Ap,1),1));


res = x(1:end-6);

min(res)
max(res)

a = find(res>10^-3)

[junk,i] = sort(res(a),'descend');

a = a(i);

found_bact_in_26 = a;


seq_database_file = fullfile(bcs_data_path, '/full16S/bac16s_full_without_ambiguous'); % Database with all sequences
S = load(seq_database_file); 
S.Sequence_uni = vec2column(S.Sequence_uni); 

fid = fopen([bcs_data_path,'Niv/SR50/sample2_res_Jan2013.txt'],'w');
for i=1:length(a)
  fprintf(fid,'frequency: %1.2f\n',junk(i));
  currSeq = find(group==group(inds(a(i))));
  %S.Header_uni{currSeq}
  for j=1:length(currSeq)
    fprintf(fid,'%s\n',S.Header_uni{currSeq(j)});
  end
  fprintf(fid,'%%%%%%%%%%\n');
end
fclose(fid);
found_bact = inds(a);
save([bcs_data_path,'Niv/SR50/sample2_res_indices_moreThan1e3_Jan2013'],'found_bact','found_bact_in_26','res')

%%%%%%%%%%

