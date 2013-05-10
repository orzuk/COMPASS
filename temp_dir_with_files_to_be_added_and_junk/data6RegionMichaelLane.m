
% check rank in each region

clear
%%%%%%%%%%%55
bcs_data_path = '~/CS/BAC/'

clear fracRelevantReads sumRelevantReads reads F tmpMat
reads = cell(6,1);
fracRelevantReads = reads;
dataIn = struct;
for i=1:6
  tmp = load([bcs_data_path,'Niv/SR50/barcode3_region',num2str(i),'_unireads_freq_norm.mat']);
  reads{i}.uniqueReads = pack_seqs(tmp.reads_uni,64);
  reads{i}.uniqueReads_length = tmp.reads_uni_freq;
  
  [fracRelevantReads{i},sumRelevantReads{i}] = currReadsFourth(reads{i}.uniqueReads,reads{i}.uniqueReads_length,values{i},0,dataIn);

  F{i} = fracRelevantReads{i}./sum(fracRelevantReads{i});
  
  tmpMat{i} = sparse(size(fracRelevantReads{i},1),6);
  tmpMat{i}(:,i) = -F{i};
end

sumRelevantReads

Ap=[M{1},tmpMat{1};M{2},tmpMat{2};M{3},tmpMat{3};M{4},tmpMat{4};M{5},tmpMat{5};M{6},tmpMat{6}];
Ap(:,1) = []; % do we need this?
[x] = solve_L2_3(Ap,zeros(size(Ap,1),1));


res = x(1:end-6);
res = res./sum(res);

min(res)
max(res)

a = find(res>10^-3)

[junk,i] = sort(res(a),'descend');

a = a(i);

found_bact_in_26 = a;

seq_database_file = fullfile(bcs_data_path, '/full16S/bac16s_full_without_ambiguous'); % Database with all sequences
S = load(seq_database_file); 
S.Sequence_uni = vec2column(S.Sequence_uni); 


fid = fopen([bcs_data_path,'Niv/SR50/sample3_res_Jan2013.txt'],'w');
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
save([bcs_data_path,'Niv/SR50/sample3_res_indices_moreThan1e3_Jan2013'],'found_bact','found_bact_in_26','res')

%%%%%%%%%%

% check reverse reads:
if 1==2
reg = 4;
tmp = load([bcs_data_path,'Niv/SR50/barcode3_region',num2str(reg),'_unireads_freq_norm.mat']);
[Header, Sequence_r] = fastaread(['~/CS/BAC/Niv/SR50/barcode3_region',num2str(reg),'_reverse_primer.fasta']);
[Header, Sequence_f] = fastaread(['~/CS/BAC/Niv/SR50/barcode3_region',num2str(reg),'_forward_primer.fasta']);

seq_for = char(ones(length(Sequence_f),26));
for i=1:length(Sequence_f)
  seq_for(i,:) = Sequence_f{i}(26:end);
end
[u_f] =  unique(seq_for,'rows');

seq_rev = char(ones(length(Sequence_r),26));
for i=1:length(Sequence_r)
  seq_rev(i,:) = seqrcomplement(Sequence_r{i}(26:end));
end
[u_r] =  unique(seq_rev,'rows');


AmitReads = load([bcs_data_path,'Niv/SR50/barcode3_region',num2str(reg),'_unireads_freq_norm.mat']);


[junkAmit,i1Amit]=intersect(int2nt(unpack_seqs(values{reg},readLength,64)),[u_f;u_r],'rows')

[junkNoam,i1Noam]=intersect([seq_f(:,26*(reg-1)+1:26*reg);seq_r(:,26*(reg-1)+1:26*reg)],[u_f;u_r],'rows');% mapped 

size(junkNoam)
size(junkAmit)

end


