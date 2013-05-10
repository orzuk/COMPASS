clear
load ~/CS/BAC/data750with4mm/database_16s_750primers_with4mm_full_without_ambiguous


a = {'364189 FM873796.1 occupant sou;393238 GQ061950.1 Topographica;237311 AM419960.1 Noma subging;50919',...
    '516705 FJ557621.1 Survey Approach Assess Endotracheal Tube Biofilms and their Origins biofilm extub',...
    '243980 AM420128.1 Noma subgingival plaque clone 401G10(oral);419923 GQ073133.1 Topographical and Te',...
    '502815 GQ006849.1 Topographical and Temporal Human Skin Microbiome skin occiput clone nbu88e07c1;45',...
    '223790 EF509061.1 Loss During Antibiotic Treatment 1 Intubated Patients Colonized Pseudomonas aerug',...
    '529136 FJ557589.1 Survey Appro;29968 AF432133.1 tongue dorsum;417057 GQ026688.1 Topographica;468282',...
    '461063 GQ058005.1 Topographical and Temporal Human;112388 AY880055.1 Prevotella sp. T05-04;571842 E',...
    '35402 AF385551.1 crevicular epithelial cells clone BE073;111426 AY880053.1 Prevotella sp. N19-22;11',...
    '519814 FJ557642.1 Survey Approach Assess Endotracheal Tube Biofilms and their Origins biofilm extub',...
    '243840 AM420123.1 Noma subgingival plaque clone 401F02(oral);130182 AY959020.1 Microbes on human va',...
    '92496 AY349397.1 human oral cavity clone GI059;576314 FJ976391.1 Interindividual healthy human oral',...
    '427242 GQ032485.1 Topographica;223881 EF510665.1 Loss During ;219179 EF510749.1 Loss During ;217497',...
    '524930 FJ558006.1 Survey Approach Assess Endotracheal Tube Biofilms and their Origins biofilm extub',...
    '380317 GQ025681.1 Topographical and Temporal Human Skin Microbiome skin antecubital fossa clone nbw',...
    '104043 AB108826.1 Prevotella salivae str. JCM 12084 EPSA11;72590 AF385512.1 tongue dorsa clone DO03',...
    '395185 GQ059406.1 Topographical and Temporal Human Skin Microbiome skin volar forearm clone nbw1227',...
    '222653 EF512000.1 Loss During Antibiotic Treatment;223076 EF511882.1 Loss During Antibiotic Treatme',...
    '420192 GQ113006.1 Topographical and Temporal Human Skin Microbiome skin antecubital fossa clone nbw',...
    '28952 AF385513.1 tongue dorsa clone DO039',...
    '220502 EF510744.1 Loss During Antibiotic Treatment 1 Intubated Patients Colonized Pseudomonas aerug',...
    '594033 EU663608.1 throughput sequencing expands microbiology spectrum brain abscess clone 9;83456 A',...
    '539347 GU470887.1 Prevotella sp. oral taxon 302 str. F0323',...
    '271710 EU137441.1 Bartonella-positive fleas: and assembly patterns Oropsylla hirsuta (prairie dog f',...
    '221666 EF511866.1 Loss During Antibiotic Treatment 1 Intubated Patients Colonized Pseudomonas aerug',...
    '257880 AM420091.1 Noma subgingival plaque clone 302E06(oral);500537 GQ009128.1 Topographical and Te',...
    '46885 AY005072.1 subgingival dental plaque clone AU126'};



F = zeros(length(a),1);
for i=1:length(a)
  i
  flag = 0;
  k = 1;
  while flag==0 && k<=length(Header_no_nonACGT)
    for j=1:length(Header_no_nonACGT{k})
      if ~isempty(findstr(Header_no_nonACGT{k}{j},a{i}))
        F(i) = [k];
        flag = 1;
      end
    end
    k = k+1;
  end
end

F = F(find(F));



S7 = load('~/CS/BAC/12Samples/Solexa/data/window_90_withNor_MergeForwardAndReverse/data_window_90_withNor_MergeForwardAndReverse_normReadsf_S7');


basicSeqNameDir = '~/CS/BAC/data750with4mm/datNoNonACGT/packed64/'
basicSeqKey = '~/CS/BAC/data750with4mm/datNoNonACGT/keyNoNonACGT_data750with4mm';
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(90,F,basicSeqNameDir,basicSeqKey);

dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(S7.uniqueReads,S7.uniqueReads_length,values,0,dataIn);



load ~/CS/BAC/12Samples/Solexa/results/window_90_withNor_MergeForwardAndReverse_data750withMM/resWeightedMatrixBased_test_window_90_withNor_MergeForwardAndReverse_data750withMM_normReadsf_S7

load ~/CS/BAC/12Samples/Solexa/results/window_90_withNor_MergeForwardAndReverse_data750withMM/sol_test_window_90_withNor_MergeForwardAndReverse_data750withMM_normReadsf_S7_noCorrection_90



for i=1:length(F)
  k = cellfun(@(x) findstr(x,Sequence_no_nonACGT{F(i)}),Sequence_no_nonACGT,'uniformoutput',0);
end

S = cell(1,length(F));
for i=1:length(S)
  i
  k = 1;
  while k<=length(Sequence_no_nonACGT)
    if ~isempty(findstr(Sequence_no_nonACGT{k},Sequence_no_nonACGT{F(i)}))
        S{i} = [S{i},k];
    end
    k = k+1;
  end
end

len = zeros(length(Sequence_no_nonACGT),1);
for i=1:length(Sequence_no_nonACGT)
  len(i) = length(Sequence_no_nonACGT{i});
end

u = unique(len);

for i=1:length(u)
  i
  c = find(len==u(i));
  [cc,i1,i2] = unique(cell2mat(Sequence_no_nonACGT(c)'),'rows');
  if length(c)~=size(cc,1)
    disp('found')
    pause
  end
end




same = sparse(length(Sequence_no_nonACGT),length(Sequence_no_nonACGT));
for i=1:length(Sequence_no_nonACGT)-1
  if mod(i,1000)==0
    i
  end
  for j=i+1:length(Sequence_no_nonACGT)
    if length(Sequence_no_nonACGT{i})==length(Sequence_no_nonACGT{j})
      if isempty(find(Sequence_no_nonACGT{i}-Sequence_no_nonACGT{j}))
        same(i,j) = 1;
        disp([i,j])
      end
    end
  end
end
  

%%%%%%%%%%%%
clear
load /home/csfaculty/shental/CS/BAC/primers750_primer_tail/bac16s_primers750_primer_tail_full_without_ambiguous

len = zeros(length(Sequence_750_tail),1);
for i=1:length(Sequence_750_tail)
  len(i) = length(Sequence_750_tail{i});
end

u = unique(len);

for i=1:length(u)
  i
  c = find(len==u(i));
  [cc,i1,i2] = unique(abs(cell2mat(Sequence_750_tail(c)')),'rows');
  if length(c)~=size(cc,1)
    disp('found')
    pause
  end
end

%%%%%%%



load ~/CS/BAC/12Samples/Solexa/results/window_90_withNor_MergeForwardAndReverse_data750withMM/res_test_window_90_withNor_MergeForwardAndReverse_data750withMM_normReadsf_S7