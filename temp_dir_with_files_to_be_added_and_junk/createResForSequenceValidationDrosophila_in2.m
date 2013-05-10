function createResForSequenceValidationDrosophila_in2(primerName)

global Header_750 Sequence_750 header_uni1to350_primers450  seq_uni1to350_primers450 Header_uni Sequence_uni

clear Sequence_amp* Header_amp*

load ~/CS/BAC/12Samples/validation/scf_files_drosophila/insilocopcr_sanger_drosophila

currSeq = ['Sequence_amp_',primerName];
currSeq = eval(currSeq);

currHeader = ['Header_amp_',primerName];
currHeader = eval(currHeader);


numbers = struct;
numbers.Full = [];
numbers.s750 = [];
numbers.s450 = [];
keyboard

clear num;k=1;
for i=1:length(currSeq)
  for j=1:length(currHeader{i})
    num(k,:) = [findNumberFromHeader(currHeader{i}{j}),i,j];
    k = k+1;
  end
end


for i=1:length(currSeq)
  for j=1:length(currHeader{i})
    [i,j]
    [s750.Header,s750.Sequence,s750.Numbers,s450.Header,s450.Sequence,s450.Numbers,Full.Header,Full.Sequence,Full.Numbers]=findPositionOfAmpSequence(currSeq{i},currHeader{i}{j});
    
    numbers.Full = [numbers.Full;[Full.Numbers i*ones(length(Full.Numbers),1) j*ones(length(Full.Numbers),1)]];
    
    if ~isempty(s750.Numbers)
      numbers.s750 = [numbers.s750;[s750.Numbers i*ones(length(s750.Numbers),1) j*ones(length(s750.Numbers),1)]];
    end
    
    if ~isempty(s450.Numbers)
      numbers.s450 = [numbers.s450;[s450.Numbers i*ones(length(s450.Numbers),1) j*ones(length(s450.Numbers),1)]];
    end   
  
  end
end

allNumbers.Numbers = unique([numbers.Full(:,1);numbers.s750(:,1);numbers.s450(:,1)]);

% for each of these  allNumbers find the sequence if appears in 750,450 and the amplified which contains it

cData = struct;
cData.Numbers = allNumbers.Numbers;
for i=1:length(allNumbers.Numbers)
  i
  clear ampIndex_Full ampIndex_s750 ampIndex_s450
  % find the Full sequence
  clear found
  flag = 0;
  j = 1;
  while flag==0 && j<length(Sequence_uni)
    for k=1
      name = Header_uni{j};
      a = find(name==' ');
      name = name(1:a(1)-1);
      name = str2num(name);
      if name==allNumbers.Numbers(i)
        flag = 1;
        found = [j];
      end
    end
    j = j+1;
  end
  if flag==0
    disp('problem');pause
  else
    length(found)
    cData.Full{i}.Sequence = Sequence_uni{found(1)};
    cData.Full{i}.Header = Header_uni{found(1)};
  end
  %%%%%%%%%%%%%%%%%%%%%%%% 
  
  ampIndex_Full = find(numbers.Full(:,1)==allNumbers.Numbers(i));
  ampIndex_s750 = find(numbers.s750(:,1)==allNumbers.Numbers(i));
  ampIndex_s450 = find(numbers.s450(:,1)==allNumbers.Numbers(i));
  if length(ampIndex_Full)>1
    disp('problem');pause
  end
  
  % write the amplified region that is connected to the current seuqence
  if isempty(ampIndex_Full)
    relevantAmpSeq = unique([numbers.s750(ampIndex_s750,2)';numbers.s450(ampIndex_s450,2)])
  else
    relevantAmpSeq = numbers.Full(ampIndex_Full,2);
    if length(relevantAmpSeq)>1
      disp('problem');pause
    end
  end
  for k=1:length(relevantAmpSeq)
      cData.amp{i}.Sequence{k} = currSeq{relevantAmpSeq(k)};
      cData.amp{i}.Header{k} = currHeader{relevantAmpSeq(k)};
  end
  %%%%%%%%%%%%% amplified
  
  % 750 find the sequence of and the 
  if ~isempty(ampIndex_s750) % find the 750 relevant sequence
    clear found
    flag = 0;
    j = 1;
    %keyboard
    while flag==0 && j<length(Sequence_750)
      for k=1:length(Header_750{j})
        name = Header_750{j}{k};
        a = find(name==' ');
        name = name(1:a(1)-1);
        name = str2num(name);
        if name==allNumbers.Numbers(i)
          flag = 1;
          found = [j,k];
        end
      end
      j = j+1;
    end
    
    if flag==0
      disp('problem');pause
    else
      cData.s750{i}.Sequence = Sequence_750{found(1)};
      cData.s750{i}.Header = Header_750{found(1)};
    end
  else
    cData.s750{i}.Sequence = [];
    cData.s750{i}.Header = [];
  end
  %%%%%%%%%%%%%%%%%% 750
  
  % 450 find the sequence of and the 
  if ~isempty(ampIndex_s450) % find the 450 relevant sequence
    clear found
    flag = 0;
    j = 1;
    while flag==0 && j<length(seq_uni1to350_primers450)
      for k=1:length(header_uni1to350_primers450{j})
        name = header_uni1to350_primers450{j}{k};
        a = find(name==' ');
        name = name(1:a(1)-1);
        name = str2num(name);
        if name==allNumbers.Numbers(i)
          flag = 1;
          found = [j,k];
        end
      end
      j = j+1;
    end
    if flag==0
      disp('problem');pause
    else
      cData.s450{i}.Sequence = seq_uni1to350_primers450{found(1)};
      cData.s450{i}.Header = header_uni1to350_primers450{found(1)};
    end
  else
    cData.s450{i}.Sequence = [];
    cData.s450{i}.Header = [];
  end
  %%%%%%%%%%%%%%%%% 450
  
  
end
% amp sequence is the only one which can have more than one sequence
%  cData.s450{i}.Header  - may have multiple entries. One of them matches the full


% align amp and 750 and 450 to Full
for i=1:length(allNumbers.Numbers)
  
  seq2 = cData.Full{i}.Sequence;
  % align amp
  for j=1:length(cData.amp{i}.Sequence)
    seq1 = cData.amp{i}.Sequence{j};
    AlignStruct = localalign(seq1, seq2,'Alphabet','NT');
    cData.amp{i}.Start(j,:) = AlignStruct.Start;
    cData.amp{i}.Stop(j,:) = AlignStruct.Stop;
    cData.amp{i}.Readme = 'the first is the amplified sequence';
  end
  
  if ~isempty(cData.s750{i}.Sequence)
    seq1 = cData.s750{i}.Sequence;
    AlignStruct = localalign(seq1, seq2,'Alphabet','NT');
    cData.s750{i}.Start = AlignStruct.Start;
    cData.s750{i}.Stop = AlignStruct.Stop;
    cData.s750{i}.Readme = 'the first is the 750 sequence';     
  end
  
  if ~isempty(cData.s450{i}.Sequence)
    seq1 = cData.s450{i}.Sequence;
    AlignStruct = localalign(seq1, seq2,'Alphabet','NT');
    cData.s450{i}.Start = AlignStruct.Start;
    cData.s450{i}.Stop = AlignStruct.Stop;
    cData.s450{i}.Readme = 'the first is the 450 sequence';     
  end
  
end

save(['~/CS/BAC/12Samples/validation/scf_files_drosophila/resModified_',primerName],'cData') 



