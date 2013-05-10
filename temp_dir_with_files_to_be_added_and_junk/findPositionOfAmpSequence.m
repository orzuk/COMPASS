function [curr_Header750,curr_Sequence750,curr_Numbers750,curr_Header450,curr_Sequence450,curr_Numbers450,curr_HeaderFull,curr_SequenceFull,curr_NumberFull]=findPositionOfAmpSequence(currSeq,currHeader)

global Header_750 Sequence_750 header_uni1to350_primers450  seq_uni1to350_primers450 Header_uni Sequence_uni


currNumber = currHeader;
a = find(currNumber==' ');
currNumber = currNumber(1:a(1)-1);
currNumber = str2num(currNumber);


%

% find header in 750
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
    if name==currNumber
      flag = 1;
      found = [j,k];
    end
  end
  j = j+1;
end

if flag==0
  curr_Header750 = [];
  curr_Sequence750 = [];
  curr_Numbers750 = [];
else
  curr_Sequence750 = Sequence_750{found(1)};
  curr_Header750 = Header_750{found(1)}{found(2)};
  curr_Numbers750 = zeros(length(Header_750{found(1)}),1);
  for k=1:length(Header_750{found(1)})
    curr_Numbers750(k) = findNumber(Header_750{found(1)}{k});
  end
  % add all numbers that share the same 750  
end
%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

% find header in 450
clear found
flag = 0;
j = 1;
while flag==0 && j<length(seq_uni1to350_primers450)
  for k=1:length(header_uni1to350_primers450{j})
    name = header_uni1to350_primers450{j}{k};
    a = find(name==' ');
    name = name(1:a(1)-1);
    name = str2num(name);
    if name==currNumber
      flag = 1;
      found = [j,k];
    end
  end
  j = j+1;
end

if flag==0
  curr_Header450 = [];
  curr_Sequence450 = [];
  curr_Numbers450 = [];
else
  curr_Sequence450 = seq_uni1to350_primers450{found(1)};
  curr_Header450 = header_uni1to350_primers450{found(1)}{found(2)};
  curr_Numbers450 = zeros(length(header_uni1to350_primers450{found(1)}),1);
  for k=1:length(header_uni1to350_primers450{found(1)})
    curr_Numbers450(k) = findNumber(header_uni1to350_primers450{found(1)}{k});
  end
end


%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%



% find in 1500
clear found
flag = 0;
j = 1;
while flag==0 && j<length(Sequence_uni)
  for k=1
    name = Header_uni{j};
    a = find(name==' ');
    name = name(1:a(1)-1);
    name = str2num(name);
    if name==currNumber
      flag = 1;
      found = [j,k];
    end
  end
  j = j+1;
end

if flag==0
  curr_HeaderFull = [];
  curr_SequenceFull = [];
  curr_NumberFull = [];
else
  curr_SequenceFull = Sequence_uni{found(1)};
  curr_HeaderFull = Header_uni{found(1)};
  curr_NumberFull = findNumber(Header_uni{found(1)});
end

%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

function currNumber=findNumber(Header)

currNumber = Header;
a = find(currNumber==' ');
currNumber = currNumber(1:a(1)-1);
currNumber = str2num(currNumber);
