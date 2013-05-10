clear
load ~/CS/BAC/primers_startAs750_length350/bac16s_1to350_primers450_HeaderForRunList
load ~/CS/BAC/primers750_primer_tail/bac16s_primers750_primer_tail_full_without_ambiguous
load('~/CS/BAC/bacteria_s16_data_uni','Header_uni') 


% A
%%%%%%%%%%
nm = [274161 58367 109865 136169 136306 159798 237436 333353 343962 345192 351827 534066 535023 551356 556801 557726 571171 572809 4418 224755 211737 ...
 159924 208693 572594 593417]

name = addNamesForMichael_in(nm,header_uni1to350_primers450);
for i=1:length(name)
  fprintf('%s\n',name{i});
end

nm = [58371];
name = addNamesForMichael_in(nm,header_uni1to350_primers450);
for i=1:length(name)
  fprintf('%s\n',name{i});
end

nm = 65495
name = addNamesForMichael_in(nm,header_uni1to350_primers450);
for i=1:length(name)
  fprintf('%s\n',name{i});
end

nm = [4452 4469 272985]
name = addNamesForMichael_in(nm,header_uni1to350_primers450);
for i=1:length(name)
  fprintf('%s\n',name{i});
end

nm = [4421 58366];
name = addNamesForMichael_in(nm,header_uni1to350_primers450);
for i=1:length(name)
  fprintf('%s\n',name{i});
end

nm = [156437 41554 27422 34622 58364 160003 209563 212408 213145 34202];
name = addNamesForMichael_in(nm,header_uni1to350_primers450);
for i=1:length(name)
  fprintf('%s\n',name{i});
end

nm = [92777 25818 38754 39658 137919 249632 278839 535203 229477];
name = addNamesForMichael_in(nm,header_uni1to350_primers450);
for i=1:length(name)
  fprintf('%s\n',name{i});
end

nm = [325754];
name = addNamesForMichael_in(nm,header_uni1to350_primers450);
for i=1:length(name)
  fprintf('%s\n',name{i});
end

%%%%%%%%%%%%%

% W
nm = [108447 173020 273831 18468 41221 52243 161967 518659 30977 52332 34685];
name = addNamesForMichael_in(nm,header_uni1to350_primers450);
for i=1:length(name)
  fprintf('%s\n',name{i});
end

nm = [564629];
name = addNamesForMichael_in(nm,header_uni1to350_primers450);
for i=1:length(name)
  fprintf('%s\n',name{i});
end

nm = [91297 91565 91667 91050];
name = addNamesForMichael_in(nm,header_uni1to350_primers450);
for i=1:length(name)
  fprintf('%s\n',name{i});
end

nm = [91789];
name = addNamesForMichael_in(nm,header_uni1to350_primers450);
for i=1:length(name)
  fprintf('%s\n',name{i});
end

nm = [170818 579340];
name = addNamesForMichael_in(nm,header_uni1to350_primers450);
for i=1:length(name)
  fprintf('%s\n',name{i});
end

nm = [130782];
name = addNamesForMichael_in(nm,header_uni1to350_primers450);
for i=1:length(name)
  fprintf('%s\n',name{i});
end

nm = [244093];
name = addNamesForMichael_in(nm,header_uni1to350_primers450);
for i=1:length(name)
  fprintf('%s\n',name{i});
end

nm = [117492];
name = addNamesForMichael_in(nm,header_uni1to350_primers450);
for i=1:length(name)
  fprintf('%s\n',name{i});
end


name = addNamesForMichael_in(nm,Header_uni);
for i=1:length(name)
  fprintf('%s\n',name{i});
end
