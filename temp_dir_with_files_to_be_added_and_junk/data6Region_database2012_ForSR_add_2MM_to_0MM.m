% add bacteria that share the same sequence as in 0MM but do not appear
clear
bcs_data_path = '~/CS/BAC/';


zer = load([bcs_data_path,'full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT'])
zerSS = load([bcs_data_path,'full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT_onlySS'])

two = load([bcs_data_path,'full16S/data_primersAmit_2MM_databaseOf2012_readLength_76_noACGT'])
twoSS = load([bcs_data_path,'full16S/data_primersAmit_2MM_databaseOf2012_readLength_76_noACGT_onlySS'])

% find the position of bact in zer in the list of two


mx = max(max(cell2mat(zer.sameSeq)),max(cell2mat(two.sameSeq)));
pos_zer = zeros(1,mx);
for i=1:length(zer.sameSeq)
  for j=1:length(zer.sameSeq{i})
    pos_zer(zer.sameSeq{i}(j)) = i;
  end
end

pos_two = zeros(1,mx);
for i=1:length(two.sameSeq)
  for j=1:length(two.sameSeq{i})
    pos_two(two.sameSeq{i}(j)) = i;
  end
end

for i=1:length(two.sameSeq)
  for j=1:length(two.sameSeq{i})
    if two.sameSeq{i}(j)==2
      [i,j]
      pause
    end
  end
end

% for each group (single) in zer - find its corresponding in 780000
[junk_zer ind_zer len_zer] = get_duplicates(pos_zer');
if junk_zer(1)==0
  junk_zer(1) = [];
  ind_zer(1) = [];
  len_zer(1) = [];
end

conversion = cell(1,length(ind_zer));
for i=1:length(ind_zer)
  clear dat
  dat(1,:) = ind_zer{i}; % name in 780000 of current group in 0
  
  % correct:
  %if length(unique(dat(1,:)))<length(dat(1,:))
  %  disp('problem')
  %  pause
  %end

  curr_two = pos_two(ind_zer{i});  
  dat(2,:) = curr_two; % 
  
  if curr_two(1)==0
    curr_two(1) = [];
  end
  
  if length(curr_two)>1 % need to find the header
    pos = zeros(1,length(zerSS.SS.Header_uni{i}));
    for l=1:length(zerSS.SS.Header_uni{i})
      currName = zerSS.SS.Header_uni{i}{l};
      for j=1:length(curr_two)
        for k=1:length(twoSS.SS.Header_uni{curr_two(j)})
          if ~isempty(findstr(twoSS.SS.Header_uni{curr_two(j)}{k},currName))
            pos(l) = j;
          end
        end
        
      end
    end
    % they should be the same
    curr_two = unique(curr_two(pos));
  end
  conversion{i} = curr_two;
end

l = cellfun(@(x) length(x),conversion,'uniformoutput',false);


% conversion2 takes the one with more regions if more mapping to more than 1
conversion2 = cell(1,length(ind_two));
for i=1:length(ind_two)
  conversion2{junk_two(i)} = unique(pos_zer(ind_two{i}));  
  if conversion2{junk_two(i)}(1)==0
    conversion2{junk_two(i)}(1) = [];
  end
  
  if length(conversion2{junk_two(i)})>1 % take the one with the largerst number of regions
    Y = conversion2{junk_two(i)};
    x = zeros(1,length(Y));
    for j=1:length(Y)
      [x(j)] = length(find(zer.seq_readLength76(Y(j),:)~=char(1))); 
    end
    [junk,curr_ind] = max(x);
    conversion2{junk_two(i)} = Y(curr_ind);
  end
end




%for i=1:length(conversion)
%  if length(conversion{i})>1 %&& ~isempty(find(conversion{i}))
%    i
%    pause
%  end
%end

11083
26215
pause
newM = zer.M;
for i=1:6
  i
  % find a row the appears in both
  [junk,i1,i2] = intersect(zer.values{i},two.values{i},'rows');
  
  for j=1:length(i1) % for each matching line
    if mod(j,1000)==0
      [i,j,length(i1)]
    end
    
    % find samples in zer
    samples_zer = find(zer.M{i}(i1(j),:));
    
    for k=1:length(samples_zer)
      ind_in_zer = find(zer.M{i}(:,samples_zer(k)));
      counter_ind_in_zer = setdiff(ind_in_zer,i1(j));clear ind_in_zer
      
      if length(counter_ind_in_zer)~=1
        disp('problem')
        pause
      end
      
      %[junk,ind_in_two] =         intersect(two.values{i},zer.values{i}(i1(j),:),'rows');
      [junk,counter_ind_in_two] = intersect(two.values{i},zer.values{i}(counter_ind_in_zer,:),'rows');
      
      samples_two_for_ind =     find(two.M{i}(i2(j),:));
      samples_two_for_counter = find(two.M{i}(counter_ind_in_two,:));
      
      added_samples_in_two = intersect(samples_two_for_ind,samples_two_for_counter);
      % the new samples must be amplified in at least one region, otherwise appear as empty
      added_samples_in_zer = unique(cell2mat(conversion2(added_samples_in_two)));
      
      if ~isempty(find(added_samples_in_zer==275098))
        newM{i}(:,275098)
        [i1(j) counter_ind_in_zer]
        pause
      end
      
      if ~isempty(added_samples_in_zer)
        [i,j,k]
        %pause
        newM{i}(i1(j),added_samples_in_zer) = 1;
        newM{i}(counter_ind_in_zer,added_samples_in_zer) = 1;
        
        %if newM{i}(28616,275098)==1
        %  1
        %  pause
        %end
        
        if ~isempty(find(sum(newM{i},1)>2))
          disp('problem')
          pause
        end
      
      end
      
    end
  end
end

newM = zer.newM;

README = ['conversion is not a function. there are lines with more than a single mapping. Hence should check results according to this - I am taking ' ...
          'the one with more regions. should always check this.']
save([bcs_data_path,'full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT_conversionFrom_2MM'],'newM','conversion','conversion2','README')
