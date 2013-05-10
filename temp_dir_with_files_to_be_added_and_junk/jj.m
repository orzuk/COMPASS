  ampData = [];
  lab750 = [];
for j=1:length(inds{i})  
    
    % amplified name and number of errors - with minor
    ampData{j} = num2str([res{inds{i}(j)}{1}.Full.Numbers length(find(align_mat(i,:)==0)) length(find(align_vec_mat(i,:)==0))]);
    ampData{j} = delSpaces(ampData{j});
    ampData{j} = ['{\color{magenta} ',ampData{j},'}'];
    data_ampData{j} = [res{inds{i}(j)}{1}.Full.Numbers length(find(align_mat(i,:)==0)) length(find(align_vec_mat(i,:)==0))];
    
    % 750
    tmp_lab750{j} = [];
    
    clear data750 data450 number750
    for k=1:length(res{inds{i}(j)})
      if ~isempty(res{inds{i}(j)}{k}.s750.Sequence)
        number750{k} = intersect(res{inds{i}(j)}{k}.s750.Numbers,tmp_ind_W_Solexa);
        data750{k} = [ overlapData{inds{i}(j)}{k}.s750.length overlapData{inds{i}(j)}{k}.s750.numErrorsMajor overlapData{inds{i}(j)}{k}.s750.numErrorsMinor]; 
      else
        number750{k} = [];
        data750{k} = [];
      end
    end
    
    taken750 = [];
    tmp_lab750 = [];
    for k=1:length(res{inds{i}(j)})
      if ~isempty(number750{k}) && isempty(intersect(taken750,number750{k}))
        tmp_lab750 = [tmp_lab750;';{\color{blue} {\bf ',num2str(number750{k}),' [',num2str(data750{k}),']}}'];
        taken750 = [taken750,number750{k}];
      end
    end
    
    lab750 = [lab750,tmp_lab750];
  
  end