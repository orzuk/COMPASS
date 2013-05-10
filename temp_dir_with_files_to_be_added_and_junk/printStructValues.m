function printStructValues(auxData,structName)

fn = fieldnames(auxData);
fprintf('%%%%%%%%%%%\n');
fprintf('%s:\n',structName);
for i=1:length(fn)
  if ~strcmp(fn{i},'ErrorStruct')& ~strcmp(fn{i},'numProcessors') 
    if isnumeric(auxData.(fn{i})) 
      fprintf('%s: %d\n',fn{i},num2str(auxData.(fn{i})));
    else
      fprintf('%s: %s\n',fn{i},auxData.(fn{i}));
    end
    
  elseif strcmp(fn{i},'ErrorStruct') % ErrorStruct
    fprintf('\n\n');
    if auxData.addNoiseFlag
      fprintf('ErrorStruct was considered\n\n\n')
    end
  end % if ErrorStruct
  
end    
fprintf('end data\n')
fprintf('%%%%%%%%%%%%\n')
