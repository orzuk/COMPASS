function auxData=checkValues(auxData,defValues)
%keyboard
fn = fieldnames(defValues);
for i=1:length(fn)
  if ~isfield(auxData,fn{i})
    auxData = setfield(auxData,fn{i},getfield(defValues,fn{i}));
  end
end
