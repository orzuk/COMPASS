function out=get_mean_posalign_Amnon_MergeForwardReverse(read, database1,database_size)


% all are forward
k=strfind(database1',read);

if isempty(k)
  disp('all should be aligned and forward. error. get_mean_posalign_Amnon_MergeForwardReverse.m')
  return
end

[I,J]=ind2sub(database_size,k);

out = struct;
out.I = I;
out.J = J;





