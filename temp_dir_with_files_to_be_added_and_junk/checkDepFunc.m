function checkDepFunc(fileName,fullPathFlag,sepNam)

a = depfun(fileName);
if ~exist('sepName')
    sepName = ' '; 
end

add = [];
for i=2:length(a)
   
    if ~isempty(findstr(a{i},'shental'));
        if fullPathFlag
            add = [add,sepName,a{i}];  
        else
            b = findstr(a{i},'\');
            add = [add,sepName,a{i}(b(end)+1:end)];  
        end
        
       
        
    end
end

if 1==2
  cd D:\shental\CS\mFilesBAC
  add_part0=checkDepFunc('distributeBAC_generalOrFourth',1);
  dos(['cp ',add_part0,' D:\shental\CS\tmp\']);
  
  cd D:\shental\CS\tmp
  add_part0_withSep=checkDepFunc('distributeBAC_generalOrFourth.m',0,' -a ');
  
  w_part0 = ['mcc -m distributeBAC_generalOrFourth.m ',add_part0_withSep(2:end)];
end


