function tmp_lab750=delSpaces(tmp_lab750)

a = find(isspace(tmp_lab750));
b = a(find(diff(a)==1));
tmp_lab750(b) = [];