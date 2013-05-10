function currNumber=findNumberFromHeader(Header)

currNumber = Header;
a = find(currNumber==' ');
currNumber = currNumber(1:a(1)-1);
currNumber = str2num(currNumber);
