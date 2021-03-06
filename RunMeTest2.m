cd 'C:\Users\Mike\Documents\MATLAB\Radius Measurement\Automation attempts\2019.01.11 revision for fully automating\lil folder'
info = dir;
a = length(info);
file_num = 3:a;
b = length(file_num);
names = cell(1,b);
for m = 1:b
    n = m+2;
    temp_struct = info(n);
    name = temp_struct.name;
    names{m} = name;
end
cd 'C:\Users\Mike\Documents\MATLAB\Radius Measurement\Automation attempts\2019.01.11 revision for fully automating'

name = names{1};
Ipack = MasterCrop(name);
wavelengthLimits = Ipack.limits;
clear Ipack

c = length(names);
indexIN = find(contains(names,'IN'));
d = length(indexIN);
indexOUT = find(contains(names,'OUT'));
e = length(indexOUT);
IN = zeros(d,2);
OUT = zeros(e,2);
for n = indexIN
    name = names{n};
    disp(name);
    position = find(indexIN == n);
    Ipack = MasterTest2(name,wavelengthLimits);
    d = Ipack.diameter_odd_even;
    num = regexp(name,'\d');
    enc = name(min(num):max(num));
    enc = str2double(enc);
    neg = strfind(name,'-');
    if isempty(neg)
    else
        enc = enc*-1;
    end
    IN(position,1) = d;
    IN(position,2) = enc;
    clear Ipack
end
for n = indexOUT
    name = names{n};
    disp(name);
    position = find(indexOUT == n);
    Ipack = MasterTest2(name,wavelengthLimits);
    d = Ipack.diameter_odd_even;
    num = regexp(name,'\d');
    enc = name(min(num):max(num));
    enc = str2double(enc);
    neg = strfind(name,'-');
    if isempty(neg)
    else
        enc = enc*-1;
    end
    OUT(position,1) = d;
    OUT(position,2) = enc;
    clear Ipack
end
