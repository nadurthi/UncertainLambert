function A=convertstrrow2doublevector(ss)
T=strsplit(ss,',');
A=zeros(1,length(T));
for i=1:length(T)
    A(i) = str2double(T{i});
end
