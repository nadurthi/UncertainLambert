function [efmoutput,L,ind]=parse_cpp_struct(filename_meta, filename_data)


fid = fopen(filename_meta);
i=1;
L=cell(1,7);
while(1)
    line = fgetl(fid);
    if line == -1
        break
    end
    L{i}=line;
    i=i+1;
    
end
if L{2} == 'NORMALIZED'
    efmoutput.normalized=true;
end
T=strsplit(L{3},',');
efmoutput.inputconfig.TFt0 = str2double(T{1}) ;
efmoutput.inputconfig.TFtf = str2double(T{2}) ;
efmoutput.inputconfig.TFgridlen = str2double(T{3}) ;

T=strsplit(L{4},',');
efmoutput.inputconfig.dept0 = str2double(T{1}) ;
efmoutput.inputconfig.deptf = str2double(T{2}) ;
efmoutput.inputconfig.depgridlen = str2double(T{3}) ;

efmoutput.inputconfig.DepTime  = linspace(efmoutput.inputconfig.dept0,efmoutput.inputconfig.deptf,efmoutput.inputconfig.depgridlen);                  % Vector of Departure Times
efmoutput.inputconfig.TimeFlights = linspace(efmoutput.inputconfig.TFt0,efmoutput.inputconfig.TFtf,efmoutput.inputconfig.TFgridlen);                % Vector of Times of Flight

T=strsplit(L{5},',');
efmoutput.inputconfig.rD =[str2double(T{1});str2double(T{2});str2double(T{3})];
efmoutput.inputconfig.vD =[str2double(T{4});str2double(T{5});str2double(T{6})];

T=strsplit(L{6},',');
efmoutput.inputconfig.rA =[str2double(T{1});str2double(T{2});str2double(T{3})];
efmoutput.inputconfig.vA =[str2double(T{4});str2double(T{5});str2double(T{6})];

T=strsplit(L{7},',');
efmoutput.inputconfig.dep_epoch = str2double(T{1}) ;
efmoutput.inputconfig.arr_epoch = str2double(T{2}) ;

fclose(fid);


%%
fid = fopen(filename_data);
ind=0;
L=cell(1,10000);
line = fgetl(fid);
while(1)
    
    if line == -1
        break
    end
    if contains(line,'INDEX = ')
       line = fgetl(fid);
       ind=ind+1;
       continue;
    end
    
    if strcmp(line,'xA')
        line = fgetl(fid);
        L{ind}.xA=line;
        line = fgetl(fid);
        continue;
    end
    
    if strcmp(line,'xB')
        line = fgetl(fid);
        L{ind}.xB=line;
        line = fgetl(fid);
        continue;
    end
    if strcmp(line,'t_des')
        line = fgetl(fid);
        L{ind}.t_des=line;
        line = fgetl(fid);
        continue;
    end
    
    if strcmp(line,'LP')
        line = fgetl(fid);
        if isfield(L{ind},'LP') 
            L{ind}.LP{2}=line;
        else
            L{ind}.LP{1}=line;
        end
        line = fgetl(fid);
        continue;
    end
    
    if strcmp(line,'LCPV')
        P={};k=1;
        while(1)
            line = fgetl(fid);
            if strcmp(line,'LCPV')==1 || strcmp(line,'LSV')==1
                break;
            end
            if isempty(line)==0
                P{k}=line;
            end
            k=k+1;
        end
        if isfield(L{ind},'LCPV') 
            L{ind}.LCPV{2}=P;
        else
            L{ind}.LCPV{1}=P;
        end
        continue;
    end
    
    if strcmp(line,'LSV')
        P={};k=1;
        while(1)
            line = fgetl(fid);
            if strcmp(line,'LSV')==1 || strcmp(line,'LTV')==1
                break;
            end
            if isempty(line)==0
                P{k}=line;
            end
            k=k+1;
        end
        if isfield(L{ind},'LSV') 
            L{ind}.LSV{2}=P;
        else
            L{ind}.LSV{1}=P;
        end
        continue;
    end
    
     if strcmp(line,'LTV')
        P={};k=1;
        while(1)
            line = fgetl(fid);
            if strcmp(line,'LTV')==1 || strcmp(line,'ENDINDEX')==1
                break;
            end
            if isempty(line)==0
                P{k}=line;
            end
            k=k+1;
        end
        if isfield(L{ind},'LTV') 
            L{ind}.LTV{2}=P;
        else
            L{ind}.LTV{1}=P;
        end
        continue;
     end
     
     
     ind
     line = fgetl(fid);
     
    
end
fclose(fid)

