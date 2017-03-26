function line=scom(fid)

line=fgetl(fid);
while (size(line,2)==0 || line(1,1)=='!')
    if feof(fid)
        error('Impossibile individuare riga non commentata')
    end
    line=fgetl(fid);
end 

