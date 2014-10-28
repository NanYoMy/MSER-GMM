function write_infile(filename,source_matrix )
%TEST Summary of this function goes here
%  write matrix to file
    fid=fopen(filename,'w+');
    [x,y]=size(source_matrix);
    for i=1:x
        for j=1:y-1
            fprintf(fid,'%d\t',source_matrix(i,j));
        end
        fprintf(fid,'%d\r\n',source_matrix(i,y));%每一行回车\n
    end
    fclose(fid);
end
