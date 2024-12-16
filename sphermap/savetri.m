function savetri(fn,VER,ITRI,meth)
% SAVETRI 
% store vertices and triangles in a .tri file
%
% savetri(filename,VER,ITRI);

% nodes now stored in 8.6f format

if nargin<4
    meth = 0;
end

f = fopen(fn, 'w');
[nver dim]=size(VER);
fprintf(f,'%d\n ',nver);

if meth
    for i=1:nver;
        fprintf(f,'%5d %8.6f %8.6f %8.6f\n',i,VER(i,1:3));
    end
else
    X = [(1:nver)' VER(:,1:3)];
    fprintf(f,'%5d %8.6f %8.6f %8.6f\n',X');
end

[ntri dim]=size(ITRI);
fprintf(f,'%d\n ',ntri);
if meth
    for i=1:ntri;
        fprintf(f,'%5d %5d %5d %5d\n',i,ITRI(i,1:3));
    end
else
    X = [(1:ntri)' ITRI(:,1:3)];
    fprintf(f,'%5d %5d %5d %5d\n',X');
end
%sprintf('\ntriangle specs written to file: %s\n',fn);
fclose(f);
