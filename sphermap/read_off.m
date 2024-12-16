function [ver, tri] = read_off(filename);
% READ_OFF
% load a triangulated surface stored in the (simplified) OFF format
%
% [VER, ITRI] = read_off(filename);

fid = fopen(filename, 'rt');
if fid == -1
    error('unable to open file');
end

fgetl(fid); % OFF

n = fscanf(fid, '%d', 3);
nv = n(1);
nt = n(2);

ver = fscanf(fid, '%f', [3, nv])';
tri = fscanf(fid, '%d', [4, nt])';
tri = tri(:, 2:end) + 1;

fclose(fid);
