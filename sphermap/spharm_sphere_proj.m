
% SPHARM_SPHERE_PROJ
% project a triangulated surface on a sphere
%
% [VERspher, poles, dateline] = spharm_sphere_proj(VER,ITRI);
% 
% input:
%   VER,ITRI:  triangulated surface (must be closed and simply connected)
% 
% outputs:
%   VERspher:  location of the nodes on the sphere (connectivity is the same)
%   poles:     indices of the two poles
%   dateline:  indices of a line from one pole to the other


scriptname = mfilename('fullpath');
slash = find(scriptname == '/');
dirname = scriptname(1:slash(end)-1);
addpath(dirname);

fid = fopen('_sphermap.config');
fname = fgetl(fid);
niter = round(str2double(fgetl(fid)));
fclose(fid);



confs.MeshGridSize = 50;
confs.MaxSPHARMDegree = 6;
confs.Tolerance = 2;
confs.Smoothing = 2;
confs.Iteration = niter;
confs.LocalIteration = 10;
confs.t_major = 'x';
confs.SelectDiagonal = 'ShortDiag';


[VER,ITRI] = loadtri(fname);
[VER2, poles, dateline] = spharm_init_param_cald(VER, ITRI);
[VER3, ITRI3, VERspher] = spharm_smooth_cald(VER, ITRI, VER2, confs);

save(fname, '-ascii', 'VERspher');

