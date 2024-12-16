%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherical Harmonic Modeling and Analysis Toolkit (SPHARM-MAT) is a 3D 
% shape modeling and analysis toolkit. 
% It is a software package developed at Shenlab in Center for Neuroimaging, 
% Indiana University (SpharmMat@gmail.com, http://www.iupui.edu/~shenlab/)
% It is available to the scientific community as copyright freeware 
% under the terms of the GNU General Public Licence.
% 
% Copyright 2009, 2010, ShenLab, Center for Neuroimaging, Indiana University
% 
% This file is part of SPHARM-MAT.
% 
% SPHARM-MAT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SPHARM-MAT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SPHARM-MAT. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vertices, faces, sph_verts] = spharm_smooth_cald(vertices, faces, sph_verts, confs)



info = [confs.MeshGridSize confs.MaxSPHARMDegree confs.Tolerance confs.Smoothing confs.Iteration ... 
    confs.LocalIteration];
reso = confs.LocalIteration; 
tmajor=confs.t_major; 

% initialization
verts = sph_verts;

% calculate relative areas of triangles on object surface net  
% vertices made a matrix
[obj_area, count] = calc_triangle_areas(vertices, faces, 'object', 'triangle');

% calculate area scaling ratio between object mesh and parameter mesh 
cal_asr(verts, faces, obj_area);

vertnum = size(vertices,1); facenum = size(faces,1);
disp(sprintf('  --- %+4d = %d vertices * 2 - %d faces ', vertnum*2-facenum, vertnum, facenum));

% remove bad areas
verts = call_locsmo(verts,vertices,faces-1,10,reso);

curinfo = info; premstch = inf; count = 0;
for j = 1:info(5)
    disp(sprintf('-------- SPHARM deg %d: %d of %d ---------',curinfo(2),j,info(5)));
    % smooth by interpolation
    [verts, mstch] = smooth(verts,faces,curinfo,obj_area,tmajor);
    cal_asr(verts, faces, obj_area);
    verts = call_locsmo(verts,vertices,faces-1,10,reso);
    disp(sprintf('mstch=%f',mstch));
    
    if (premstch-mstch)>0.001
        count=0; premstch=mstch;
    else
        count=count+1;
    end
    
    if (mstch<1.001 | count>=3) break; end
end

sph_verts = verts;

metric = get_metric('_sphermap.metric',11);
metric(end+1)=mstch;
metric(end+1)=j;

measure.vals = metric';
measure.name = {'bad_areas','sqerrobj','sqerrprm','adcobj','adcprm',...
    'L2obj','L2prm','L2bobj','L2bprm','Linf','Linfb','meanstch','iter'};


return;
endfunction

%
% Get constraints
%

function angles = get_angles(vs,faces)

fnum = size(faces,1);
d = size(vs);

angles = [];
for j = 1:4
    A = vs(faces(:,j),:);
    B = vs(faces(:,mod(j,4)+1),:);
    C = vs(faces(:,mod(j-2,4)+1),:);
    y = A(:,1).*B(:,2).*C(:,3) - A(:,1).*B(:,3).*C(:,2) + ...
        A(:,2).*B(:,3).*C(:,1) - A(:,2).*B(:,1).*C(:,3) + ...
        A(:,3).*B(:,1).*C(:,2) - A(:,3).*B(:,2).*C(:,1);
    x = B(:,1).*C(:,1) + B(:,2).*C(:,2) + B(:,3).*C(:,3) - ...
       (A(:,1).*C(:,1) + A(:,2).*C(:,2) + A(:,3).*C(:,3)).* ...
       (A(:,1).*B(:,1) + A(:,2).*B(:,2) + A(:,3).*B(:,3));
    angles(:,j) = atan2(y,x);
end
% fix angle range (0 - 2pi)
ind = find(angles<0);
angles(ind) = angles(ind) + 2*pi;

return;
endfunction

%**************************************************************************
%
% get metric from file metric
%

function metric = get_metric(fname,n)

% read from metric file
disp(sprintf('read verts from %s ...',fname));
fid = fopen(fname,'r');
metric = fread(fid, n, 'double'); 
fclose(fid);

if exist(fname)
    delete(fname);
end

return;
endfunction

%**************************************************************************
%
% call LocalSmoothing
%

function verts = call_locsmo(param_vs,obj_vs,faces,extents,reso)

infilename = tempname;
outfilename = tempname;

vnum = size(param_vs,1); fnum = size(faces,1);
% save to infile
disp(sprintf('save param_vs, obj_vs, faces, extents, reso to infile ...'));
fid = fopen(infilename,'wb+');
fwrite(fid, [reso vnum fnum], 'int');
fwrite(fid, extents, 'double');
fwrite(fid, param_vs, 'double');
fwrite(fid, obj_vs, 'double');
fwrite(fid, faces, 'int');
fclose(fid);

scriptname = mfilename('fullpath');
slash = find(scriptname == '/');
dirname = scriptname(1:slash(end));

cmd = [dirname, 'spharm_local_smoothing ', infilename, ' ', outfilename];
unix(cmd);
%spharm_local_smoothing(infilename, outfilename); 

% read from outfile
disp(sprintf('read verts from outfile ...'));
fid = fopen(outfilename,'r');
verts = fread(fid, vnum*3, 'double'); 
verts = reshape(verts,size(param_vs));
fclose(fid);

% remove remeshout
if exist(outfilename,'file')
    delete(outfilename);
end
if exist(infilename,'file')
    delete(infilename);
end

return;
endfunction

%**************************************************************************

%
% smooth spherical mesh defined by verts and faces
%

function [verts, mstch] = smooth(verts,faces,info,obj_areas,tmajor)

meshsize = info(1); shd = info(2); tole = info(3); smft = info(4);

% area scaling ratio
[org_areas, count] = calc_triangle_areas(verts, faces, 'parameter', 'triangle');
org_areas = org_areas./obj_areas;

% smooth areas a little bit
areas = org_areas.^(1/smft);

% deal with negative areas and super large areas (not done yet)
ind = find(areas<=0);
if ~isempty(ind)
    disp([sprintf('!!! === !!! %d bad relative areas:',length(ind)) sprintf(' %f',areas(ind))]);
end

% find the center of each face, assign radius for it (using area of the face)
cents = (verts(faces(:,1),:) + verts(faces(:,2),:) + verts(faces(:,3),:))/3;

% create spharm descriptor for a shape drived from centers plus their radiuses
d = shd; 
[fvec, d, Z, name_temp] = create_SPHARM_des_LSF(areas, [], cents, d, '', '');

ind = 1:(d+1)^2; radius = real(Z(:,ind)*fvec(ind)); radius = radius.^smft;

[ma, ma_ind] = max(abs(log(radius)));
ma_ind = ma_ind(1); ma_cent = mean(verts(faces(ma_ind,:),:));

[phi,theta] = cart2sph(ma_cent(1), ma_cent(2), ma_cent(3));

Ra = rotate_mat(0, pi/2-theta, 0)*rotate_mat(0, 0, phi); % rotate
Rb = rotate_mat(0, 0, -phi)*rotate_mat(0, -pi/2+theta, 0); % rotate back

verts = verts*Ra'; 

% find the center of each face, assign radius for it (using area of the face)
cents = (verts(faces(:,1),:) + verts(faces(:,2),:) + verts(faces(:,3),:))/3;

% create spharm descriptor for a shape drived from centers plus their radiuses
d = shd; 
%[fvec,d,Z] = spharm_vec(areas, cents, d);
[fvec, d, Z, name_temp] = create_SPHARM_des_LSF(areas, [], cents, d, '', '');

% create mesh for interpolation
gsize = 2/meshsize; ti = -1:(gsize/2):1; n = length(ti);
[PHIs,THETAs] = meshgrid(ti,ti*pi); PHIs = asin(PHIs)+pi/2; THETAs = THETAs+pi;

% each grid consists of 4 grids generated before so that the centroid can be easily located
gind = 2:2:length(ti); % length(ti) should be an odd number (grid is indexed by its center)
gPHIs = PHIs(gind,gind); gTHTs = THETAs(gind,gind);
Zm = spharm_basis(d,gPHIs(:),gTHTs(:));
garea = real(Zm(:,ind)*fvec(ind)); garea = reshape(garea,size(gPHIs)); garea = garea.^smft;

% get each grid height: gheight*garea=1
gheight = ones(size(garea))./garea;
gmin = min(gheight(:)); bigind = find(gheight>gmin*tole); gheight(bigind) = gmin*tole;
disp(sprintf('Adjust %d big values to %d * %f (gmin) = %f',length(bigind),tole,gmin,gmin*tole));

% prepare for interpolation
mind = 1:2:length(ti); X = PHIs(mind,mind); Y = THETAs(mind,mind);
gheight_0 = zeros(size(gheight)+1);
gheight_0(2:end,2:end) = gheight; gheight = gheight_0;

% x major
tt = sum(gheight(:));
tvals = cumsum(sum(gheight,1))/tt; tvals = tvals(ones(1,length(tvals)),:);
lucumsum = cumsum(cumsum(gheight,1),2)/tt;
pvals = zeros(size(lucumsum)); nzind = find(tvals~=0);
pvals(nzind) = lucumsum(nzind)./tvals(nzind);
pvals(:,1) = pvals(:,2); pvals = (pvals-0.5)*pi*2;
tvals = real(asin(tvals*2-1))+pi/2; % tvals need to be adjusted

stretch = gheight(find(gheight~=0)); stretch = max(stretch,1./stretch); mstch = mean(stretch(:));
titstr = sprintf('spharm stretch: mean %f std %f',mstch,std(stretch(:))); disp(titstr);

% do interpolation for each vertices
[ps, ts] = cart2sph(verts(:,1),verts(:,2),verts(:,3));
ts = pi/2-ts;
in = find(ps<0); ps(in) = ps(in)+2*pi; 

new_ts = interp2(X,Y,tvals,ts,ps);
new_ps = ps;

[verts(:,1),verts(:,2),verts(:,3)] = sph2cart(new_ps, pi/2-new_ts,1);

verts = verts*Rb';

return;
endfunction


%**************************************************************************

%
% calculate area scaling ratio between object mesh and parameter mesh
%

function cal_asr(verts, faces, obj_area)

% calculate relative areas of spherical triangles in parameter space
[par_area, count] = calc_triangle_areas(verts, faces, 'parameter', 'triangle');

asr = par_area./obj_area; asr2 = 1./asr;
stretch = max(asr,asr2);

disp(sprintf('==> reduce %f, enlarge %f, objmean %f, prmmean %f',...
              max(1./asr),max(asr),sum(stretch.*obj_area),sum(stretch.*par_area)));

return;
endfunction


% ============================================
%
% Goal: Create canonical spherical harmonic bases
%
% Li Shen 
% 04/11/2002 - create
% 10/15/2002 - rename and modify
% 11/03/2008 - renamed by Sungeun Kim.

function Z = calculate_SPHARM_basis(vs, degree)

[PHI,THETA] = cart2sph(vs(:,1),vs(:,2),vs(:,3));
ind = find(PHI<0);
PHI(ind) = PHI(ind)+2*pi;
THETA = pi/2-THETA;
vertnum = size(THETA,1);

Z = spharm_basis(degree,THETA,PHI); 

return;
endfunction


function Z = spharm_basis(max_degree,theta,phi)

Z = []; vnum = size(theta,1);

% save calculations for efficiency
for k = 0:(2*max_degree)
    fact(k+1) = factorial(k);
end
for m = 0:max_degree
    exp_i_m_phi(:,m+1) = exp(i*m*phi);
    sign_m(m+1) = (-1)^(m);
end

for n = 0:max_degree

	% P = legendre(n,X) computes the associated Legendre functions of degree n and 
	% order m = 0,1,...,n, evaluated at X. Argument n must be a scalar integer 
	% less than 256, and X must contain real values in the domain -1<=x<=1.
	% The returned array P has one more dimension than X, and each element
	% P(m+1,d1,d2...) contains the associated Legendre function of degree n and order
	% m evaluated at X(d1,d2...).

    Pn = legendre(n,cos(theta'))';
    
    posi_Y = [];
    nega_Y = [];
    
    m= 0:n;
    v = sqrt(((2*n+1)/(4*pi))*(fact(n-m+1)./fact(n+m+1)));
    v = v(ones(1,vnum),:).*Pn(:,m+1).*exp_i_m_phi(:,m+1);
    posi_Y(:,m+1) = v; % positive order;
    nega_Y(:,n-m+1) = sign_m(ones(1,vnum),m+1).*conj(v); % negative order
    
    Z(:,end+1:end+n) = nega_Y(:,1:n);
    Z(:,end+1:end+n+1) = posi_Y(:,1:(n+1));
end

return;
endfunction




%
% create spherical harmonic descriptor
%

function [fvec, deg, Z, new_name] = create_SPHARM_des_LSF(vertices, faces, sph_verts ,maxDeg, filename, OutDirectory)

if isempty(vertices) | isempty(sph_verts)
    if ~isempty(filename)
        load(filename);
    else
        disp('There is no useful information');
        return;
    end
end

[pa, na, ex] = fileparts(filename);
new_name = '';

if ~exist('vertices', 'var') | ~exist('sph_verts', 'var')
    disp('One or more of vertices, and spherical vertices are missing');
    return;
end

vertnum = size(sph_verts,1);

max_d = maxDeg;
% Note that degree 'd' we want to use depends on the vertnum 
% The total number of unknowns is (d+1)*(d+1)
% The total number of equations is vertnum
% We want equ_num >= unk_num
deg = max(1, floor(sqrt(vertnum)*1/2));
deg = min(deg, max_d);
disp(sprintf('Use spharm up to %d degree (vec_len=%d).',deg,(deg+1)^2));

Z = calculate_SPHARM_basis(sph_verts, deg); 

[x,y] = size(Z);
disp(sprintf('Least square for %d equations and %d unknowns',x,y));

% Least square fitting
% fvec = Z\vertices;   %This does not work as it is expected to work in certain environment
for i=1:size(vertices,2)
    fvec(:,i) = Z\vertices(:,i);
end

if ~isempty(filename)
    if ~isempty(OutDirectory)
        new_name = sprintf('%s/%s_LSF_des.mat', OutDirectory, na(1:end-4));
    else
        new_name = sprintf('%s/%s_LSF_des.mat', pa, na(1:end-4));
    end
    if exist(new_name, 'file')
        prompt = {'Enter new filename:'};
        dlg_title = 'New File Name';
        num_lines = 1;
        def = {new_name};
        answer = inputdlg(prompt,dlg_title,num_lines,def);    
        new_name = answer{1};
    end
    save(new_name, 'vertices', 'faces', 'sph_verts', 'fvec');
end

return;
endfunction


%
% Rotate around x, y, z (counterclockwise when looking towards the origin)
%

function R = rotate_mat(x, y, z)
% It is assumed that this matrix is muplied from the hind, rotating objects
% along x-axis, y-axis, and z-axis in this order.

Rx = [      1       0       0; ...
            0  cos(x) -sin(x); ...
            0  sin(x)  cos(x)];

Ry = [ cos(y)       0  sin(y); ...
            0       1       0; ...
      -sin(y)       0  cos(y)];

Rz = [ cos(z) -sin(z)       0; ...
       sin(z)  cos(z)       0; ...
            0      0        1];

R = Rz*Ry*Rx;

return;
endfunction



%
% calculate relative areas of triangles on object surface net
%

function obj_area = cal_obj_area(vertices,faces)

A = faces(:,1); B = faces(:,2); C = faces(:,3);
a = sqrt(sum(((vertices(A,:)-vertices(B,:)).^(2))'))';
b = sqrt(sum(((vertices(B,:)-vertices(C,:)).^(2))'))';
c = sqrt(sum(((vertices(C,:)-vertices(A,:)).^(2))'))';
s = (a+b+c)/2;
obj_area = sqrt(s.*(s-a).*(s-b).*(s-c));
obj_area = obj_area/sum(obj_area);

return;
endfunction


%
% calculate relative areas of spherical triangles in parameter space
%

function par_area = cal_par_area(vs,faces)

angles = [];
for j = 1:3
    % note that the order of A B C is clockwise (see 08-22-02.htm notes)
    A = vs(faces(:,j),:);
    B = vs(faces(:,mod(j,3)+1),:);
    C = vs(faces(:,mod(j-2,3)+1),:);
    y = A(:,1).*B(:,2).*C(:,3) - A(:,1).*B(:,3).*C(:,2) + ...
        A(:,2).*B(:,3).*C(:,1) - A(:,2).*B(:,1).*C(:,3) + ...
        A(:,3).*B(:,1).*C(:,2) - A(:,3).*B(:,2).*C(:,1);
    x = B(:,1).*C(:,1) + B(:,2).*C(:,2) + B(:,3).*C(:,3) - ...
       (A(:,1).*C(:,1) + A(:,2).*C(:,2) + A(:,3).*C(:,3)).* ...
       (A(:,1).*B(:,1) + A(:,2).*B(:,2) + A(:,3).*B(:,3));
    angles(:,j) = atan2(y,x); 
end
ind = find(angles<0);
angles(ind) = angles(ind) + 2*pi;
par_area = sum(angles')' - pi;
par_area = par_area/(4*pi);

return;
endfunction



function [areas, count] = calc_triangle_areas(verts, faces, space, focus)

if strcmp(space, 'object')
    tri_areas = cal_obj_area(verts, faces);
elseif strcmp(space, 'parameter')
    tri_areas = cal_par_area(verts, faces);
end

if strcmp(focus, 'vertex')
    areas = zeros(length(verts),1);
    count = zeros(length(verts),1);

% find all incident triangles upon each vertex and accumulate the areas of the triangles   
    for i=1:size(faces,2)
        [u1, m1, n1] = unique(faces(:,i));
        areas(u1) = areas(u1) + accumarray(n1, tri_areas);
        count(u1) = count(u1) + accumarray(n1, ones(length(n1),1));
    end
    
elseif strcmp(focus, 'triangle')
    areas = tri_areas;
    count = ones(length(tri_areas),1);
end

return;
endfunction