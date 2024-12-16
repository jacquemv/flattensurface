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

function [sph_verts, poles, dateline] = spharm_init_param_cald(vertices, faces)



vertnum = size(vertices,1);
facenum = size(faces,1);

uvertnum=length(unique(vertices,'rows'));

% find a subset of unique vertices
if (uvertnum ~= vertnum)
    vs = vertices;
    vertices = unique(vs, 'rows');
    vertnum = size(vertices, 1);
    fs = faces;
    faces = zeros(size(fs));
    for i=1:size(fs,2)
        vf=vs(fs(:,i),:);
        [tf, lc] = ismember(vf,vertices, 'rows');
        faces(:,i) = lc;
    end
    faces = unique(faces,'rows');
    facenum = size(faces,1);
    clear vs vf tf lc fs;
end

% create edges
edges = create_edges(faces);
disp(sprintf('%+4d = %d (%d) vertices * 2 - %d faces; %d edges', ...
              vertnum*2-facenum, vertnum, length(unique(faces)), facenum,size(edges,1)));

[AM, WAM] = adjacency_mat(vertices,edges);

[obj_area, count] = calc_triangle_areas(vertices, faces, 'object', 'triangle');
obj_area = obj_area/sum(obj_area);    

[poles, direct2nd] = determine_poles(vertices); % north 1 south 2
    
[verts1,landmarks,dateline] = init_param_wt(vertices,faces,WAM,AM,poles,direct2nd); % initial parameterization
[area_orig, count_orig] = calc_triangle_areas(vertices, faces, 'object', 'triangle');
[area_dist1, count_dist1] = calc_triangle_areas(verts1, faces, 'parameter', 'triangle');
ca = area_dist1./area_orig;
expansion = max(ca, 1./ca);

adc1 = sum(expansion.*area_dist1)/sum(area_dist1); % average case

verts2 = verts1*rotate_mat(pi/2, 0, 0)';
verts2 = recal_longitude(verts2,faces,edges,WAM,AM);
verts2 = verts2*rotate_mat(-pi/2, 0, 0)';
[area_dist2, count_dist2] = calc_triangle_areas(verts2, faces, 'parameter', 'triangle');
ca = area_dist2./area_orig;
expansion = max(ca, 1./ca);

adc2 = sum(expansion.*area_dist2)/sum(area_dist2); % average case

verts3 = verts2*rotate_mat(0, 0, pi/2)';
verts3 = recal_longitude(verts3,faces,edges,WAM,AM);
verts3 = verts3*rotate_mat(0, 0, -pi/2)';
[area_dist3, count_dist3] = calc_triangle_areas(verts3, faces, 'parameter', 'triangle');

ca = area_dist3./area_orig;
expansion = max(ca, 1./ca);

adc3 = sum(expansion.*area_dist3)/sum(area_dist3); % average case

[minADC, idxMin] = min([adc1 adc2 adc3]);

if idxMin == 1
    sph_verts = verts1;
elseif idxMin == 2
    sph_verts = verts2;
else
    sph_verts = verts3;
end
    
sph_verts = verts1;

return;
endfunction

%
% create edges based on faces
%

function edges = create_edges(faces)

fnum = size(faces,1);

edges(1:fnum,:)          = [faces(:,1),faces(:,2)];
edges((1:fnum)+fnum,:)   = [faces(:,2),faces(:,3)];
edges((1:fnum)+fnum*2,:) = [faces(:,3),faces(:,1)];
ind = find((edges(:,1)-edges(:,2))<0);
edges = edges(ind,:);

return;
endfunction

%
% calculate adjacency matrix and also weighted one
%
%     AM = n x n node-node adjacency matrix
%         (Note: AM(i,j) = 0   => Arc (i,j) does not exist;
%                AM(i,j) = 1   => Arc (i,j) exists)
%                AM(i,j) = NaN => i==j)
%     WAM = n x n node-node weighted adjacency matrix of arc lengths
%         (Note: WAM(i,j) = 0   => Arc (i,j) does not exist;
%                WAM(i,j) = NaN => Arc (i,j) exists with 0 weight)
%

function [AM, WAM] = adjacency_mat(vertices,edges)

n = size(vertices,1);

% edge weights matrix (weighted by 1 or the length of the edge)
AM = sparse(n,n); WAM = sparse(n,n); d = [n,n];
ind = sub2ind(d,edges(:,1),edges(:,2)); 
AM(ind) = 1; 
AM = max(AM,AM');
WAM(ind) = sqrt(sum((vertices(edges(:,1),:)-vertices(edges(:,2),:))'.^2)); WAM = max(WAM,WAM');

% assign diagonal
ind = sub2ind(d,1:n,1:n); 
AM(ind) = NaN; WAM(ind) = NaN; % diagonal

return;
endfunction

%
% use PCA to determine north pole (big z) and south pole (small z)
%

function [poles, direc] = determine_poles(vertices)

pca_dim = 3;
[pca_ps, pca_b, var_amt, latent] = do_pca(vertices,pca_dim);

pc1v = pca_ps(:,1); % first PC value
min_i = find(pc1v==min(pc1v)); n = length(min_i);
max_i = find(pc1v==max(pc1v))'; m = length(max_i);
disp(sprintf('after pca min_max_number: %d %d',n,m));

% pick the pair which is the furthest away from each other
min_i = min_i(:,ones(1,m));
max_i = max_i(ones(1,n),:);
dist = sum((vertices(min_i(:),:) - vertices(max_i(:),:)).^2,2);
[maxd, maxdi] = max(dist); % pick the max distance
min_i = min_i(maxdi); max_i = max_i(maxdi);

if vertices(min_i,3) > vertices(max_i,3)
    poles = [min_i max_i];
else
    poles = [max_i min_i];
end

direc1st = vertices(poles(2),:)-vertices(poles(1),:);
direc1st = direc1st./norm(direc1st);

pc2v = pca_ps(:,2); % second PC value
pc3v = pca_ps(:,3);
IDX = find(pc3v==0);
if isempty(IDX)
    IDX = find(abs(pc3v)==min(abs(pc3v)));
end
max_j = find(pc2v==max(pc2v(IDX)))'; max_j = max_j(1);
sdirect = vertices(max_j,:)-vertices(poles(1),:);
sdirect = sdirect./norm(sdirect);

direc3rd = cross(direc1st, sdirect);
direc2nd = cross(direc3rd,direc1st);

direc = [direc2nd;direc3rd];

return;
endfunction


%
% initial parameterization (each edge weighted by its length)
%   W is the edge weight matrix (spaital case: adjacency matrix, each edge has a weight of 1)
%

function [sph_verts,landmarks,dateline] = init_param_wt(vertices,faces,WAM,AM,landmarks,direct2nd)

%%%%%%%%%%%%%%%%%%%%%%%%%
% Latitude from Diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up matrix B
vertnum = size(vertices,1); d = [vertnum, vertnum];
B = sparse(vertnum,vertnum);
%ind = find(AM>0); B(ind) = -1./AM(ind); % edges
ind = find(WAM>0); B(ind) = -1./WAM(ind); % edges

ind = sub2ind(d,1:vertnum,1:vertnum);  
B(ind) = -sum(B); % diagonal

% Set up matrix A
A = B;
% nouth pole
A(landmarks(1),:) = 0;
A(landmarks(1),landmarks(1)) = 1;
% sorth pole
A(landmarks(2),:) = 0;
A(landmarks(2),landmarks(2)) = 1;

% Set up constant vector b
b = sparse(vertnum,1);
b(landmarks(1)) = pi;

disp('Solving simultaneous Linear Equations for latitude ...');
theta = A\b;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select dateline
%%%%%%%%%%%%%%%%%%%%%%%%%%
dateline = []; here = landmarks(1);

direct_org = vertices(landmarks(2),:) - vertices(landmarks(1),:);
direct_org = direct_org./norm(direct_org);

% Approach 2. Global bi-secting plane, defined by two principal axes

direct_2nd = direct2nd(1,:);
normal2surf = cross(direct_org,direct_2nd);
normal2surf = normal2surf ./ norm(normal2surf);

while (here~=landmarks(2))
    dateline(end+1) = here;
    
    %(1) Find all neighbor vertices to the current vertex
    nbs = find(AM(here,:)==1); 
    
    %(2) Among all neighbor vertices, find vertices whose latitude are
    %smaller than the current latitude in parametric space
    nbs2 = nbs(find(theta(nbs)<theta(here)));
    
    %(3) Among the selected neighbor vertices from the previous step, find
    %ones which are closest located to a planed, defined by the first two
    %PCA directions
    nb_verts = vertices(nbs2,:)-repmat(vertices(landmarks(1),:),[length(nbs2),1]);
    leng_nb_verts = sqrt(sum((nb_verts.^2),2)); nb_verts = nb_verts ./ leng_nb_verts(:,ones(1,3));
    
    if length(dateline) > 1
        concord = abs(sum(nb_verts.*repmat(normal2surf,[length(nbs2),1]),2));
        [vv,next] = min(concord);  next = find(concord==vv);
        if (length(next) > 1)
            nb_verts2 = nb_verts(next,:);
            concord2 = abs(sum(nb_verts2.*repmat(direct_org,[length(next),1]),2));
            [vvv,next2]=max(concord2);
            next = next(next2);
        end
        here = nbs2(next(1));
    else
        concord = sum(nb_verts.*repmat(direct_org,[length(nbs2),1]),2);
        [vv,next]=max(concord); next = find(concord==vv);
        if (length(next) > 1)
            nb_verts2 = nb_verts(next,:);
            concord2 = abs(sum(nb_verts2.*repmat(normal2surf,[length(next),1])));
            [vvv,next2]=min(concord2);
            next = next(next2);
        end
        
        here = nbs2(next(1));
    end

end
dateline(end+1) = here;

% Straighten the dateline as much as possible
dateline2=[];
curIDX = length(dateline);
while curIDX >1
    curr=dateline(curIDX);
    dateline2(end+1) = curr;
    
    %(1) Find all neighbor vertices to the current vertex
    nbs = find(AM(curr,:)==1);   
    
    ii = find(ismember(nbs,dateline(1:curIDX-1))==1);
    if length(ii) == 0
        disp('Error in dateline');
    elseif length(ii) == 1
        curIDX = curIDX-1;
    else
        [commons, ia, ib] = intersect(dateline(1:curIDX-1),nbs(ii));
        curIDX = min(ia);
    end
end
dateline2(end+1) = dateline(1);
dateline = dateline2(end:-1:1);


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Longitude from Diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up matrix A
A = B;
% Cut link to pole
for i=landmarks([1,2])
    nbs = find(AM(i,:)==1);
    if ~isempty(nbs)
        ind = sub2ind(d,nbs,nbs);  
        ind2 = sub2ind(d,nbs,i*ones(size(nbs)));  
        A(ind) = A(ind)+B(ind2);
    end
end


% walk on date line (including north/south poles)
for nd=dateline
    A(nd,:) = 0;
    A(nd,nd) = 1;
end

% Set up constant vector b
b = sparse(vertnum,1);

% walk on date line (excluding north/south poles)
for i=2:length(dateline)-1
    curr = dateline(i);
    ii1 = find(faces(:,1)==curr); ii2 = find(faces(:,2)==curr); ii3 = find(faces(:,3)==curr);
    related_faces = [faces(ii1,:);faces(ii2,[2 3 1]);faces(ii3,[3 1 2])];
    pts = [];
    while ~isempty(related_faces)
        ii = find(ismember(related_faces(:,2),[dateline(1:i-1) pts]));
        switch length(ii)
            case 0
                break;
            case 1
                vv = related_faces(ii,3);
                if ismember(vv,dateline(i+1:end))
                   break
                end
                pts(end+1) = vv;
                related_faces = [related_faces(1:(ii-1),:);related_faces((ii+1):end,:)];
            otherwise
                disp(sprintf('There are %d choices. Investigate!',length(ii))); break;
        end
    end
    pts = pts';

    ind2 = sub2ind(d,pts,dateline(i)*ones(size(pts)));  %gen_utils('con_2d_to_1d',pts,dateline(i),d);
    b(pts) = b(pts)-2*pi*B(ind2);
end

disp('Solving simultaneous Linear Equations for longitude ...');
phi = A\b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% before this point
%   theta: latitude (0,pi) (from south to north)
%   phi: longitude (0,2pi)
% after this point
%   theta: longitude (-pi,pi)
%   phi: latitude (-pi/2,pi/2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert the spherical coordinates of vertices to cartesian
% coordiantes
theta2 = phi;
phi = theta-pi/2;
theta = theta2;
ind = find(theta>pi);
theta(ind) = theta(ind)-2*pi;
clear theta2;

landmarks = locate_landmarks(theta,phi,landmarks);
landmarks = landmarks';

[sph_verts(:,1),sph_verts(:,2),sph_verts(:,3)] = sph2cart(theta,phi,1);
sph_verts = full(sph_verts);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reduce the distortion of vertices near pole ares by applying STPS 
% to latitude and longitude separately.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply STPS to latitude first

% % Calculate the length of dateline from the north pole to the south pole
total_dist = 0;
lat = zeros(length(dateline),1);
for i=2:length(dateline)
    total_dist = total_dist+WAM(dateline(i),dateline(i-1));
    lat(i) = total_dist;
end
theta_dateline = pi*(1-lat./total_dist);

% By using the neighbor vertices to two poles as template landmarks, apply
% STPS method to distort spherical mesh

%(1) Find the template landmarks 
nbs1 = find(AM(dateline(1),:)==1)';
theta_nbs1 = (pi-theta_dateline(2))*lat(2)./WAM(dateline(1),nbs1)';
theta_nbs1 = pi-theta_nbs1;

nbs2 = find(AM(dateline(end),:)==1)';
theta_nbs2 = theta_dateline(end-1)*(WAM(dateline(end),nbs2)./WAM(dateline(end),dateline(end-1)))';

nbs3 = dateline(3:end-2)';
theta_nbs3 = theta_dateline(3:end-2);

ldmk_points=[dateline(1);nbs1;nbs2;nbs3;dateline(end)];
theta_ldmk = [pi; theta_nbs1;theta_nbs2;theta_nbs3;0];

ldmk_target = sph_verts(ldmk_points,:);
[ldmks(:,1), ldmks(:,2), ldmks(:,3)] = cart2sph(ldmk_target(:,1),ldmk_target(:,2),ldmk_target(:,3));
indp2 = find(ldmks(:,1)<0); ldmks(indp2,1) = ldmks(indp2,1)+2*pi;
[ldmks_temp, sortIDX] = sort(ldmks(:,1));
ldmk_target = ldmk_target(sortIDX,:);


theta_ldmk = theta_ldmk(sortIDX);
phi_ldmk = ldmks_temp;

% Convert the spherical coordinates of template landmarks to cartesian coordiantes
theta2_ldmk = phi_ldmk;
phi_ldmk = theta_ldmk-pi/2;
theta_ldmk = theta2_ldmk;
ind = find(theta_ldmk>pi);
theta_ldmk(ind) = theta_ldmk(ind)-2*pi;
clear theta2_ldmk;

[ldmk_template(:,1),ldmk_template(:,2),ldmk_template(:,3)] = sph2cart(theta_ldmk,phi_ldmk,ones(length(theta_ldmk),1));

% Apply STPS to change latitudes
[new_ldmk_target, sph_verts] = calc_STPS(ldmk_template, ldmk_target, sph_verts);

return;
endfunction


%
% locate landmarks like center of east/west hemisphere, north/south pole
%

function landmarks = locate_landmarks(theta, phi, landmarks)

disp(sprintf('North Pole : %05d (%f, %f) pi',landmarks(1),0,-0.5));
disp(sprintf('South Pole : %05d (%f, %f) pi',landmarks(2),0,0.5));

values = phi.^2 + (theta+pi/2).^2;
[value, east] = min(values);
landmarks(3) = east; % center of east hemisphere
%disp(sprintf('East Center: %05d (%f, %f) pi',east,double(theta(east))/pi,double(phi(east))/pi));

values = phi.^2 + (theta-pi/2).^2;
[value, west] = min(values);
landmarks(4) = west; % center of west hemisphere
%disp(sprintf('West Center: %05d (%f, %f) pi',west,double(theta(west))/pi,double(phi(west))/pi));

return;
endfunction


%
% recalculate longitude
%

function verts = recal_longitude(verts,faces,edges,WAM,AM)

[phi,theta] = cart2sph(verts(:,1),verts(:,2),verts(:,3));

ix = find(phi==-pi|phi==pi|theta==-pi/2|theta==pi/2);
if ~isempty(ix)
    disp(sprintf('%d longitudes == 0, do something',length(ix))); return;
end

phi2 = phi; ix = find(phi2<0);       phi2(ix) = phi2(ix)+pi*2;
ix = find(abs(phi (edges(:,1))-phi (edges(:,2)))>pi ...
        | abs(phi2(edges(:,1))-phi2(edges(:,2)))>pi);

% Set up matrix B
vn = size(verts,1); d = [vn, vn];
B = sparse(vn,vn);
ind = find(WAM>0); B(ind) = -1./WAM(ind); % edges
B(sub2ind(d,1:vn,1:vn)) = -sum(B); % diagonal

% Set up constant vector b
b = sparse(vn,1);

vix = edges(ix,:); vix = unique(vix);
for k=1:length(vix)
    i = vix(k);
    B(i,:) = 0; B(i,i) = 1; b(i) = phi(i);
end

disp('Solving simultaneous Linear Equations for longitude ...');
phi = B\b;

[verts(:,1), verts(:,2), verts(:,3)] = sph2cart(phi,theta,1);

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


function [coeff, score, latent, tsquare] = princomp(x,econFlag)
%PRINCOMP Principal Components Analysis.
%   COEFF = PRINCOMP(X) performs principal components analysis on the N-by-P
%   data matrix X, and returns the principal component coefficients, also
%   known as loadings.  Rows of X correspond to observations, columns to
%   variables.  COEFF is a P-by-P matrix, each column containing coefficients
%   for one principal component.  The columns are in order of decreasing
%   component variance.
%
%   PRINCOMP centers X by subtracting off column means, but does not
%   rescale the columns of X.  To perform PCA with standardized variables,
%   i.e., based on correlations, use PRINCOMP(ZSCORE(X)).  To perform PCA
%   directly on a covariance or correlation matrix, use PCACOV.
%
%   [COEFF, SCORE] = PRINCOMP(X) returns the principal component scores,
%   i.e., the representation of X in the principal component space.  Rows
%   of SCORE correspond to observations, columns to components.
%
%   [COEFF, SCORE, LATENT] = PRINCOMP(X) returns the principal component
%   variances, i.e., the eigenvalues of the covariance matrix of X, in
%   LATENT.
%
%   [COEFF, SCORE, LATENT, TSQUARED] = PRINCOMP(X) returns Hotelling's
%   T-squared statistic for each observation in X.
%
%   When N <= P, SCORE(:,N:P) and LATENT(N:P) are necessarily zero, and the
%   columns of COEFF(:,N:P) define directions that are orthogonal to X.
%
%   [...] = PRINCOMP(X,'econ') returns only the elements of LATENT that are
%   not necessarily zero, i.e., when N <= P, only the first N-1, and the
%   corresponding columns of COEFF and SCORE.  This can be significantly
%   faster when P >> N.
%
%   See also BARTTEST, BIPLOT, CANONCORR, FACTORAN, PCACOV, PCARES, ROTATEFACTORS.

%   References:
%     [1] Jackson, J.E., A User's Guide to Principal Components,
%         Wiley, 1988.
%     [2] Jolliffe, I.T. Principal Component Analysis, 2nd ed.,
%         Springer, 2002.
%     [3] Krzanowski, W.J., Principles of Multivariate Analysis,
%         Oxford University Press, 1988.
%     [4] Seber, G.A.F., Multivariate Observations, Wiley, 1984.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 2.9.2.9 $  $Date: 2006/10/02 16:35:01 $

% When X has more variables than observations, the default behavior is to
% return all the pc's, even those that have zero variance.  When econFlag
% is 'econ', those will not be returned.
if nargin < 2, econFlag = 0; end

[n,p] = size(x);
if isempty(x)
    pOrZero = ~isequal(econFlag, 'econ') * p;
    coeff = zeros(p,pOrZero); coeff(1:p+1:end) = 1;
    score = zeros(n,pOrZero);
    latent = zeros(pOrZero,1);
    tsquare = zeros(n,1);
    return
end

% Center X by subtracting off column means
x0 = x - repmat(mean(x,1),n,1);
r = min(n-1,p); % max possible rank of X0

% The principal component coefficients are the eigenvectors of
% S = X0'*X0./(n-1), but computed using SVD.
[U,sigma,coeff] = svd(x0,econFlag); % put in 1/sqrt(n-1) later

if nargout < 2
    % When econFlag is 'econ', only (n-1) components should be returned.
    % See comment below.
    if (n <= p) && isequal(econFlag, 'econ')
        coeff(:,n) = [];
    end

else
    % Project X0 onto the principal component axes to get the scores.
    if n == 1 % sigma might have only 1 row
        sigma = sigma(1);
    else
        sigma = diag(sigma);
    end
    score = U .* repmat(sigma',n,1); % == x0*coeff
    sigma = sigma ./ sqrt(n-1);

    % When X has at least as many variables as observations, eigenvalues
    % n:p of S are exactly zero.
    if n <= p
        % When econFlag is 'econ', nothing corresponding to the zero
        % eigenvalues should be returned.  svd(,'econ') won't have
        % returned anything corresponding to components (n+1):p, so we
        % just have to cut off the n-th component.
        if isequal(econFlag, 'econ')
            sigma(n,:) = []; % make sure this shrinks as a column
            coeff(:,n) = [];
            score(:,n) = [];

        % Otherwise, set those eigenvalues and the corresponding scores to
        % exactly zero.  svd(,0) won't have returned columns of U
        % corresponding to components (n+1):p, need to fill those out.
        else
            sigma(n:p,1) = 0; % make sure this extends as a column
            score(:,n:p) = 0;
        end
    end

    % The variances of the pc's are the eigenvalues of S = X0'*X0./(n-1).
    latent = sigma.^2;

    % Hotelling's T-squared statistic is the sum of squares of the
    % standardized scores, i.e., Mahalanobis distances.  When X appears to
    % have column rank < r, ignore components that are orthogonal to the
    % data.
    if nargout == 4
        if n > 1
            q = sum(sigma > max(n,p).*eps(sigma(1)));
            if q < r
                warning('stats:princomp:colRankDefX', ...
                        ['Columns of X are linearly dependent to within machine precision.\n' ...
                         'Using only the first %d components to compute TSQUARED.'],q);
            end
        else
            q = 0;
        end
        tsquare = (n-1) .* sum(U(:,1:q).^2,2); % == sum((score*diag(1./sigma)).^2,2)
    end
end

endfunction


% ============================================
% do_pca.m
%
% Principle Component Analysis
%   Points are transformed by translation (zero-mean) and dimension-reduced-rotation
%
% Li Shen 
% 07/12/2002 - create

function [pca_ps, pca_b, var_amt, latent] = do_pca(points,pca_dim)

% Reduce the number of eigenvectors to be N-c (subjects minus classes).
% Then project onto the new basis, i.e. calculate the dot product of each existing
% subject vector with each of the eigenvectors.

pca_dim = min(pca_dim,size(points,2));
[pca_b, score, latent, tsquare] = princomp(points);
pca_b = pca_b(:,1:pca_dim);
pca_ps = points*pca_b;
var_amt = sum(latent(1:pca_dim))/sum(latent);
neg_ev_ind = find(latent<0);
if (~isempty(neg_ev_ind))
    disp(['Negative eigenvalues:', sprintf(' %f',latent(neg_ev_ind))]);
end 
latent = latent(1:pca_dim);

return;
endfunction


%
% move spherical vertices
%

function new_verts = move_spherical_vertices(vs, ldmks, c_longitude, d_longitude, c_colatitude, d_colatitude)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate new positions of sample points (vertices)
% Convert Cartesian coordinate into spherical coordinate
[longitude_vs, latitude_vs] = cart2sph(vs(:,1), vs(:,2), vs(:,3));

ind = find(longitude_vs<0); longitude_vs(ind) = longitude_vs(ind)+2*pi;
colatitude_vs = pi/2-latitude_vs;

R_vs = calculate_R_matrix(vs, ldmks);
R_vs = R_vs';

% Calculate displacement in longitude and colatitude directions
u_longitude_vs = c_longitude'*R_vs+d_longitude;
u_colatitude_vs = c_colatitude'*R_vs+d_colatitude;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate new spherical coordinates
new_longitude_vs = longitude_vs+u_longitude_vs';
new_colatitude_vs = colatitude_vs+u_colatitude_vs';

new_latitude_vs = pi/2-new_colatitude_vs;

% Limit the value of latitude between pi/2 and -pi/2
lg_LT=find(new_latitude_vs > (pi/2)); new_latitude_vs(lg_LT) = pi/2;
sm_LT=find(new_latitude_vs < (-1*pi/2)); new_latitude_vs(sm_LT) = -1*pi/2;

% Convert spherical coordinate into Cartesian coordinate
[x, y, z] = sph2cart(new_longitude_vs, new_latitude_vs, ones(size(new_longitude_vs,1),1));

new_verts = [x y z];

return;
endfunction


%
% Calculate c vector and a scalar d in displacement mapping function u
%

function [c_longitude, d_longitude, c_colatitude, d_colatitude] = calculate_displacement(ldmk1, ldmk2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate spherical coordinates of two sets of landmarks in [longitude,
% colatitude]
% Spherical coordinates of landmark 1
[longitude1, latitude1] = cart2sph(ldmk1(:,1), ldmk1(:,2), ldmk1(:,3)); 
indp1 = find(longitude1<0); longitude1(indp1) = longitude1(indp1)+2*pi;
colatitude1 = pi/2-latitude1;

% Spherical coordinates of landmark 2
[longitude2, latitude2] = cart2sph(ldmk2(:,1), ldmk2(:,2), ldmk2(:,3));  
indp2 = find(longitude2<0); longitude2(indp2) = longitude2(indp2)+2*pi;
colatitude2 = pi/2-latitude2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate mapping function u for longitude and colatitude
% u_longitude(p) = displacement in longitude
% u_colatitude(p) = displacement in colatitude
% u = sum(c*R)+d

% First, calculate c and d from the known points (landmarks)
Rn = calculate_R_matrix(ldmk2,ldmk2);
iRn =  inv(Rn);

T = ones(size(ldmk2,1),1);
I = diag(ones(size(ldmk2,1), 1));

z_longitude = longitude1-longitude2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT : Minimize longituidinal displacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lg_IDX = find(z_longitude >= pi);
z_longitude(lg_IDX) = z_longitude(lg_IDX)-2*pi;

sm_IDX = find(z_longitude <= -1*pi);
z_longitude(sm_IDX) = z_longitude(sm_IDX)+2*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z_colatitude = colatitude1-colatitude2;

d_longitude = (inv(T'*iRn*T)*T'*iRn)*z_longitude;
c_longitude = iRn*(z_longitude-d_longitude*T);

d_colatitude = (inv(T'*iRn*T)*T'*iRn)*z_colatitude;
c_colatitude = iRn*(z_colatitude-d_colatitude*T);

return;
endfunction

%
% calculate R matrix
%

function R = calculate_R_matrix(p1, p2)
[x1, y1] = size(p1);
[x2, y2] = size(p2);

R = zeros(x1, x2);

Z = p1*(p2');  % Inner product matrix (cos(r(p1,p2)))
W = (1-Z)/2; 

ind = find(W==0);

Q = (log(1+1./sqrt(W)).*(12*(W.^2)-4*(W))-12*(W.^(3/2))+6*(W)+1)/2;  
Q(ind) = 0.5;
R = real((Q-1/3)/(4*pi));  

return;
endfunction



function [new_ldmk2, new_verts] = calc_STPS(ldmk1, ldmk2, verts)

% Calculate c and d for calculating displacement in longitudinal and
% colatitudinal directions
[c_longitude, d_longitude, c_colatitude, d_colatitude] = calculate_displacement(ldmk1, ldmk2);

% Move spherical vertices to new locations
new_verts = move_spherical_vertices(verts, ldmk2, c_longitude, d_longitude, c_colatitude, d_colatitude);

% Calculate new landmark positions
new_ldmk2 = move_spherical_vertices(ldmk2, ldmk2, c_longitude, d_longitude, c_colatitude, d_colatitude);

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