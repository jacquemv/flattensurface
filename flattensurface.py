import os
import subprocess
import tempfile
import numpy as np
  
#------------------------------------------------------------------------------
def global_parameterization(vertices, triangles, method='meanvalue', 
                            rigidity=1000):
    """Triangulated surface parameterization onto a circular or free-border
    two-dimensional parameter space. The surface must be oriented and must have
    at least one boundary.
    https://doc.cgal.org/latest/Surface_mesh_parameterization

    The library CGAL must be installed.

    Args:
        vertices (nv-by-3 float array): vertex positions
        triangles (nt-by-3 int array): vertex indices for each triangle
        method (str): 
            'authalic' (discrete authalic, Desbrun et al.),
            'barycentric' (Tutte barycentric mapping),
            'conformal' (discrete conformal map, Eck et al.),
            'meanvalue' (mean value coordinates, Floater et al.)
            'freerigid' (as rigid as possible, free border, Liu et al.),
            'freeconformal' (least square conformal map, free border,
             Levy et al.)
        rigidity (float): parameter for 'freerigid' method (default: 1000);
            if 0 the method is equivalent to 'freeconformal'; larger values
            give more importance to shape preservation

    Returns:
        nv-by-2 array: parameters u (first column) and v (second column)
    """
    method = {'authalic': 'A', 'barycentric': 'B', 'conformal': 'C',
              'meanvalue': 'M', 'freerigid': 'R', 'freeconformal': 'L'}[method]
    path = os.path.abspath(__file__)
    dir_path = os.path.dirname(path)
    cmd = [dir_path + '/cgal/parameterize', method]
    if method == 'R':
        cmd += str(rigidity)
    with tempfile.NamedTemporaryFile() as file:
        write_off(file.name, vertices, triangles)
        proc = subprocess.run(cmd, stdin=open(file.name), 
                              stdout=subprocess.PIPE)
        out = proc.stdout.decode()
    if method == 'R':
        out = out[out.find('-----')+5:] # delete iteration output
    return np.fromstring(out, sep=' ').reshape((-1, 2))

#------------------------------------------------------------------------------
def local_polar_parameterization(vertices, triangles, pole, rmax=None):
    """Triangulated surface parameterization using local polar coordinates
    and geodesics

    The parameterization cannot be extended in regions where two shortest
    paths from the pole exist

    Algorithm "Discrete geodesic polar coordinates" by Melvaer & Reimers.
    The OpenMesh library is required.

    Args:
        vertices (nv-by-3 float array): vertex positions
        triangles (nt-by-3 int array): vertex indices for each triangle
        pole (int): index of the vertex representing the pole
        rmax (float): maximum distance from the pole
            (default: try to parameterize all the mesh);
            parameters at vertices with radii larger than rmax are set to
            NaN

    Returns:
        nv-by-2 array: parameters radius (first column) and angle (second
        column)
    """
    path = os.path.abspath(__file__)
    dir_path = os.path.dirname(path)
    cmd = [dir_path + '/openmesh/dgpc', '-', str(int(pole))]
    if rmax is not None:
        cmd.append(str(float(rmax)))
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, 
                            stdout=subprocess.PIPE, universal_newlines=True)
    write_obj(proc.stdin, vertices, triangles)
    out, err = proc.communicate()
    uv = np.fromstring(out, sep=' ').reshape((-1, 2))
    uv[uv[:, 0] > 1e300] = np.nan
    return uv

#------------------------------------------------------------------------------
def spherical_parameterization(vertices, triangles, niter=50):
    """Project a triangulated surface onto a sphere using the Spherical
    Harmonic Modeling and Analysis Toolkit
    https://www.med.upenn.edu/shenlab/spharm-mat.html
    Shen L, and Makedon FS. Spherical mapping for processing of 3-D closed 
    surfaces. Image and Vision Computing, 24(7): 743-761, 2006.

    Requires octave and lapack; creates temporary files _sphermap* in the
    working directory

    Args:
        vertices (nv-by-3 float array): vertex positions
        triangles (nt-by-3 int array): vertex indices for each triangle
        niter (int): number of smoothing iterations (default: 50)

    Returns:
        nv-by-3 array: cartesian coordinates of the vertices on the sphere
    """
    nv = vertices.shape[0]
    temptri = '_sphermap.off'
    write_off(temptri, vertices, triangles)
    with open('_sphermap.config', 'wt') as file:
        file.write(temptri + '\n' + str(niter))
    path = os.path.abspath(__file__)
    dir_path = os.path.dirname(path)
    cmd = ['octave', dir_path+'/sphermap/spharm_sphere_proj.m']
    subprocess.run(cmd)
    xyz = np.fromfile(temptri, sep=' ').reshape((-1, 3))[:nv]
    if os.path.exists(temptri):
        os.remove(temptri)
    if os.path.exists('_sphermap.config'):
        os.remove('_sphermap.config')
    return xyz

#------------------------------------------------------------------------------
def write_off(file, vertices, triangles):
    """Export mesh in a simplified OFF format (ASCII)

    Args:
        file (str or file): file name or file object
        vertices (nv-by-3 float array): vertex positions
        triangles (nt-by-3 int array): vertex indices for each triangle
    """
    vertices = np.array(vertices, dtype=np.float64)
    triangles = np.array(triangles, dtype=np.int32)
    closefile = False
    if isinstance(file, str):
        file = open(file, 'wt')
        closefile = True
    nv, nt = vertices.shape[0], triangles.shape[0]
    print(f'OFF\n{nv} {nt} 0', file=file)
    np.savetxt(file, vertices, delimiter=' ', fmt='%.6f')
    np.savetxt(file, np.hstack((np.full((nt, 1), 3), triangles)),
               delimiter=' ', fmt='%d')
    if closefile:
        file.close()

#------------------------------------------------------------------------------
def read_off(file):
    """Import mesh in a simplified OFF format (ASCII)
    does not implement the full specification; used only for testing

    Args:
        file (str or file): file name or file object

    Returns:
        vertices (nv-by-3 float array): vertex positions
        triangles (nt-by-3 int array): vertex indices for each triangle
    """
    closefile = False
    if isinstance(file, str):
        file = open(file, 'rt')
        closefile = True
    file.readline()
    nv, nt, _ = [int(x) for x in file.readline().split()]
    ver = np.fromfile(file, sep=' ', dtype=np.float64, count=nv*3)
    tri = np.fromfile(file, sep=' ', dtype=np.int32, count=nt*4)
    if closefile:
        file.close()
    return ver.reshape((-1, 3)), tri.reshape((-1, 4))[:, 1:]

#---------------------------------------------------------------------------
def write_obj(file, vertices, triangles):
    """Export mesh in OBJ format (ASCII)

    Args:
        file (str or file): file name or file object
        vertices (nv-by-3 float array): vertex positions
        triangles (nt-by-3 int array): vertex indices for each triangle
    """
    vertices = np.array(vertices, dtype=np.float64)
    triangles = np.array(triangles, dtype=np.int32)
    closefile = False
    if isinstance(file, str):
        file = open(file, 'wt')
        closefile = True
    np.savetxt(file, vertices, delimiter=' ', fmt='v %.6f %.6f %.6f')
    np.savetxt(file, triangles+1, delimiter=' ', fmt='f %d %d %d')
    if closefile:
        file.close()
