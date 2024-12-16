import numpy as np
import flattensurface as flat
import pytest


@pytest.mark.parametrize("method", [
    'authalic', 'barycentric', 'conformal', 'meanvalue', 'freeconformal'
])
def test_global_parametrization(method):
    ver, tri = flat.read_off('box.off')
    uv = flat.global_parameterization(ver, tri, method=method)
    assert uv.shape == (ver.shape[0], 2)
    assert np.all(np.isfinite(uv))


def test_freerigid_parametrization():
    ver, tri = flat.read_off('torso.off')
    uv = flat.global_parameterization(ver, tri[:-1], method='freerigid')
    assert uv.shape == (ver.shape[0], 2)
    assert np.all(np.isfinite(uv))
    

def test_spherical_parameterization():
    ver, tri = flat.read_off('torso.off')
    xyz = flat.spherical_parameterization(ver, tri)
    assert xyz.shape == (ver.shape[0], 3)
    assert np.all(np.isfinite(xyz))
