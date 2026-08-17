import numpy as np
import pytest

import linkcell


def test_version():
    assert linkcell.__version__ == "0.3.1"


def test_periodic_image():
    xyz = np.ascontiguousarray([[0.2, 0.0, 0.0], [9.4, 0.0, 0.0]], dtype=np.float64)
    cell = np.ascontiguousarray([10.0, 10.0, 10.0], dtype=np.float64)
    nn = np.from_dlpack(linkcell.knearest(xyz, cell, 1))
    assert nn.shape == (2, 1)
    assert int(nn[0, 0]) == 1
    assert int(nn[1, 0]) == 0


def test_lattice_rows():
    xyz = np.ascontiguousarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
    cell = np.ascontiguousarray(
        [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]], dtype=np.float64
    )
    nn = np.from_dlpack(linkcell.knearest(xyz, cell, 1))
    assert int(nn[0, 0]) == 1


def test_mask_drops_source():
    xyz = np.ascontiguousarray(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float64
    )
    cell = np.ascontiguousarray([10.0, 10.0, 10.0], dtype=np.float64)
    mask = np.ascontiguousarray([1, 0, 1], dtype=np.int32)
    nn = np.from_dlpack(linkcell.knearest(xyz, cell, 1, mask=mask))
    assert int(nn[1, 0]) == -1
    assert int(nn[0, 0]) == 2
    assert int(nn[2, 0]) == 0


def test_k_zero():
    xyz = np.zeros((2, 3), dtype=np.float64)
    cell = np.array([10.0, 10.0, 10.0], dtype=np.float64)
    with pytest.raises(ValueError, match="k must"):
        linkcell.knearest(xyz, cell, 0)


def test_gpu_available_cpu_wheel():
    assert linkcell.gpu_available() is False
