import numpy as np
import pytest

import linkcell

torch = pytest.importorskip("torch")


@pytest.mark.skipif(not torch.cuda.is_available(), reason="no CUDA")
def test_torch_cuda_matches_host():
    xyz_h = np.ascontiguousarray([[0.2, 0.0, 0.0], [9.4, 0.0, 0.0]], dtype=np.float64)
    cell = np.ascontiguousarray([10.0, 10.0, 10.0], dtype=np.float64)
    host = np.from_dlpack(linkcell.knearest(xyz_h, cell, 1)).copy()
    xyz_d = torch.from_numpy(xyz_h).cuda()
    out = torch.from_dlpack(linkcell.knearest(xyz_d, cell, 1))
    assert out.device.type == "cuda"
    assert out.shape == (2, 1)
    assert int(out[0, 0].item()) == int(host[0, 0])
    assert int(out[1, 0].item()) == int(host[1, 0])


def test_gpu_available_is_bool():
    assert isinstance(linkcell.gpu_available(), bool)
