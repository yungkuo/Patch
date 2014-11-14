"""
Tests for bsearch modules.

To run tests install pytest and type:

    py.tests

"""

import pytest
import numpy as np

from bsearch_py import burstsearch_py
from bsearch_c import burstsearch_c


## Burst search parameters
bsearch_params = dict(m=10, threshold=0.1)

def get_data():
    np.random.seed(1)
    signal = np.random.randn(10000)
    return signal

@pytest.fixture(scope="module")
def data(request):
    return get_data()


## Tests functions
def test_bsearch(data):
    """Smoke test bsearch_py
    """
    bursts_py = burstsearch_py(data, **bsearch_params)
    bursts_c = burstsearch_c(data, **bsearch_params)
    assert np.allclose(bursts_c.start, bursts_py.start)
    assert np.allclose(bursts_c.stop, bursts_py.stop)
    assert np.allclose(bursts_c.score, bursts_py.score)


if __name__ == '__main__':
    pytest.main("-x -v search/test_bsearch.py")
