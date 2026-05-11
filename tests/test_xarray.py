import numpy as np

from fullpos.grids import octahedral_pl
from fullpos.reduced import pack_reduced, unpack_reduced


def test_reduced_pack_roundtrip_on_matching_rows() -> None:
    pl = octahedral_pl(2)
    values = np.ones(int(pl.sum()), dtype=np.float32)
    regular = unpack_reduced(values, pl, int(pl.max()))
    packed = pack_reduced(regular, pl)
    np.testing.assert_allclose(packed, values)
