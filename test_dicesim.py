import math

import numpy as np
import pytest

pytest.importorskip("pybullet")
pytest.importorskip("scipy")

import dicesim as ds
import dicetray as dt


@pytest.fixture(scope="module")
def tray():
    return dt.generate(dt.Params(diameter=200, sectors=240, res=2.0))


@pytest.fixture()
def sim(tray):
    s = ds.Sim(tray, kinds=("d6", "d20"))
    yield s
    s.close()


def spec(kind, pos, quat=(0, 0, 0, 1), lin=(0, 0, 0), ang=(0, 0, 0)):
    return dict(kind=kind, pos_mm=pos, quat=quat, lin_mps=lin, ang_rps=ang)


def test_d6_rests_at_half_its_edge_above_the_floor(sim, tray):
    floor = tray.V[0, 2]
    sim.set_materials(0.4, 0.2)
    (r,) = sim.run([spec("d6", (0, 0, floor + 8 + 6))], t_max=2.0)
    assert r["settled"]
    assert r["z"] == pytest.approx(floor + 8.0, abs=0.4)     # guards the unit scaling / contact margin


@pytest.mark.parametrize("mu,slides", [(0.6, False), (0.15, True)])
def test_coulomb_threshold_on_the_cone(sim, tray, mu, slides):
    # tan(18 deg) = 0.325: a die lying flat on the cone stays for mu = 0.6 and slides for mu = 0.15
    pr = tray.prof
    r0 = 0.5 * (pr["r_b"] + pr["r_s"])
    zs = float(np.interp(r0, pr["R"][:pr["i_s"] + 1], pr["Z"][:pr["i_s"] + 1]))
    a = math.radians(tray.info["min_descent_deg"])
    n = np.array([-math.sin(a), 0.0, math.cos(a)])
    # pybullet quaternion for a rotation about +y by -a, so the die's bottom face lies on the cone
    q = sim.pb.getQuaternionFromEuler([0, -a, 0])
    pos = np.array([r0, 0.0, zs]) + n * 8.05
    sim.set_materials(mu, 0.0)
    (r,) = sim.run([spec("d6", tuple(pos), q)], t_max=0.8)
    moved = math.hypot(r["x"] - r0, r["y"])
    # a staying die still settles ~2.5 mm into contact from the hand-placed pose; a sliding one runs away
    assert (moved > 20.0) if slides else (moved < 4.0)


def test_same_throw_gives_identical_results(sim, tray):
    sp = ds.random_throw(7, ["d6", "d20"], tray.prof)
    sim.set_materials(0.4, 0.3)
    a = sim.run(sp, t_max=2.0)
    b = sim.run(sp, t_max=2.0)
    assert [(round(x["x"], 6), round(x["y"], 6)) for x in a] == [(round(x["x"], 6), round(x["y"], 6)) for x in b]


def test_throws_end_up_with_finite_positions_and_a_status(sim, tray):
    sim.set_materials(0.3, 0.3)
    rows = []
    for t in range(3):
        rows += sim.run(ds.random_throw(100 + t, ["d6", "d20"], tray.prof), t_max=3.0)
    assert all(np.isfinite([r["x"], r["y"], r["z"]]).all() for r in rows)
    assert {r["status"] for r in rows} <= {"flat", "blend", "slope", "out"}


@pytest.mark.parametrize("kind", list(ds.DICE_MM))
def test_dice_geometry_is_convex_with_positive_volume(kind):
    v = ds.dice_vertices(kind, ds.DICE_MM[kind])
    assert ds.hull_volume(v) > 100.0
    assert np.allclose(v.mean(axis=0), 0.0, atol=1e-9)


def test_d10_kites_are_planar():
    v = ds.dice_vertices("d10", 8.5)
    top, up, low = v[0], v[2::2], v[3::2]
    for k in range(5):
        pts = np.array([top, up[k], low[k], up[(k + 1) % 5]])
        n = np.cross(pts[1] - pts[0], pts[3] - pts[0])
        assert abs(np.dot(pts[2] - pts[0], n)) < 1e-9
