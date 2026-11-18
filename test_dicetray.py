import math

import numpy as np
import pytest

import dicetray as dt


def analytic_volume(prof):
    # Pappus on the closed profile polygon (axis side closes it)
    R, Z = prof["R"], prof["Z"]
    r1, z1 = np.roll(R, -1), np.roll(Z, -1)
    return abs(np.sum(math.pi * (R * R + R * r1 + r1 * r1) / 3.0 * (z1 - Z)))


@pytest.mark.parametrize("d", [150, 200, 220])
def test_smooth_mesh_is_watertight_and_matches_analytic_volume(d):
    res = dt.generate(dt.Params(diameter=d))
    vol, wt = dt.mesh_stats(res.V, res.F)
    assert wt and vol > 0
    assert vol == pytest.approx(analytic_volume(res.prof), rel=2e-3)
    assert np.hypot(res.V[:, 0], res.V[:, 1]).max() == pytest.approx(d / 2, rel=1e-6)
    assert res.V[:, 2].min() == pytest.approx(0.0, abs=1e-9)


def test_random_is_deterministic_per_seed_and_differs_between_seeds():
    p = dict(variant="random", diameter=200, sectors=360)
    a = dt.generate(dt.Params(seed=3, **p)).V
    b = dt.generate(dt.Params(seed=3, **p)).V
    c = dt.generate(dt.Params(seed=4, **p)).V
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)


@pytest.mark.parametrize("d", [150, 200, 220])
@pytest.mark.parametrize("seed", range(1, 9))
def test_random_never_creates_a_trap(d, seed):
    res = dt.generate(dt.Params(variant="random", diameter=d, seed=seed, sectors=360))
    assert res.info["watertight"]
    assert res.info["min_descent_deg"] >= 8.0
    # independent re-check on the final vertices: z must rise with r on every cone column
    prof = res.prof
    M = 360
    n_rows = prof["i_s"] - prof["i_b"] + 1
    r = np.hypot(res.V[1:-1, 0], res.V[1:-1, 1]).reshape(-1, M)[prof["i_b"] - 1:prof["i_s"]]
    z = res.V[1:-1, 2].reshape(-1, M)[prof["i_b"] - 1:prof["i_s"]]
    assert r.shape[0] == n_rows
    assert (np.diff(z, axis=0) > 0).all()


def test_check_catches_oversized_features_and_autofit_repairs_them():
    big = dict(variant="random", diameter=200, amp=6.0, sigma=4.0, seed=2, sectors=360)
    with pytest.raises(dt.DescentError):
        dt.generate(dt.Params(autofit=False, **big))
    res = dt.generate(dt.Params(autofit=True, **big))
    assert res.info["amp_scale"] < 1.0 and res.info["min_descent_deg"] >= 8.0


def test_features_stay_inside_the_slope_band():
    res = dt.generate(dt.Params(variant="random", diameter=220, seed=5, sectors=360))
    r = np.hypot(res.rand[:, 0], res.rand[:, 1])
    assert (r >= res.prof["r_b"] + 3).all() and (r <= res.prof["r_s"] - 3).all()


def test_displacement_is_zero_outside_the_cone_and_at_its_ends():
    res = dt.generate(dt.Params(variant="random", diameter=200, seed=1, sectors=360))
    pr = res.prof
    for r in (pr["r_b"], pr["r_s"], pr["r_f"], pr["R_in"]):
        x = np.array([r])
        assert dt.displacement(pr, res.rand, res.stubs, res.scale, x, np.zeros(1))[0] == 0.0


def test_stubs_only_warn_and_mesh_stays_watertight():
    res = dt.generate(dt.Params(variant="random", diameter=220, stubs=6, seed=1, sectors=360))
    assert res.info["watertight"] and res.info["n_stubs"] > 0


def test_felt_recess_is_watertight():
    res = dt.generate(dt.Params(diameter=150, felt=50, floor=2.0, sectors=360))
    vol, wt = dt.mesh_stats(res.V, res.F)
    assert wt and vol > 0
    assert res.V[:, 2].min() == pytest.approx(0.0, abs=1e-9)


def test_too_small_diameter_is_rejected():
    with pytest.raises(ValueError):
        dt.generate(dt.Params(diameter=100))


def test_stl_roundtrip_triangle_count(tmp_path):
    res = dt.generate(dt.Params(diameter=150, sectors=120))
    f = tmp_path / "t.stl"
    dt.write_stl(f, res.V, res.F)
    assert f.stat().st_size == 84 + 50 * len(res.F)
    assert int(np.frombuffer(f.read_bytes()[80:84], "<u4")[0]) == len(res.F)


def test_default_floor_is_five_layers_at_0p2mm():
    res = dt.generate(dt.Params(diameter=150, sectors=120))
    assert res.V[0, 2] == pytest.approx(1.0)        # top pole = centre of the flat
    assert res.prof["Z"][res.prof["i_b"] - 5] >= 1.0


@pytest.mark.parametrize("kw", [{}, {"zone_base_r": 45.0}, {"zone_wall_z": 26.0}, {"zone_base_r": 41.0, "zone_wall_z": 30.0}])
def test_zones_are_watertight_and_add_up_to_the_whole_body(kw):
    p = dt.Params(variant="random", diameter=200, seed=2, sectors=360, **kw)
    res = dt.generate(p)
    parts = dt.build_zones(res, p)
    total = 0.0
    for name, (V, F) in parts.items():
        vol, wt = dt.mesh_stats(V, F)
        assert wt and vol > 0, name
        total += vol
    whole, _ = dt.mesh_stats(res.V, res.F)
    assert total == pytest.approx(whole, rel=1e-9)


def test_zone_parameters_are_validated():
    p = dt.Params(diameter=200, sectors=120)
    res = dt.generate(p)
    for kw in ({"zone_base_r": 10.0}, {"zone_base_r": 90.0}, {"zone_wall_z": 5.0}, {"zone_wall_z": 100.0}):
        with pytest.raises(ValueError):
            dt.build_zones(res, dt.Params(diameter=200, sectors=120, **kw))


def test_felt_recess_needs_thicker_floor():
    with pytest.raises(ValueError):
        dt.generate(dt.Params(diameter=150, felt=50))   # floor 1.0 - depth 1.0 < 0.6


def test_3mf_structure_roundtrip(tmp_path):
    import re
    import xml.etree.ElementTree as ET
    import zipfile

    p = dt.Params(variant="random", diameter=150, seed=1, sectors=180)
    res = dt.generate(p)
    parts = dt.build_zones(res, p)
    f = tmp_path / "t.3mf"
    dt.write_3mf(f, parts, bed=256, filaments=(1, 2, 4))
    z = zipfile.ZipFile(f)
    assert z.namelist() == ["[Content_Types].xml", "_rels/.rels", "3D/3dmodel.model", "Metadata/model_settings.config"]
    ns = {"m": "http://schemas.microsoft.com/3dmanufacturing/core/2015/02"}
    root = ET.fromstring(z.read("3D/3dmodel.model"))
    objs = root.findall("m:resources/m:object", ns)
    assert len(objs) == 4
    total = 0.0
    for o in objs[:3]:
        V = np.array([[float(v.get(a)) for a in "xyz"] for v in o.findall("m:mesh/m:vertices/m:vertex", ns)])
        F = np.array([[int(t.get(a)) for a in ("v1", "v2", "v3")] for t in o.findall("m:mesh/m:triangles/m:triangle", ns)])
        vol, wt = dt.mesh_stats(V, F)
        assert wt and vol > 0
        total += vol
    assert total == pytest.approx(dt.mesh_stats(res.V, res.F)[0], rel=1e-5)
    comp = [c.get("objectid") for c in objs[3].findall("m:components/m:component", ns)]
    assert comp == ["1", "2", "3"]
    item = root.find("m:build/m:item", ns)
    assert item.get("objectid") == "4" and item.get("transform").endswith("128 128 0")
    cfg = ET.fromstring(z.read("Metadata/model_settings.config"))
    assert cfg.find("object").get("id") == "4"
    got = {pt.get("id"): {m.get("key"): m.get("value") for m in pt.findall("metadata")} for pt in cfg.iter("part")}
    assert [got[i]["name"] for i in "123"] == ["base", "slope", "wall"]
    assert [got[i]["extruder"] for i in "123"] == ["1", "2", "4"]


def test_3mf_loads_in_an_independent_reader(tmp_path):
    trimesh = pytest.importorskip("trimesh")
    p = dt.Params(diameter=150, sectors=180)
    res = dt.generate(p)
    f = tmp_path / "t.3mf"
    dt.write_3mf(f, dt.build_zones(res, p))
    sc = trimesh.load(str(f))
    meshes = list(sc.geometry.values())
    assert len(meshes) == 3 and all(m.is_watertight for m in meshes)
    assert sum(m.volume for m in meshes) == pytest.approx(dt.mesh_stats(res.V, res.F)[0], rel=1e-5)


# ---- deflectors -------------------------------------------------------------------------

@pytest.mark.parametrize("mode", ["steps", "diamonds", "both"])
def test_deflector_meshes_are_watertight_and_stay_inside_the_cone_band(mode):
    base = dt.generate(dt.Params(diameter=200, sectors=360))
    res = dt.generate(dt.Params(diameter=200, sectors=360, deflectors=mode))
    vol, wt = dt.mesh_stats(res.V, res.F)
    assert wt and vol > 0
    r = np.hypot(res.V[:, 0], res.V[:, 1])
    prof = res.prof
    far = (r <= prof["r_b"] + 1e-6) | (r >= prof["r_s"] - 1e-6)
    # identical to the plain tray outside the cone: no height step at blend or fillet
    ref = dt.generate(dt.Params(diameter=200, sectors=360))
    assert np.allclose(res.V[far], ref.V[far])
    assert res.V[:, 2].max() == pytest.approx(base.V[:, 2].max())      # wall height unchanged


def test_steps_keep_the_descent_rule_and_autofit_shrinks_tall_steps():
    ok = dt.generate(dt.Params(diameter=200, deflectors="steps"))
    assert ok.info["min_descent_deg"] >= 8.0 and ok.info["n_risers"] >= 2 and ok.info["amp_scale"] == 1.0
    tall = dt.generate(dt.Params(diameter=200, deflectors="steps", step_height=5.0))
    assert tall.info["amp_scale"] < 1.0 and tall.info["min_descent_deg"] >= 8.0
    with pytest.raises(dt.DescentError):
        dt.generate(dt.Params(diameter=200, deflectors="steps", step_height=5.0, autofit=False))


def test_bars_are_seeded_and_do_not_move_the_random_bumps():
    a = dt.generate(dt.Params(variant="random", seed=5, diameter=200, sectors=180, deflectors="diamonds"))
    b = dt.generate(dt.Params(variant="random", seed=5, diameter=200, sectors=180, deflectors="diamonds"))
    c = dt.generate(dt.Params(variant="random", seed=6, diameter=200, sectors=180, deflectors="diamonds"))
    plain = dt.generate(dt.Params(variant="random", seed=5, diameter=200, sectors=180))
    assert np.array_equal(a.V, b.V) and not np.array_equal(a.V, c.V)
    assert np.array_equal(a.rand, plain.rand)


def test_bar_rings_respect_coverage_and_the_bar_height():
    p = dt.Params(diameter=200, sectors=720, deflectors="diamonds", bar_coverage=0.4)
    res = dt.generate(p)
    plain = dt.generate(dt.Params(diameter=200, sectors=720))
    lift = res.V[:, 2] - plain.V[:, 2]
    assert lift.max() == pytest.approx(p.bar_height, abs=0.15) and lift.min() >= -1e-9
    for r in res.defl["rings"]:
        bars = [b for b in res.defl["bars"] if abs(math.hypot(b[0], b[1]) - r) < 1e-6]
        half = res.defl["bars"][0][3]
        ext = 2 * half * math.sin(math.radians(p.bar_angle)) + 2 * dt.BAR_W
        assert len(bars) * ext <= 0.4 * 2 * math.pi * r + 1e-6


def test_pocket_depth_finds_a_basin_and_ignores_a_plain_slope():
    r = np.linspace(0, 10, 11)[:, None]
    Z = np.repeat(r, 36, axis=1)                 # rises outwards, row 0 is the outlet
    assert dt.pocket_depth(Z).max() == pytest.approx(0.0)
    Z2 = Z.copy()
    Z2[4, :] += 3.0                              # a complete ring dam: water ponds uphill of it up to the crest
    assert dt.pocket_depth(Z2).max() == pytest.approx(2.0)
    Z2[4, 7] -= 3.0                              # one gap in the ring lets the water out again
    assert dt.pocket_depth(Z2).max() == pytest.approx(0.0)


def test_closed_ring_of_tangential_bars_is_repaired():
    p = dt.Params(diameter=200, sectors=360, deflectors="diamonds", bar_coverage=1.0,
                  bar_angle=90.0, bar_length=40.0)
    res = dt.generate(p)
    assert res.info["bars_removed"] >= 1 and res.info["max_pocket_mm"] <= dt.POCKET_TOL
    assert any("removed" in w for w in res.info["warnings"])


def test_unknown_deflector_mode_is_rejected():
    with pytest.raises(ValueError):
        dt.generate(dt.Params(deflectors="spikes"))


def test_zones_with_deflectors_still_sum_to_the_body():
    p = dt.Params(diameter=200, sectors=240, deflectors="both")
    res = dt.generate(p)
    parts = dt.build_zones(res, p)
    total = sum(dt.mesh_stats(*parts[n])[0] for n in parts)
    assert all(dt.mesh_stats(*parts[n])[1] for n in parts)
    assert total == pytest.approx(dt.mesh_stats(res.V, res.F)[0], rel=1e-5)
