#!/usr/bin/env python3
# === CONTEXT ===
# Purpose:     Generate a printable STL of a round dice tray whose sloped ring funnels
#              thrown dice into a flat centre. Variant "smooth" = plain cone,
#              variant "random" = same cone plus seeded Gaussian bumps/dimples that
#              deflect dice sideways without ever trapping them.
# Decisions:   - Constant-slope cone (default 18 deg) instead of a gentle 8-10 deg bowl:
#                a die stays put on an incline if mu > tan(angle). Plastic-on-plastic mu
#                is assumed ~0.3-0.5 (NOT measured), so the slope must be steeper than that.
#              - Randomisers keep their steepest flank below the base slope, so the radial
#                derivative dz/dr stays > 0 everywhere => no local minimum => no trap.
#                The rule is checked on the finished vertices (check_descent) and the
#                amplitude is scaled down automatically (autofit) instead of trusting
#                the closed-form estimate; envelopes add gradient the estimate ignores.
#              - Mesh = revolved profile on a polar grid, not a heightmap + boolean ops:
#                exact blend arcs and a vertical wall, watertight by construction, numpy only.
#              - Output is a solid body; hollowing is left to the slicer (infill), because a
#                sloped underside would need supports.
#              - STL for a single-colour body. 3MF (--format 3mf) added on request for
#                multi-colour: three touching parts (base, slope, wall) in one object, each with
#                its own extruder in Metadata/model_settings.config. Layout follows the part and
#                extruder keys read by the bbs_3mf.cpp loader of OrcaSlicer and ElegooSlicer; it
#                has not been opened in ElegooSlicer itself. Zone cuts are planes/cylinders of
#                the revolve, so each part is a watertight solid and the volumes add up exactly.
#                --zone-stls writes the same three parts as STLs for a manual multi-part import.
#              - floor = 1.0 mm (5 layers at 0.2 mm): the centre is nearly flush with the bed
#                layer, which saves material and time. A felt recess needs a thicker --floor.
# Constraints: Python >= 3.10, numpy. matplotlib only for --preview.
#              Bed of the Elegoo Centauri Carbon 2 is 256 x 256 mm -> warn above 250 mm.
# Invariants:  Profile order is top pole -> outwards -> rim -> outer wall -> bottom pole;
#              the triangle winding in build_mesh depends on it (volume must come out > 0).
#              Poles are single vertices. Displacement is applied to cone rows only and is
#              exactly 0 on the first/last cone row; displacing blend/fillet rows would
#              create height steps and catch dice.
# Related:     DESIGN-v1.md (repo root; also dice-tray/ in the 3DPrinting project) holds the design
#              rationale and print settings; dicesim.py imports generate() and the profile radii;
#              test_dicetray.py checks mesh, zone and descent rules.
# Tests:       python3 -I -m pytest test_dicetray.py
# Limitations: Dice motion is not simulated; mu and the 18 deg default need a test print.
#              --stubs (steep pegs) can break the no-trap guarantee, so the descent check
#              only warns for them. No sector/wedge export. 3MF zones: --zone-wall-z must lie on
#              the vertical wall and --zone-base-r between flat edge and blend end (undisplaced rows).
# === END CONTEXT ===
import argparse
import math
import sys
from dataclasses import dataclass, field

import numpy as np


@dataclass
class Params:
    variant: str = "smooth"          # "smooth" | "random"
    diameter: float = 200.0          # outer diameter, mm
    flat: float | None = None        # flat centre diameter; None -> 0.41 * diameter
    slope: float = 18.0              # cone angle, deg
    blend_radius: float = 30.0       # radius of the arc between flat and cone, mm
    rim: float = 4.0                 # wall thickness, mm
    wall_fillet: float = 6.0         # radius of the arc between cone and wall, mm
    wall_height: float = 12.0        # straight wall above the fillet, mm
    floor: float = 1.0               # material under the flat centre, mm (5 layers at 0.2 mm)
    felt: float = 0.0                # diameter of an optional felt recess, 0 = none
    felt_depth: float = 1.0
    chamfer: float = 0.6             # bottom edge chamfer (elephant foot), mm
    seed: int = 1
    features: int | None = None      # None -> about one per 600 mm2 of the feature band
    amp: float = 1.2                 # max bump height / dimple depth, mm
    sigma: float = 6.5               # feature width, mm
    dimple_ratio: float = 0.5
    stubs: int = 0
    min_slope: float = 8.0           # required steepest-descent angle everywhere on the cone, deg
    autofit: bool = True
    res: float = 0.75                # profile sampling step, mm
    sectors: int = 720
    zone_base_r: float | None = None  # 3MF: outer radius of the base zone; None -> end of blend
    zone_wall_z: float | None = None  # 3MF: height where the wall zone starts; None -> top of fillet


class DescentError(RuntimeError):
    pass


@dataclass
class Result:
    V: np.ndarray
    F: np.ndarray
    info: dict = field(default_factory=dict)


def smoothstep(t):
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def build_profile(p: Params):
    """Radial profile (r, z). Returns dict with arrays R, Z and key indices/radii."""
    flat = p.flat if p.flat is not None else round(0.41 * p.diameter)
    r_f, R_out = flat / 2.0, p.diameter / 2.0
    R_in = R_out - p.rim
    a = math.radians(p.slope)
    Rb, Rw = p.blend_radius, p.wall_fillet
    if p.felt and p.felt / 2.0 >= r_f:
        raise ValueError("felt recess must be smaller than the flat centre")
    if p.felt and p.floor - p.felt_depth < 0.6:
        raise ValueError("felt recess too deep for this floor: raise --floor or lower --felt-depth")

    pts = []
    if p.felt:
        rfelt = p.felt / 2.0
        pts += [(0.0, p.floor - p.felt_depth), (rfelt, p.floor - p.felt_depth), (rfelt, p.floor)]
    else:
        pts.append((0.0, p.floor))
    pts.append((r_f, p.floor))

    nb = max(3, math.ceil(Rb * a / p.res))
    for ph in np.linspace(0.0, a, nb + 1)[1:]:
        pts.append((r_f + Rb * math.sin(ph), p.floor + Rb * (1 - math.cos(ph))))
    i_b = len(pts) - 1
    r_b, z_b = pts[-1]

    r_s = R_in - Rw * (1 - math.sin(a))
    if r_s - r_b < 15.0:
        raise ValueError("slope region too short: reduce flat/blend/fillet or raise diameter")
    z_s = z_b + (r_s - r_b) * math.tan(a)
    ncone = max(2, math.ceil((r_s - r_b) / math.cos(a) / p.res))
    for t in np.linspace(0.0, 1.0, ncone + 1)[1:]:
        pts.append((r_b + t * (r_s - r_b), z_b + t * (z_s - z_b)))
    i_s = len(pts) - 1

    nf = max(3, math.ceil(Rw * (math.pi / 2 - a) / p.res))
    for ph in np.linspace(a, math.pi / 2, nf + 1)[1:]:
        pts.append((r_s + Rw * (math.sin(ph) - math.sin(a)), z_s + Rw * (math.cos(a) - math.cos(ph))))
    z_wall = pts[-1][1]
    i_w = len(pts) - 1
    z_top = z_wall + p.wall_height
    pts += [(R_in, z_top), (R_out, z_top), (R_out, p.chamfer), (R_out - p.chamfer, 0.0), (0.0, 0.0)]

    A = np.array(pts)
    return dict(R=A[:, 0], Z=A[:, 1], i_b=i_b, i_s=i_s, i_w=i_w, r_b=r_b, r_s=r_s, r_f=r_f,
                R_in=R_in, R_out=R_out, z_top=z_top, z_wall=z_wall, flat=flat)


def _dart_throw(rng, n, lo, hi, dmin, draw):
    """Seeded dart throwing in the annulus [lo, hi]; draw() returns (A, sigma)."""
    placed = []
    for _ in range(20000):
        if len(placed) >= n:
            break
        r = math.sqrt(rng.uniform(lo * lo, hi * hi))
        th = rng.uniform(0.0, 2 * math.pi)
        x, y = r * math.cos(th), r * math.sin(th)
        if all(math.hypot(x - q[0], y - q[1]) >= dmin for q in placed):
            A, s = draw()
            placed.append((x, y, A, s))
    return np.array(placed).reshape(-1, 4)


def make_features(p: Params, prof, rng):
    lo, hi = prof["r_b"] + 3.0, prof["r_s"] - 3.0
    n = p.features if p.features is not None else int(math.pi * (hi * hi - lo * lo) / 600.0)

    def draw():
        sign = -1.0 if rng.random() < p.dimple_ratio else 1.0
        return sign * p.amp * rng.uniform(0.5, 1.0), p.sigma * rng.uniform(0.8, 1.2)

    rand = _dart_throw(rng, n, lo, hi, 2.4 * p.sigma, draw)
    stubs = np.zeros((0, 4))
    if p.stubs:
        stubs = _dart_throw(rng, p.stubs, prof["r_b"] + 8.0, prof["r_s"] - 8.0, 14.0, lambda: (3.5, 3.0))
    return rand, stubs


def _gauss(X, Y, feats):
    d = np.zeros_like(X)
    for x0, y0, A, s in feats:
        d += A * np.exp(-((X - x0) ** 2 + (Y - y0) ** 2) / (2 * s * s))
    return d


def displacement(prof, rand, stubs, scale, X, Y):
    """Height offset on the cone; exactly 0 at r<=r_b and r>=r_s."""
    r = np.hypot(X, Y)
    r_b, r_s = prof["r_b"], prof["r_s"]
    ramp = 0.5 * (r_s - r_b)
    env = smoothstep((r - r_b) / ramp) * (1.0 - smoothstep((r - (r_s - 8.0)) / 8.0))
    edge = smoothstep((r - r_b) / 6.0) * (1.0 - smoothstep((r - (r_s - 6.0)) / 6.0))
    d = scale * env * _gauss(X, Y, rand) + edge * _gauss(X, Y, stubs)
    return np.where((r > r_b) & (r < r_s), d, 0.0)


def check_descent(prof, Zc):
    """Smallest radial descent angle (deg) over all cone columns; Zc rows = cone profile points."""
    Rc = prof["R"][prof["i_b"]:prof["i_s"] + 1]
    slope = np.diff(Zc, axis=0) / np.diff(Rc)[:, None]
    return math.degrees(math.atan(slope.min()))


def build_mesh(R, Z, M):
    """Revolve rows (R[k], Z[k, :]) around z. R[0] == 0 and R[-1] == 0 mark single-vertex
    poles (solid touching the axis); with no poles the polygon is closed (last row joins first).
    Rows must run clockwise in the (r, z) plane for outward-facing triangles."""
    n = len(R)
    top, bot = R[0] == 0.0, R[-1] == 0.0
    assert top == bot, "rows must have either two poles or none"
    lo, hi, off = (1, n - 1, 1) if top else (0, n, 0)
    th = 2 * np.pi * np.arange(M) / M
    V = np.empty((off + (hi - lo) * M + (1 if bot else 0), 3))
    if top:
        V[0] = (0, 0, Z[0, 0])
        V[-1] = (0, 0, Z[-1, 0])
    for k in range(lo, hi):
        s = off + (k - lo) * M
        V[s:s + M, 0] = R[k] * np.cos(th)
        V[s:s + M, 1] = R[k] * np.sin(th)
        V[s:s + M, 2] = Z[k]
    j = np.arange(M)
    j1 = (j + 1) % M
    idx = lambda k, jj: off + (k - lo) * M + jj
    F = []
    if top:
        F.append(np.stack([np.zeros(M, int), idx(lo, j), idx(lo, j1)], 1))
        F.append(np.stack([idx(hi - 1, j1), idx(hi - 1, j), np.full(M, len(V) - 1)], 1))
    pairs = [(k, k + 1) for k in range(lo, hi - 1)] + ([] if top else [(hi - 1, lo)])
    for k, k2 in pairs:
        a, b, c, d = idx(k, j), idx(k2, j), idx(k2, j1), idx(k, j1)
        F.append(np.stack([a, b, c], 1))
        F.append(np.stack([a, c, d], 1))
    return V, np.concatenate(F)


def mesh_stats(V, F):
    v0, v1, v2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    vol = float(np.einsum("ij,ij->i", v0, np.cross(v1, v2)).sum() / 6.0)
    e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]).astype(np.int64)
    N = len(V)
    key, rkey = e[:, 0] * N + e[:, 1], e[:, 1] * N + e[:, 0]
    watertight = bool(np.unique(key).size == key.size and np.array_equal(np.sort(key), np.sort(rkey)))
    return vol, watertight


def generate(p: Params) -> Result:
    prof = build_profile(p)
    R, Z0 = prof["R"], prof["Z"]
    n, M = len(R), p.sectors
    Z = np.repeat(Z0[:, None], M, axis=1)
    i_b, i_s = prof["i_b"], prof["i_s"]
    warnings, scale = [], 1.0
    rand = stubs = np.zeros((0, 4))

    th = 2 * np.pi * np.arange(M) / M
    rows = np.arange(i_b, i_s + 1)
    X, Y = R[rows, None] * np.cos(th)[None, :], R[rows, None] * np.sin(th)[None, :]

    if p.variant == "random":
        rand, stubs = make_features(p, prof, np.random.default_rng(p.seed))
        if p.features is not None and len(rand) < p.features:
            warnings.append(f"only {len(rand)} of {p.features} features fit")
        while True:
            D = displacement(prof, rand, stubs, scale, X, Y)
            D[0], D[-1] = 0.0, 0.0
            ms = check_descent(prof, Z0[rows, None] + D)
            if ms >= p.min_slope or not p.autofit or scale < 0.05:
                break
            scale *= 0.9
        Z[rows] += D
    else:
        ms = check_descent(prof, Z0[rows, None] + np.zeros_like(X))

    if ms < p.min_slope:
        msg = f"descent check failed: steepest-descent angle {ms:.1f} deg < {p.min_slope:.1f} deg"
        if len(stubs):
            warnings.append(msg + " (stubs present, trap-free guarantee void)")
        else:
            raise DescentError(msg)

    V, F = build_mesh(R, Z, M)
    vol, wt = mesh_stats(V, F)
    info = dict(diameter=2 * prof["R_out"], flat=prof["flat"], height=float(V[:, 2].max()),
                slope_run=prof["r_s"] - prof["r_b"], wall_base_z=prof["z_wall"], volume_cm3=vol / 1000.0,
                triangles=len(F), watertight=wt, min_descent_deg=ms, amp_scale=scale,
                n_features=len(rand), n_stubs=len(stubs), warnings=warnings)
    if info["diameter"] > 250:
        warnings.append("diameter > 250 mm: tight on a 256 mm bed")
    res = Result(V, F, info)
    res.prof, res.rand, res.stubs, res.scale, res.Z = prof, rand, stubs, scale, Z
    return res


def write_stl(path, V, F):
    t = V[F]
    n = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    n = n / np.where(ln == 0, 1, ln)
    rec = np.zeros(len(F), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    rec["n"], rec["v"] = n, t
    with open(path, "wb") as f:
        f.write(b"dicetray".ljust(80, b" "))
        f.write(np.uint32(len(F)).tobytes())
        f.write(rec.tobytes())


def build_zones(res: Result, p: Params):
    """Split the body into three watertight solids for multi-colour printing.
    base  = r <= zone_base_r (flat centre + blend), full height down to z=0
    slope = the ring outside it, up to zone_wall_z
    wall  = everything above zone_wall_z (upper wall and rim)
    Returns {name: (V, F)}; the volumes add up to the volume of the whole body."""
    prof, M, Zf = res.prof, p.sectors, res.Z
    R, i_b, i_w = prof["R"], prof["i_b"], prof["i_w"]
    r_base = prof["r_b"] if p.zone_base_r is None else p.zone_base_r
    z_c = prof["z_wall"] if p.zone_wall_z is None else p.zone_wall_z
    if not (prof["r_f"] - 1e-9 <= r_base <= prof["r_b"] + 1e-9):
        raise ValueError(f"zone base radius must be within {prof['r_f']:.1f}..{prof['r_b']:.1f} mm")
    if not (prof["z_wall"] - 1e-9 <= z_c <= prof["z_top"] - 1.0):
        raise ValueError(f"zone wall height must be within {prof['z_wall']:.1f}..{prof['z_top'] - 1:.1f} mm")
    r_base = min(r_base, prof["r_b"])
    R_in, R_out, z_top = prof["R_in"], prof["R_out"], prof["z_top"]
    const = lambda z: np.full(M, float(z))

    k0 = max(k for k in range(i_b + 1) if R[k] <= r_base + 1e-9)
    base = [(R[k], Zf[k]) for k in range(k0 + 1)]
    if R[k0] < r_base - 1e-9:        # cut falls inside the blend arc: add a vertex on the chord
        t = (r_base - R[k0]) / (R[k0 + 1] - R[k0])
        x = (r_base, Zf[k0] + t * (Zf[k0 + 1] - Zf[k0]))
        base.append(x)
        slope, start = [x], k0 + 1
    else:
        slope, start = [], k0
    slope += [(R[k], Zf[k]) for k in range(start, i_w + 1)]
    if z_c > prof["z_wall"] + 1e-9:
        slope.append((R_in, const(z_c)))
    slope += [(R_out, const(z_c)), (R_out, const(p.chamfer)), (R_out - p.chamfer, const(0)), (r_base, const(0))]
    base += [(r_base, const(0)), (0.0, const(0))]
    wall = [(R_in, const(z_c)), (R_in, const(z_top)), (R_out, const(z_top)), (R_out, const(z_c))]

    def mesh(rows):
        r = np.array([a for a, _ in rows], dtype=float)
        z = np.stack([np.broadcast_to(np.asarray(b, dtype=float), (M,)) for _, b in rows])
        return build_mesh(r, z, M)

    return {"base": mesh(base), "slope": mesh(slope), "wall": mesh(wall)}


def write_3mf(path, parts, title="dice tray", bed=256.0, filaments=(1, 2, 3)):
    """One 3MF object made of touching parts, each assigned to its own filament.
    Layout: 3MF core model (parts as objects + one component object) plus
    Metadata/model_settings.config with per-part name/extruder, as read by OrcaSlicer-family slicers."""
    import zipfile

    def mesh_xml(V, F):
        vs = "\n".join(f'     <vertex x="{x:.5f}" y="{y:.5f}" z="{z:.5f}"/>' for x, y, z in V.tolist())
        ts = "\n".join(f'     <triangle v1="{a}" v2="{b}" v3="{c}"/>' for a, b, c in F.tolist())
        return f"    <mesh>\n     <vertices>\n{vs}\n     </vertices>\n     <triangles>\n{ts}\n     </triangles>\n    </mesh>"

    names = list(parts)
    obj_id = len(names) + 1
    objs = "\n".join(f'   <object id="{i}" type="model">\n{mesh_xml(*parts[n])}\n   </object>'
                     for i, n in enumerate(names, 1))
    comps = "\n".join(f'     <component objectid="{i}"/>' for i in range(1, len(names) + 1))
    model = ('<?xml version="1.0" encoding="UTF-8"?>\n'
             '<model unit="millimeter" xml:lang="en-US" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02">\n'
             f' <metadata name="Title">{title}</metadata>\n <metadata name="Application">dicetray.py</metadata>\n'
             f' <resources>\n{objs}\n   <object id="{obj_id}" type="model">\n    <components>\n{comps}\n    </components>\n   </object>\n </resources>\n'
             f' <build>\n  <item objectid="{obj_id}" transform="1 0 0 0 1 0 0 0 1 {bed / 2:g} {bed / 2:g} 0" printable="1"/>\n </build>\n</model>\n')
    ident = "1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1"
    part_cfg = "\n".join(
        f'    <part id="{i}" subtype="normal_part">\n'
        f'      <metadata key="name" value="{n}"/>\n'
        f'      <metadata key="matrix" value="{ident}"/>\n'
        f'      <metadata key="extruder" value="{filaments[i - 1]}"/>\n'
        f'    </part>' for i, n in enumerate(names, 1))
    cfg = ('<?xml version="1.0" encoding="UTF-8"?>\n<config>\n'
           f'  <object id="{obj_id}">\n    <metadata key="name" value="{title}"/>\n'
           f'    <metadata key="extruder" value="{filaments[0]}"/>\n{part_cfg}\n  </object>\n</config>\n')
    types = ('<?xml version="1.0" encoding="UTF-8"?>\n'
             '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">\n'
             ' <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>\n'
             ' <Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/>\n'
             ' <Default Extension="config" ContentType="text/xml"/>\n</Types>\n')
    rels = ('<?xml version="1.0" encoding="UTF-8"?>\n'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">\n'
            ' <Relationship Target="/3D/3dmodel.model" Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/>\n'
            '</Relationships>\n')
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", types)
        z.writestr("_rels/.rels", rels)
        z.writestr("3D/3dmodel.model", model)
        z.writestr("Metadata/model_settings.config", cfg)


def write_preview(path, res: Result, p: Params):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource

    prof = res.prof
    g = np.linspace(-prof["R_in"], prof["R_in"], 600)
    X, Y = np.meshgrid(g, g)
    r = np.hypot(X, Y)
    up = prof["i_s"] + 1
    z = np.interp(r, prof["R"][:up], prof["Z"][:up]) + displacement(prof, res.rand, res.stubs, res.scale, X, Y)
    z = np.where(r <= prof["r_s"], z, np.nan)
    fig, ax = plt.subplots(1, 2, figsize=(11, 5.2))
    rgb = LightSource(azdeg=315, altdeg=40).shade(np.nan_to_num(z), plt.cm.cividis, vert_exag=6, blend_mode="overlay")
    rgb[np.isnan(z)] = 1.0
    ax[0].imshow(rgb, extent=[g[0], g[-1], g[0], g[-1]], origin="lower")
    ax[0].set_title(f"{p.variant}, d={p.diameter:g} mm, seed {p.seed} (hillshade, z exaggerated)")
    d = displacement(prof, res.rand, res.stubs, res.scale, X, Y)
    im = ax[1].imshow(np.where(r <= prof["r_s"], d, np.nan), extent=[g[0], g[-1], g[0], g[-1]],
                      origin="lower", cmap="RdBu_r", vmin=-p.amp, vmax=p.amp)
    ax[1].set_title("height offset vs. plain cone [mm]")
    fig.colorbar(im, ax=ax[1], shrink=0.8)
    for a_ in ax:
        a_.set_xlabel("x [mm]")
    fig.tight_layout()
    fig.savefig(path, dpi=110)


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Round dice tray generator (STL / 3MF).")
    d = Params()
    ap.add_argument("--variant", choices=["smooth", "random"], default=d.variant)
    ap.add_argument("--diameter", type=float, default=d.diameter)
    ap.add_argument("--flat", type=float, default=None, help="flat centre diameter (default 0.41 x diameter)")
    ap.add_argument("--slope", type=float, default=d.slope)
    ap.add_argument("--blend-radius", type=float, default=d.blend_radius)
    ap.add_argument("--rim", type=float, default=d.rim)
    ap.add_argument("--wall-fillet", type=float, default=d.wall_fillet)
    ap.add_argument("--wall-height", type=float, default=d.wall_height)
    ap.add_argument("--floor", type=float, default=d.floor)
    ap.add_argument("--felt", type=float, default=d.felt, help="felt recess diameter, 0 = off")
    ap.add_argument("--felt-depth", type=float, default=d.felt_depth)
    ap.add_argument("--seed", type=int, default=d.seed)
    ap.add_argument("--features", type=int, default=None)
    ap.add_argument("--amp", type=float, default=d.amp)
    ap.add_argument("--sigma", type=float, default=d.sigma)
    ap.add_argument("--dimple-ratio", type=float, default=d.dimple_ratio)
    ap.add_argument("--stubs", type=int, default=d.stubs)
    ap.add_argument("--min-slope", type=float, default=d.min_slope)
    ap.add_argument("--no-autofit", action="store_true")
    ap.add_argument("--res", type=float, default=d.res)
    ap.add_argument("--sectors", type=int, default=d.sectors)
    ap.add_argument("--format", choices=["stl", "3mf", "both"], default="stl",
                    help="3mf = three colour parts (base, slope, wall) in one object")
    ap.add_argument("--zone-base-r", type=float, default=None, help="base zone outer radius, mm")
    ap.add_argument("--zone-wall-z", type=float, default=None, help="wall zone start height, mm")
    ap.add_argument("--filaments", type=lambda s: tuple(int(x) for x in s.split(",")), default=(1, 2, 3),
                    help="filament numbers for base,slope,wall (default 1,2,3)")
    ap.add_argument("--bed", type=float, default=256.0, help="bed size for centring the 3MF, mm")
    ap.add_argument("--zone-stls", action="store_true", help="also write the three zones as separate STLs")
    ap.add_argument("-o", "--output", help="output name; the .stl/.3mf suffix is set by --format")
    ap.add_argument("--preview", help="write a PNG preview (needs matplotlib)")
    return ap.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    p = Params(variant=a.variant, diameter=a.diameter, flat=a.flat, slope=a.slope,
               blend_radius=a.blend_radius, rim=a.rim, wall_fillet=a.wall_fillet,
               wall_height=a.wall_height, floor=a.floor, felt=a.felt, felt_depth=a.felt_depth,
               seed=a.seed, features=a.features, amp=a.amp, sigma=a.sigma,
               dimple_ratio=a.dimple_ratio, stubs=a.stubs, min_slope=a.min_slope,
               autofit=not a.no_autofit, res=a.res, sectors=a.sectors,
               zone_base_r=a.zone_base_r, zone_wall_z=a.zone_wall_z)
    if len(a.filaments) != 3:
        print("error: --filaments needs three numbers, e.g. 1,2,3", file=sys.stderr)
        return 2
    try:
        res = generate(p)
        parts = build_zones(res, p) if (a.format != "stl" or a.zone_stls) else None
    except (ValueError, DescentError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    stem = a.output or f"dicetray_{p.variant}_d{p.diameter:g}" + (f"_s{p.seed}" if p.variant == "random" else "")
    if stem.lower().endswith((".stl", ".3mf")):
        stem = stem[:-4]
    i = res.info
    if a.format in ("stl", "both"):
        write_stl(stem + ".stl", res.V, res.F)
        print(f"wrote {stem}.stl")
    if a.format in ("3mf", "both"):
        write_3mf(stem + ".3mf", parts, bed=a.bed, filaments=a.filaments)
        print(f"wrote {stem}.3mf (parts base/slope/wall -> filaments {','.join(map(str, a.filaments))}, "
              + ", ".join(f"{n} {mesh_stats(*parts[n])[0] / 1000:.0f} cm3" for n in parts) + ")")
    if a.zone_stls:
        for n, (V, F) in parts.items():
            write_stl(f"{stem}_{n}.stl", V, F)
            print(f"wrote {stem}_{n}.stl")
    print(f"  diameter {i['diameter']:.0f} mm, flat centre {i['flat']:.0f} mm, height {i['height']:.1f} mm, "
          f"slope run {i['slope_run']:.1f} mm")
    print(f"  solid volume {i['volume_cm3']:.0f} cm3, {i['triangles']} triangles, watertight {i['watertight']}")
    print(f"  min descent angle {i['min_descent_deg']:.1f} deg (required {p.min_slope:g}), "
          f"amplitude scale {i['amp_scale']:.2f}, features {i['n_features']}, stubs {i['n_stubs']}")
    for w in i["warnings"]:
        print(f"  warning: {w}", file=sys.stderr)
    if a.preview:
        write_preview(a.preview, res, p)
        print(f"wrote {a.preview}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
