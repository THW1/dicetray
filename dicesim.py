#!/usr/bin/env python3
# === CONTEXT ===
# Purpose:     Rigid-body simulation (pybullet) of polyhedral dice thrown into the trays made by
#              dicetray.py: how many dice end up on the flat centre, on the slope, or leave the
#              tray, per tray variant and per friction value.
# Decisions:   - The friction coefficient of dice on a printed tray is unknown, so --mu takes a
#                list and every result is a sweep, not a single number. Bullet multiplies the
#                two contact coefficients, so each body gets sqrt(mu) and mu is the effective value.
#              - Throws are aimed into the tray (landing point at rim height drawn within --aim x inner
#                radius); an unaimed random velocity mostly flung dice over the wall, which measured
#                the throw model instead of the tray (wall height made no difference).
#              - The same seeded throws are replayed for every variant and mu, so differences
#                come from the tray, not from luck.
#              - Lengths are scaled to 0.1 sim unit per mm (gravity scaled with them, so time is
#                real time). Bullet's contact margins are tuned for metre-sized bodies; at 1 unit per
#                mm they would be as large as the dice. The rest-height test guards this.
#              - The tray enters as the dicetray.py mesh (static concave triangle mesh), at a
#                coarser sampling (--sim-res) than the printable one; the surface is the same
#                because features depend only on the seed, not on the sampling.
#              - Dice are convex hulls of ideal polyhedra, no rounded edges, no pips, mass from
#                an acrylic density; inertia is that of a sphere of equal volume.
#              - Result categories use the centre of mass: flat = r <= flat radius, blend = up to the
#                end of the blend arc, slope = on the cone, out = beyond the inner wall (rim or outside).
#              - Tumbling is measured by the up-face: every TRACK_EVERY steps (12.5 ms) the face whose
#                normal points most upwards is recorded; a change between two samples in which the die
#                touches the cone (r > flat edge, within a die size of the surface) counts as a tumble.
#                "slid" = ended on the flat after touching the slope without a single change, i.e. the
#                unwanted glide that deflectors are meant to prevent. Fast double flips between two
#                samples are undercounted.
# Constraints: pybullet, numpy, scipy (hull volume); matplotlib only for --plot. Single thread,
#              deterministic for a given pybullet build.
# Invariants:  Timestep 1/960 s: a falling die moves <3 mm per step, less than half the smallest
#              die, so it cannot tunnel through the single-layer triangle mesh.
#              Linear/angular damping must stay 0; pybullet's default 0.04 would brake the dice.
# Related:     dicetray.py::generate supplies the mesh and the radii (profile dict).
# Tests:       python3 -I -m pytest test_dicesim.py
# Limitations: Dice sizes, density, restitution, throw speeds are illustrative assumptions.
#              Bullet's friction model is Coulomb with a pyramid cone; real dice on PETG/PLA
#              were not measured. The numbers rank trays and show trends; they are not
#              predictions for a specific dice set.
# === END CONTEXT ===
import argparse
import math
import sys

import numpy as np

import dicetray as dt

SCALE = 0.1                       # sim units per mm
G = 9810.0 * SCALE                # gravity in sim units/s^2
DT = 1.0 / 960.0
SOLVER_ITERS = 60
CHECK_EVERY = 48                  # steps between rest checks (0.05 s)
TRACK_EVERY = 12                  # steps between orientation samples (12.5 ms); divides CHECK_EVERY
PHI = (1 + math.sqrt(5)) / 2

# typical edge/ring dimensions in mm (illustrative)
DICE_MM = {"d4": 19.0, "d6": 16.0, "d8": 10.5, "d10": 8.5, "d12": 8.5, "d20": 13.5}


def dice_vertices(kind, size, scale=1.0):
    s = size * scale
    if kind == "d4":
        a = s / (2 * math.sqrt(2))
        v = [(a, a, a), (a, -a, -a), (-a, a, -a), (-a, -a, a)]
    elif kind == "d6":
        h = s / 2
        v = [(x, y, z) for x in (-h, h) for y in (-h, h) for z in (-h, h)]
    elif kind == "d8":
        a = s / math.sqrt(2)
        v = [(a, 0, 0), (-a, 0, 0), (0, a, 0), (0, -a, 0), (0, 0, a), (0, 0, -a)]
    elif kind == "d10":          # pentagonal trapezohedron: ring offset c = h*tan(18 deg)^2 makes the kites planar
        r = s
        h = r
        c = h * math.tan(math.radians(18)) ** 2
        v = [(0, 0, h), (0, 0, -h)]
        for k in range(5):
            v.append((r * math.cos(math.radians(72 * k)), r * math.sin(math.radians(72 * k)), c))
            v.append((r * math.cos(math.radians(72 * k + 36)), r * math.sin(math.radians(72 * k + 36)), -c))
    elif kind == "d12":
        k = s / (2 / PHI)
        v = [(x, y, z) for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
        for a in (-1, 1):
            for b in (-1, 1):
                v += [(0, a / PHI, b * PHI), (a / PHI, b * PHI, 0), (a * PHI, 0, b / PHI)]
        v = [(x * k, y * k, z * k) for x, y, z in v]
    elif kind == "d20":
        k = s / 2
        v = []
        for a in (-1, 1):
            for b in (-1, 1):
                v += [(0, a, b * PHI), (a, b * PHI, 0), (a * PHI, 0, b)]
        v = [(x * k, y * k, z * k) for x, y, z in v]
    else:
        raise ValueError(f"unknown die {kind}")
    return np.array(v, dtype=float)


def face_normals(v):
    """Unit outward normals of the faces of the convex hull of v (coplanar triangles merged)."""
    from scipy.spatial import ConvexHull
    n = ConvexHull(v).equations[:, :3]
    uniq = {}
    for row in n:
        uniq.setdefault(tuple(np.round(row, 3)), row)
    return np.array(list(uniq.values()))


def up_face(normals, quat):
    """Index of the face pointing most upwards for orientation quat (x, y, z, w)."""
    x, y, z, w = quat
    # third row of the rotation matrix: z component of R @ n
    rz = np.array([2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)])
    return int(np.argmax(normals @ rz))


def hull_volume(v):
    from scipy.spatial import ConvexHull
    return ConvexHull(v).volume


class Sim:
    def __init__(self, tray: dt.Result, dice_scale=1.0, kinds=tuple(DICE_MM)):
        import pybullet as pb
        self.pb = pb
        self.cid = pb.connect(pb.DIRECT)
        c = self.cid
        pb.setGravity(0, 0, -G, physicsClientId=c)
        pb.setPhysicsEngineParameter(fixedTimeStep=DT, numSolverIterations=SOLVER_ITERS, physicsClientId=c)
        self.prof = tray.prof
        V = (tray.V * SCALE).tolist()
        shape = pb.createCollisionShape(pb.GEOM_MESH, vertices=V, indices=tray.F.ravel().tolist(),
                                        flags=pb.GEOM_FORCE_CONCAVE_TRIMESH, physicsClientId=c)
        self.tray = pb.createMultiBody(0, shape, physicsClientId=c)
        plane = pb.createCollisionShape(pb.GEOM_PLANE, physicsClientId=c)
        self.ground = pb.createMultiBody(0, plane, physicsClientId=c)
        self.shapes, self.mass, self.inertia, self.normals = {}, {}, {}, {}
        for k in kinds:
            v = dice_vertices(k, DICE_MM[k], dice_scale)
            self.shapes[k] = pb.createCollisionShape(pb.GEOM_MESH, vertices=(v * SCALE).tolist(), physicsClientId=c)
            self.normals[k] = face_normals(v)
            vol = hull_volume(v * SCALE)                       # sim units^3
            self.mass[k] = vol * 1.2                           # acrylic ~1.2 g/cm3, arbitrary mass unit
            r_eq = (3 * vol / (4 * math.pi)) ** (1 / 3)
            self.inertia[k] = 0.4 * self.mass[k] * r_eq ** 2
        self.set_materials(0.4, 0.3)

    def set_materials(self, mu, restitution):
        self.mu, self.e = mu, restitution
        for b in (self.tray, self.ground):
            self._material(b)

    def _material(self, body, inertia=None):
        kw = dict(lateralFriction=math.sqrt(self.mu), restitution=math.sqrt(self.e), rollingFriction=0.0,
                  spinningFriction=0.0, linearDamping=0.0, angularDamping=0.0, physicsClientId=self.cid)
        if inertia is not None:
            kw["localInertiaDiagonal"] = [inertia] * 3
        self.pb.changeDynamics(body, -1, **kw)

    def run(self, specs, t_max=4.0):
        """specs: list of dict(kind, pos_mm, quat, lin_mps (3,), ang_rps (3,)).
        Returns per die dict(kind, x, y, z [mm], status, settled)."""
        pb, c = self.pb, self.cid
        bodies = []
        for s in specs:
            b = pb.createMultiBody(self.mass[s["kind"]], self.shapes[s["kind"]],
                                   basePosition=(np.asarray(s["pos_mm"]) * SCALE).tolist(),
                                   baseOrientation=list(s["quat"]), physicsClientId=c)
            self._material(b, self.inertia[s["kind"]])
            pb.resetBaseVelocity(b, (np.asarray(s["lin_mps"]) * 1000.0 * SCALE).tolist(),
                                 list(s["ang_rps"]), physicsClientId=c)
            bodies.append(b)
        prof = self.prof
        calm, settled, r_in = 0, False, prof["R_in"]
        nb = len(bodies)
        face = [None] * nb                # last up-face while touching the cone
        tumbles, touches = [0] * nb, [0] * nb
        zs_R, zs_Z = prof["R"][:prof["i_w"] + 1], prof["Z"][:prof["i_w"] + 1]
        for step in range(1, int(t_max / DT) + 1):
            pb.stepSimulation(physicsClientId=c)
            if step % TRACK_EVERY:
                continue
            check = step % CHECK_EVERY == 0
            moving = False
            for i, b in enumerate(bodies):
                p_, q_ = pb.getBasePositionAndOrientation(b, physicsClientId=c)
                pos = np.array(p_) / SCALE
                r = math.hypot(pos[0], pos[1])
                kind = specs[i]["kind"]
                above = pos[2] - float(np.interp(r, zs_R, zs_Z))
                if prof["r_f"] < r <= r_in and above <= DICE_MM[kind]:
                    f = up_face(self.normals[kind], q_)
                    touches[i] += 1
                    if face[i] is not None and f != face[i]:
                        tumbles[i] += 1
                    face[i] = f
                else:
                    face[i] = None
                if check:
                    lin, ang = pb.getBaseVelocity(b, physicsClientId=c)
                    on_table = r > prof["R_out"] + 1.0
                    if not on_table and (np.linalg.norm(lin) / SCALE > 10.0 or np.linalg.norm(ang) > 1.0):
                        moving = True
            if check:
                calm = 0 if moving else calm + 1
                if calm >= 5:
                    settled = True
                    break
        out = []
        for i, (s, b) in enumerate(zip(specs, bodies)):
            pos = np.array(pb.getBasePositionAndOrientation(b, physicsClientId=c)[0]) / SCALE
            r = math.hypot(pos[0], pos[1])
            status = ("out" if r > r_in else "flat" if r <= prof["r_f"]
                      else "blend" if r <= prof["r_b"] else "slope")
            out.append(dict(kind=s["kind"], x=pos[0], y=pos[1], z=pos[2], r=r, status=status, settled=settled,
                            tumbles=tumbles[i], touches=touches[i],
                            slid=bool(status == "flat" and touches[i] >= 3 and tumbles[i] == 0)))
            pb.removeBody(b, physicsClientId=c)
        return out

    def close(self):
        self.pb.disconnect(self.cid)


def random_throw(seed, kinds, prof, aim=0.8, drop=(20.0, 60.0)):
    """Dice start staggered in height above the tray centre (first die `drop` mm above the rim,
    24 mm more per further die), with random orientation and spin. The horizontal velocity is
    chosen so the die crosses rim height over a random point within `aim` x inner radius, i.e.
    the throw is aimed into the tray; leaving it afterwards takes a bounce."""
    rng = np.random.default_rng(seed)
    specs = []
    for i, k in enumerate(kinds):
        r = 0.35 * prof["R_in"] * math.sqrt(rng.uniform())
        th = rng.uniform(0, 2 * math.pi)
        q = rng.normal(size=4)
        q /= np.linalg.norm(q)
        h = rng.uniform(*drop) + 24 * i                            # height above the rim, mm
        rt = aim * prof["R_in"] * math.sqrt(rng.uniform())
        tt = rng.uniform(0, 2 * math.pi)
        t_fall = math.sqrt(2 * h / 9810.0)                         # s, from release to rim height
        vx = (rt * math.cos(tt) - r * math.cos(th)) / t_fall / 1000.0   # m/s
        vy = (rt * math.sin(tt) - r * math.sin(th)) / t_fall / 1000.0
        specs.append(dict(kind=k,
                          pos_mm=(r * math.cos(th), r * math.sin(th), prof["z_top"] + h),
                          quat=q, lin_mps=(vx, vy, 0.0),
                          ang_rps=tuple(rng.uniform(-25, 25, size=3))))
    return specs


def summarize(rows):
    n = len(rows)
    c = {s: sum(r["status"] == s for r in rows) for s in ("flat", "blend", "slope", "out")}
    uns = sum(not r["settled"] for r in rows)
    return n, c, uns


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Dice-throw simulation for dicetray.py trays.")
    ap.add_argument("--variants", default="smooth,random")
    ap.add_argument("--deflectors", default="none", help="comma separated list of none,diamonds,steps,both; every variant is run with each")
    ap.add_argument("--deflector-pitch", type=float, default=dt.Params.deflector_pitch)
    ap.add_argument("--step-height", type=float, default=dt.Params.step_height)
    ap.add_argument("--bar-length", type=float, default=dt.Params.bar_length)
    ap.add_argument("--bar-height", type=float, default=dt.Params.bar_height)
    ap.add_argument("--bar-angle", type=float, default=dt.Params.bar_angle)
    ap.add_argument("--bar-coverage", type=float, default=dt.Params.bar_coverage)
    ap.add_argument("--diameter", type=float, default=150.0)
    ap.add_argument("--flat", type=float, default=None)
    ap.add_argument("--slope", type=float, default=18.0)
    ap.add_argument("--floor", type=float, default=1.0)
    ap.add_argument("--blend-radius", type=float, default=30.0)
    ap.add_argument("--wall-height", type=float, default=12.0)
    ap.add_argument("--seed", type=int, default=1, help="tray feature seed (random variant)")
    ap.add_argument("--amp", type=float, default=1.2)
    ap.add_argument("--sigma", type=float, default=6.5)
    ap.add_argument("--features", type=int, default=None)
    ap.add_argument("--throws", type=int, default=20)
    ap.add_argument("--throw-seed", type=int, default=100)
    ap.add_argument("--dice", default=",".join(DICE_MM))
    ap.add_argument("--dice-scale", type=float, default=1.0)
    ap.add_argument("--mu", default="0.2,0.4,0.6", help="effective friction coefficients, comma separated")
    ap.add_argument("--restitution", type=float, default=0.3)
    ap.add_argument("--aim", type=float, default=0.8, help="dice are aimed at random points within this fraction of the inner radius")
    ap.add_argument("--drop", default="20,60", help="release height range of the lowest die above the rim, mm")
    ap.add_argument("--t-max", type=float, default=3.0)
    ap.add_argument("--sim-res", type=float, default=2.5, help="profile step of the simulated tray, mm (use ~1 with deflectors)")
    ap.add_argument("--sim-sectors", type=int, default=240, help="angular sectors of the simulated tray (use ~480 with deflectors)")
    ap.add_argument("--by-type", action="store_true")
    ap.add_argument("--plot", help="PNG with final resting positions")
    return ap.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    kinds = [k.strip() for k in a.dice.split(",")]
    mus = [float(x) for x in a.mu.split(",")]
    drop = tuple(float(x) for x in a.drop.split(","))
    results = {}
    trays = {}
    configs = [(v.strip(), d.strip()) for v in a.variants.split(",") for d in a.deflectors.split(",")]
    for variant, defl in configs:
        label = variant if defl == "none" else f"{variant}+{defl}"
        p = dt.Params(variant=variant, diameter=a.diameter, flat=a.flat, slope=a.slope, floor=a.floor,
                      wall_height=a.wall_height, blend_radius=a.blend_radius, seed=a.seed, amp=a.amp, sigma=a.sigma,
                      features=a.features, res=a.sim_res, sectors=a.sim_sectors, deflectors=defl,
                      deflector_pitch=a.deflector_pitch, step_height=a.step_height, bar_length=a.bar_length,
                      bar_height=a.bar_height, bar_angle=a.bar_angle, bar_coverage=a.bar_coverage)
        try:
            tray = dt.generate(p)
        except (ValueError, dt.DescentError) as e:
            print(f"error: {e}", file=sys.stderr)
            return 2
        sim = Sim(tray, a.dice_scale, kinds)
        trays[label] = tray
        for mu in mus:
            sim.set_materials(mu, a.restitution)
            rows = []
            for t in range(a.throws):
                rows += sim.run(random_throw(a.throw_seed + t, kinds, tray.prof, a.aim, drop), a.t_max)
            results[(label, mu)] = rows
        sim.close()

    print(f"tray d={a.diameter:g} mm, slope {a.slope:g} deg, {a.throws} throws x {len(kinds)} dice, "
          f"restitution {a.restitution:g}, aim {a.aim:g} x inner radius, drop {drop[0]:g}-{drop[1]:g} mm")
    print(f"{'tray':16s} {'mu':>5s} {'n':>5s} {'flat%':>6s} {'blend%':>7s} {'slope%':>7s} {'out%':>6s} "
          f"{'tumble%':>8s} {'slid%':>6s} {'unsettled%':>11s}")
    for (label, mu), rows in results.items():
        n, c, uns = summarize(rows)
        tum = sum(r["tumbles"] > 0 for r in rows)
        slid = sum(r["slid"] for r in rows)
        print(f"{label:16s} {mu:5.2f} {n:5d} {100 * c['flat'] / n:6.1f} {100 * c['blend'] / n:7.1f} {100 * c['slope'] / n:7.1f} "
              f"{100 * c['out'] / n:6.1f} {100 * tum / n:8.1f} {100 * slid / n:6.1f} {100 * uns / n:11.1f}")
    if a.by_type:
        for title, f in (("flat%", lambda r: r["status"] == "flat"), ("tumble% (>= 1 up-face change on the cone)", lambda r: r["tumbles"] > 0),
                         ("slid% (flat after touching the cone without tumbling)", lambda r: r["slid"])):
            print(f"\n{title} by die type")
            print(f"{'tray':16s} {'mu':>5s} " + " ".join(f"{k:>5s}" for k in kinds))
            for (label, mu), rows in results.items():
                cells = []
                for k in kinds:
                    sub = [r for r in rows if r["kind"] == k]
                    cells.append(f"{100 * sum(f(r) for r in sub) / len(sub):5.0f}")
                print(f"{label:16s} {mu:5.2f} " + " ".join(cells))
    if a.plot:
        write_plot(a.plot, results, trays, mus)
        print(f"wrote {a.plot}")
    return 0


def write_plot(path, results, trays, mus):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    variants = list(trays)
    fig, axs = plt.subplots(len(mus), len(variants), figsize=(4.2 * len(variants), 4.2 * len(mus)), squeeze=False)
    col = {"flat": "tab:green", "blend": "gold", "slope": "tab:orange", "out": "tab:red"}
    for i, mu in enumerate(mus):
        for j, v in enumerate(variants):
            ax, pr = axs[i][j], trays[v].prof
            for r, ls in ((pr["r_f"], "-"), (pr["r_s"], ":"), (pr["R_in"], "-")):
                ax.add_patch(plt.Circle((0, 0), r, fill=False, ls=ls, color="gray"))
            for s, c in col.items():
                pts = [(r["x"], r["y"]) for r in results[(v, mu)] if r["status"] == s and r["r"] < 1.4 * pr["R_out"]]
                if pts:
                    ax.scatter(*zip(*pts), s=10, c=c, label=s)
            ax.set_aspect("equal")
            ax.set_xlim(-1.4 * pr["R_out"], 1.4 * pr["R_out"])
            ax.set_ylim(-1.4 * pr["R_out"], 1.4 * pr["R_out"])
            ax.set_title(f"{v}, mu={mu:g}")
    axs[0][0].legend(loc="upper left", fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=100)


if __name__ == "__main__":
    sys.exit(main())
