# Round dice tray with centre funnel — design document v1.3

Hobby project: a round, 3D-printed dice tray whose sloped ring drives thrown dice into a flat centre. Two variants: **smooth** (plain cone) and **random** (cone plus soft bumps and dimples that add twists but never trap a die). Output as STL, or as a 3MF with three colour parts.

Date: 2026-11-15 · Printer: Elegoo Centauri Carbon 2 (build volume 256 × 256 × 256 mm, per project FAQ) · Slicer: ElegooSlicer · Code: `dicetray.py` (generator), `dicesim.py` (dice simulation), repo https://github.com/THW1/dicetray

Changes since v1.2: optional deflectors (tilted bars, concentric steps, or both) and a tumble metric in the simulation (section 6).

Changes since v1.1: centre floor reduced to 1.0 mm (5 layers at 0.2 mm); 3MF export with three colour zones; rigid-body dice simulation, which shows that the 18° default is too gentle for friction values of 0.4 and above (section 5).

## 1. Geometry

Cross-section (revolved), from the centre outwards:

1. Flat centre, **1.0 mm thick** (5 layers at 0.2 mm). Optional felt recess; it needs a thicker `--floor`.
2. Blend arc (radius 30 mm), tangent to the flat, so there is no edge where the slope starts.
3. Cone at constant 18° (default).
4. Fillet (radius 6 mm) into the wall, so no sharp corner can hold a die.
5. Vertical wall, 12 mm above the fillet, rim 4 mm wide, 0.6 mm bottom chamfer.

Computed by the generator (defaults, flat diameter = 0.41 × outer diameter):

| Outer ⌀ | Flat ⌀ | Cone run | Cone rise | Height | Solid volume |
|---|---|---|---|---|---|
| 150 mm | 61 mm | 27 mm | 8.8 mm | 29.0 mm | 148 cm³ |
| 200 mm | 82 mm | 42 mm | 13.5 mm | 33.7 mm | 319 cm³ |
| 220 mm | 90 mm | 48 mm | 15.5 mm | 35.6 mm | 416 cm³ |

The solid volume is the full body; the real filament use is lower and depends on infill. Start with 150 mm and 200 mm; 220 mm uses the same defaults.

## 2. Slope and friction

A die stays on an incline if the friction coefficient μ is above tan(angle). At 18°, tan = 0.325. Plastic on plastic is assumed at μ ≈ 0.3–0.5 [UNVERIFIED, not measured], so 18° is at the edge of what works. The simulation (section 5) confirms this. A steeper cone (`--slope 25` to `30`) and a smaller blend radius (`--blend-radius 10`) gave clearly better collection. These are not yet the defaults; decide after a test print.

Felt belongs in the flat centre only. Felt on the slope adds friction and defeats the funnel.

## 3. Random variant: bumps and dimples without traps

Features are Gaussian bumps (+) and dimples (−), h = A·exp(−d²/2σ²), placed with seeded dart throwing (minimum spacing 2.4σ) in the slope band, roughly 50/50.

| Parameter | Default |
|---|---|
| Amplitude A | up to 1.2 mm (each feature 50–100 % of it) |
| Width σ | 6.5 mm (each feature ±20 %) |
| Density | about one feature per 600 mm² of the band |
| Envelope | smoothstep from 0 at the blend end to full at half the cone run, fading to 0 again 8 mm before the fillet |

**No-trap rule.** If the radial derivative dz/dr is greater than 0 everywhere on the cone, the surface has no local minimum. The steepest flank of a Gaussian is 0.607·A/σ. The envelope adds gradient and overlapping features add up, so the generator does not rely on the estimate:

- It measures the smallest descent angle on the finished vertices and requires at least 8° (`--min-slope`).
- If a seed fails, it scales all amplitudes down in 10 % steps (autofit) until the rule holds. Over seeds 1–20 at 150, 200 and 220 mm, the measured minimum descent was 8.0–9.9° and autofit never went below 0.73 of the amplitude.
- If autofit is off or cannot fix it, the run stops with an error and no file is written.

Dimples are wide and shallow so a d4 tip or d20 corner does not lodge in them. `--stubs N` adds steeper pegs; they can hold a die back, so the descent check only warns for them.

The simulation shows no clear difference in collection between smooth and random (differences of up to 8 points, consistent with sampling noise). The bumps change the paths, not the end result.

## 4. Multi-colour 3MF

`--format 3mf` (or `both`) writes one object made of three touching parts. Each part gets its own filament number (`--filaments base,slope,wall`, default 1,2,3).

| Part | Region | Volume at 150 / 200 mm |
|---|---|---|
| base | flat centre and blend arc, r ≤ end of blend, full height down to the bed | 6 / 9 cm³ |
| slope | ring outside the base, up to the top of the fillet; includes the lower outer wall | 120 / 280 cm³ |
| wall | everything above the top of the fillet: upper wall and rim | 22 / 30 cm³ |

The borders are adjustable: `--zone-base-r` (between flat edge and blend end) and `--zone-wall-z` (on the vertical wall, from the top of the fillet to 1 mm below the rim). The cuts are planes and cylinders of the revolve, so each part is a watertight solid and the three volumes add up to the whole body exactly (tested).

With the 1 mm floor, the base part is a 5-layer disc: one colour for the first five layers and the centre.

**Format and status.** The file is a 3MF core model (three mesh objects plus one object with components, centred on the bed) with `Metadata/model_settings.config` holding part name, matrix and `extruder` per part. This layout was checked against the loader source of OrcaSlicer and ElegooSlicer (`bbs_3mf.cpp`), and the file loads in an independent reader (trimesh: three watertight meshes). It has **not** been opened in ElegooSlicer. If the import shows one grey object, use `--zone-stls` and import the three STLs as one multi-part object. The slicer keeps a part's extruder only if the project has that many filaments; with the Canvas set to four colours, 1–3 are available.

## 5. Dice simulation

`dicesim.py` throws d4, d6, d8, d10, d12 and d20 (typical sizes, ideal convex polyhedra, mass from an acrylic density) into the tray with pybullet. Throws are aimed at random points inside the tray, released 20–60 mm above the rim (first die; each further die 24 mm higher), with random orientation and spin. The same seeded throws are replayed for every variant and every μ. The classification uses the centre of mass at the end (3 s): **flat** = on the flat centre, **blend** = on the blend arc, **cone** = on the slope, **out** = beyond the inner wall.

Checks of the engine setup (tests): a die rests within 0.4 mm of half its edge above the floor; on an 18° incline it stays at μ = 0.6 and slides at μ = 0.15; the same throw gives identical results twice. In a probe sweep with μ from 0.3 to 1.5, a die stayed put from μ = 0.4 upward and crept down at 0.3, consistent with tan 18° = 0.325 within the engine's accuracy.

Result: share of dice in % (flat / on the cone), 12 throws × 6 dice = 72 dice per cell, 150 mm unless noted, blend radius 30 mm unless noted. The standard error of one cell is about ±6 points. No die left the tray in any run.

| Tray | μ = 0.2 | μ = 0.4 | μ = 0.6 |
|---|---|---|---|
| 18°, smooth | 72 / 3 | 49 / 24 | 33 / 36 |
| 18°, random | 68 / 1 | 46 / 24 | 38 / 33 |
| 25°, smooth | 82 / 0 | 51 / 1 | 54 / 15 |
| 25°, random | 79 / 0 | 51 / 4 | 49 / 17 |
| 30°, smooth | 81 / 3 | 60 / 1 | 46 / 4 |
| 30°, random | 82 / 0 | 51 / 1 | 42 / 6 |
| 18°, smooth, 200 mm | 65 / 8 | 40 / 35 | 39 / 38 |
| 18°, random, 200 mm | 65 / 4 | 42 / 29 | 39 / 31 |
| 25°, smooth, blend 10 mm | – | 65 / 18 | 57 / 32 |
| 30°, smooth, blend 10 mm | – | 75 / 10 | 64 / 22 |

Reading the table:

- Friction dominates. At μ = 0.2 about two thirds to four fifths of the dice reach the flat; at 0.4–0.6 an 18° cone leaves 24–38 % of them on the slope.
- 25–30° removes most dice from the cone, but with the 30 mm blend arc many stop on the blend, where the slope falls to zero at the flat edge (blend share 17–53 % in those rows).
- A small blend radius (10 mm) with 30° gave the best result in these runs: 75 % flat at μ = 0.4, 64 % at 0.6, against 60 % and 46 % with the 30 mm blend. Why some dice still stay on the cone at 30° was not examined.
- The 200 mm tray was not better than 150 mm at 18°.
- Check of the 3 s snapshot: 150 mm, 18°, smooth, μ = 0.4 with 10 s simulated gave 51 / 19 against 49 / 24, so the numbers are stable. Per cell, 25–83 % of the dice belonged to throws that had not fully come to rest at 3 s (dice creeping on marginal slopes); in the 10 s check it was 17 %.

Limits: sizes, density, restitution (0.3) and the throw model are illustrative; μ on printed plastic is unknown, which is why μ is a sweep. Use the numbers to compare trays and see trends, not as a prediction for a specific dice set. A confirmed test print with your dice beats the simulation.

## 6. Deflectors (optional)

Question: can obstacles like those on a roulette wheel make d4/d6 tumble instead of gliding down the cone? A plain bar across the flow is a dam and holds dice, so two shapes were built (`--deflectors`):

- **steps**: concentric risers of 2 mm. Moving inwards, a die drops off each riser. The treads between them rise outwards at the cone angle minus the share used by the risers, so the radial descent stays ≥ `--min-slope` and the no-trap rule of section 3 holds exactly. Steps and bumps share the autofit scale.
- **diamonds**: short rounded bars (14 mm long, 4 mm wide at the base, 3 mm high) in staggered rings, 25° off the downhill direction, one tilt per ring, alternating between rings. At most half of a ring is covered, so there is always a lane. Bars are dams on purpose, so the radial rule cannot hold for them. The generator floods the finished height field from the inner edge (priority flood) and removes bars that close a basin deeper than 0.5 mm. This guards the frictionless limit only. Whether a die sticks against a bar depends on μ and the angle; the geometric condition I derived for sliding along a bar, tan(slope)·(cos φ − μ_bar·sin φ) > μ, is my own and was not checked against a source.
- **both**: bars on top of steps, risers midway between the bar rings.

Chevron (V) bars were rejected because the funnel between two bars is narrower than a d6, and a free random bar layout because neighbouring bars close pockets. The bumps of a seed are the same with and without deflectors.

The simulation got two metrics: **tumble** = the up-face changed at least once between samples (every 12.5 ms) while the die touched the cone, and **slid** = the die ended on the flat after touching the cone without any change. Setup: smooth tray, d4 and d6 only, 40 throws each = 80 dice per cell (standard error about ±5 points per cell, ±8 per die type), same throws for every tray, 6 s simulated, finer tray mesh (1 mm, 480 sectors). Cells show flat % (d4 / d6) and, after the slash, the share left on the cone.

| Tray | μ | none | steps | diamonds | both |
|---|---|---|---|---|---|
| 200 mm, 30°, blend 10 | 0.3 | 61 (48/75) / 8 | 69 (50/88) / 8 | 83 (65/100) / 14 | 76 (52/100) / 18 |
| 200 mm, 30°, blend 10 | 0.5 | 46 (22/70) / 11 | 51 (22/80) / 24 | 64 (38/90) / 26 | 64 (32/95) / 20 |
| 200 mm, 18°, blend 30 | 0.2 | 43 (42/42) / 1 | 46 (45/48) / 3 | 59 (45/72) / 30 | 66 (52/80) / 18 |
| 200 mm, 18°, blend 30 | 0.3 | 35 (32/38) / 9 | 24 (18/30) / 15 | 35 (22/48) / 36 | 40 (20/60) / 28 |
| 150 mm, 30°, blend 10 | 0.3 | 78 (68/88) / 4 | – | 85 (70/100) / 11 | 81 (62/100) / 10 |
| 150 mm, 30°, blend 10 | 0.5 | 56 (40/72) / 5 | – | 64 (40/88) / 14 | 71 (55/88) / 14 |

Reading the table:

- Bars raised the share of d6 on the flat in all six comparisons (+10 to +30 points, up to 100 % at μ = 0.3 and 30°). For d4 the gain appeared only at 200 mm / 30° (+16 to +17 points); elsewhere it was +3, −10, +2 and 0, which is within noise.
- The price is more dice staying on the cone: 1–11 % on the plain tray against 10–36 % with bars. At 18° the bars take dice mostly from the blend arc (56 % on the blend with the plain tray, 11–33 % with bars), so the net gain on the flat is smaller than the d6 numbers suggest.
- **The metrics do not show forced tumbling.** The up-face changes for 56–79 % of the dice on every tray, including the plain cone, because the throw itself and the first bounce already turn the die. The share of dice that glided down without any change (`slid`) was 2–18 % everywhere, and the differences between trays are within noise. So the data support better collecting of d6, not the claim that the bars make the dice tumble.
- Steps alone did nothing measurable (−11 to +8 points, within noise); at 18° and μ = 0.3 the result was worse. Both combined was not better than diamonds alone (differences of −6 to +7 points). I see no reason to prefer `both` on this evidence; it only adds the risers.
- In 17–57 % of the cells' throws some die was still moving at 6 s, so the classification is a snapshot, as in section 5.
- Not examined: other bar angles, heights and lengths, or whether a die can rest across two bars. Defaults are first guesses.

Suggestion: first test print without deflectors, then a 200 mm tray with `--deflectors diamonds` and see whether d6 and d4 collect better with your dice and your surface.

## 7. Printing

- Print upright, flat bottom down, no supports.
- PETG or PLA, 0.4 mm nozzle. Walls 3–4, gyroid infill 8–10 %. The 1 mm centre is five solid layers: set at least 5 bottom and 5 top layers so it does not turn into infill.
- Slope surface: variable layer height 0.08–0.12 mm reduces visible steps on the shallow blend; the flat centre can use 0.2 mm.
- Optional felt disc in the centre: `--felt 50 --floor 2` cuts a 1 mm recess, leaving a 1 mm floor.

## 8. Usage

```
python3 dicetray.py --variant smooth --diameter 150
python3 dicetray.py --variant random --diameter 200 --seed 7 --preview preview.png
python3 dicetray.py --variant random --diameter 150 --format 3mf --filaments 1,2,3
python3 dicetray.py --diameter 150 --slope 30 --blend-radius 10 --format both --zone-stls
python3 dicesim.py --diameter 150 --throws 12 --mu 0.2,0.4,0.6 --by-type --plot sim.png
python3 dicetray.py --diameter 200 --deflectors diamonds --preview defl.png
python3 dicesim.py --variants smooth --deflectors none,diamonds,both --diameter 200 --slope 30 --blend-radius 10 \
  --dice d4,d6 --throws 40 --mu 0.3,0.5 --t-max 6 --sim-res 1.0 --sim-sectors 480 --by-type
```

Every generator run prints diameter, height, volume, watertight status, measured minimum descent angle and the autofit scale. Run `dicesim.py` without `-I` so the script directory is importable. A 12-throw sweep over two variants and three μ values takes about 5–8 minutes on two cores; the deflector sweep above (3–4 trays, two μ, 80 dice per cell) about 15–25 minutes.

## 9. Test plan

1. Print the 150 mm smooth tray with `--slope 30 --blend-radius 10`, and one with the 18° default if you want the comparison. Throw your usual dice set 50 times from about 30 cm height and count the dice that end outside the flat centre.
2. If many dice stop on the slope, raise the slope or polish the slope surface (it lowers friction). If they leave the tray, raise `--wall-height`.
3. Try the random variant with the settings that worked.
4. Open the 3MF in ElegooSlicer, check that three parts with filaments 1–3 appear; if not, use `--zone-stls`.
5. Then 200 mm, and 220 mm if needed.

## 10. Assumptions and limits

- μ range, dice dimensions and the throw model are assumptions until a test print says otherwise.
- The no-trap rule covers the radial direction on the cone and a die treated as a point.
- The 3MF is untested in ElegooSlicer (section 4).
- A sector or wedge export for quick tests is not implemented; use the small diameter.

## 11. Files

- `dicetray.py` — generator (STL, 3MF, preview); the file header holds the design decisions and invariants.
- `dicesim.py` — pybullet dice simulation.
- `test_dicetray.py` — 54 tests: watertight mesh, volume against the analytic value, seed determinism, no-trap check, autofit repair, felt recess, STL size, 5-layer floor, zones watertight and adding up to the body, 3MF structure and independent reload, deflectors (watertight, unchanged outside the cone, step descent and autofit, bars seeded and coverage-limited, flood check, repair of a closed ring).
- `test_dicesim.py` — 18 tests: rest height, friction threshold on the cone, determinism, dice geometry, face counts and up-face detection, tumble/slid flags.
