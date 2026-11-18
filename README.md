# dicetray

Generator for a round, 3D-printable dice tray. The sloped ring funnels thrown dice into a flat centre. Output is STL, or a 3MF with three colour parts (base, slope, wall). A pybullet simulation throws polyhedral dice at the generated trays.

Two variants:

- `smooth`: flat centre, blend arc, constant 18° cone, fillet, wall.
- `random`: the same, plus seeded soft bumps and dimples on the slope. The generator measures the steepest-descent angle on the finished surface and refuses (or scales the features down) if a die could come to rest on the slope.

Design rationale, numbers and test plan: [DESIGN-v1.md](DESIGN-v1.md). The header of each Python file holds the decisions and invariants for that file.

## Install

```
pip install numpy scipy pybullet matplotlib pytest
```

`numpy` is enough for `dicetray.py`. `scipy` and `pybullet` are for `dicesim.py`, `matplotlib` for the PNG previews.

## Generate

```
python3 dicetray.py --variant smooth --diameter 150
python3 dicetray.py --variant random --diameter 200 --seed 7 --preview preview.png
python3 dicetray.py --variant random --diameter 150 --format 3mf --filaments 1,2,3
python3 dicetray.py --diameter 200 --format both --zone-stls
python3 dicetray.py --diameter 200 --slope 20 --floor 2 --felt 50
```

- The centre floor is 1.0 mm by default (5 layers at 0.2 mm). A felt recess needs a thicker `--floor`.
- `--format 3mf` writes one object with three parts and per-part filament numbers (`--filaments base,slope,wall`). `--zone-base-r` and `--zone-wall-z` move the colour borders.
- Use the 3MF layout as a starting point: it follows the part/extruder keys read by OrcaSlicer-family loaders and was not opened in ElegooSlicer. `--zone-stls` writes the same parts as separate STLs.

## Generator options

All lengths are in mm, angles in degrees. Defaults are in brackets.

**Shape**

- `--variant smooth|random` (smooth): plain cone, or cone plus bumps and dimples.
- `--diameter` (200): outer diameter. The bed of the Centauri Carbon 2 is 256 mm; the generator warns above 250.
- `--flat` (0.41 × diameter): diameter of the flat centre where dice should end up.
- `--slope` (18): cone angle. A die stays on an incline if the friction coefficient is above tan(angle), so a steeper cone collects better but throws dice around more. See DESIGN-v1.md, section 5, for simulation results.
- `--blend-radius` (30): radius of the concave arc between the flat centre and the cone. A larger radius makes a longer, gentler transition where dice can stop; a smaller one makes a sharper corner. The arc is only radius × sin(slope) long horizontally (9 mm for 30 mm at 18°).
- `--wall-fillet` (6): radius of the arc between the cone and the wall, so no sharp corner holds a die.
- `--wall-height` (12): straight wall above the fillet. Raise it if dice leave the tray.
- `--rim` (4): wall thickness.
- `--floor` (1.0): material under the flat centre. 1.0 mm is 5 layers at 0.2 mm.
- `--felt` (0 = off) and `--felt-depth` (1.0): diameter and depth of a recess for a felt disc in the centre. At least 0.6 mm of floor must remain, so use e.g. `--felt 50 --floor 2`.

**Bumps and dimples (random variant only)**

- `--seed` (1): the same seed gives the same tray.
- `--features` (about one per 600 mm² of the slope): number of bumps and dimples.
- `--amp` (1.2): maximum height of a bump or depth of a dimple. Each feature gets 50–100 % of it.
- `--sigma` (6.5): width of a feature. Each feature varies by ±20 %.
- `--dimple-ratio` (0.5): fraction of features that are dimples.
- `--min-slope` (8): the steepest-descent angle that must hold everywhere on the cone, so no feature can hold a die.
- `--no-autofit`: by default, a seed that breaks `--min-slope` is scaled down in 10 % steps until it passes. With this flag the run stops with an error instead.
- `--stubs` (0): number of steeper rounded pegs in the middle of the slope. They can hold a die back, so the check only warns for them.

**Deflectors (both variants)**

Off by default. They break up a glide so that fast d4/d6 tumble. Results from the simulation are in DESIGN-v1.md, section 6; read them before printing, the gain is in collecting dice, not in a proven tumble.

- `--deflectors none|diamonds|steps|both` (none): `diamonds` = short rounded bars (louvers) in staggered rings, tilted against the downhill direction, with the tilt alternating from ring to ring. `steps` = concentric drops on the cone. `both` = bars on top of steps, with the risers midway between the bar rings.
- `--deflector-pitch` (10): radial distance between rings and risers. The number of rings follows from the width of the slope, so a 150 mm tray at 18° fits one or two rings and a 200 mm tray three.
- `--bar-length` (14), `--bar-height` (3.0), `--bar-angle` (25): size of a bar (length includes the rounded ends; the base is 4 mm wide) and its tilt against the downhill direction.
- `--bar-coverage` (0.5): largest share of a ring's circumference that its bars may block, so lanes stay open.
- `--step-height` (2.0): drop of one riser. The treads between risers keep the `--min-slope` descent, and a step height that breaks it is scaled down like the bumps.
- Bars are dams by design, so the no-trap rule of the bumps cannot hold for them. Instead the generator floods the finished height field from the inner edge and removes bars that close a basin deeper than 0.5 mm. This covers the frictionless case only; a die can still stick against a bar.
- The bumps of a seed do not change when deflectors are switched on.

**Output**

- `--format stl|3mf|both` (stl): `3mf` writes one object with three colour parts (base, slope, wall).
- `--filaments` (1,2,3): filament numbers for base, slope and wall.
- `--zone-base-r` (end of the blend arc): outer radius of the base part. It can range from the flat edge to the end of the blend.
- `--zone-wall-z` (top of the fillet): height where the wall part starts. It can range from the top of the fillet to 1 mm below the rim.
- `--zone-stls`: also write the three parts as separate STLs, for a manual multi-part import.
- `--bed` (256): bed size, used to centre the 3MF.
- `-o, --output`: output name; the suffix follows `--format`. Default: `dicetray_<variant>_d<diameter>[_s<seed>]`.
- `--preview file.png`: top-view image of the surface and the height offset of the features.
- `--res` (0.75) and `--sectors` (720): sampling along the profile and around the circle. Raise them only if facets show.

## Simulate

Run it without `-I`, so the script directory is on the import path.

```
python3 dicesim.py --diameter 150 --throws 12 --mu 0.2,0.4,0.6 --by-type --plot sim.png
```

The friction coefficient on printed plastic is not measured, so `--mu` is a sweep. Dice sizes, mass and restitution are illustrative. Use the numbers to compare trays and see trends, not as a prediction for a specific dice set.

## Simulation options

The tray options `--diameter` (150 here), `--flat`, `--slope`, `--blend-radius`, `--floor`, `--wall-height`, `--seed`, `--amp`, `--sigma` and `--features` work as above.

- `--variants` (smooth,random): trays to compare.
- `--deflectors` (none): comma-separated list from `none,diamonds,steps,both`; every variant runs with each entry, e.g. `--deflectors none,diamonds,both`. The bar and step sizes use the generator options `--deflector-pitch`, `--step-height`, `--bar-length`, `--bar-height`, `--bar-angle` and `--bar-coverage`.
- `--throws` (20): number of throws per tray and friction value. Each throw drops all selected dice. Ten throws of six dice give 60 dice per cell and a standard error of about 6 points.
- `--throw-seed` (100): the same seed replays the same throws for every tray, so differences come from the tray.
- `--dice` (d4,d6,d8,d10,d12,d20) and `--dice-scale` (1.0): which dice, and a size factor.
- `--mu` (0.2,0.4,0.6): effective friction coefficients, one run each. Dice on printed plastic are not measured, so keep this a sweep.
- `--restitution` (0.3): effective bounciness.
- `--aim` (0.8): each die is aimed at a random point within this fraction of the inner radius, so throws land inside the tray.
- `--drop` (20,60): release height range of the lowest die above the rim, in mm. Each further die starts 24 mm higher.
- `--t-max` (3.0): simulated seconds per throw. The result is the position at that time; a longer time lets creeping dice settle.
- `--sim-res` (2.5) and `--sim-sectors` (240): mesh sampling of the tray in the simulation. A coarser mesh runs faster. With deflectors use about 1.0 and 480, otherwise the 4 mm bars are only a few samples wide.
- `--by-type`: also print flat, tumble and slid shares per die type.
- `--plot file.png`: final resting positions per variant and friction value.

The table shows the share of dice per zone: flat centre, blend arc, cone (slope), and out (beyond the inner wall). Two more columns look at the way down: `tumble%` is the share of dice whose up-face changed at least once while on the cone (sampled every 12.5 ms), `slid%` the share that ended on the flat after touching the cone without any change. A throw also changes the up-face on impact, so `tumble%` is high even on a plain cone and does not separate trays well; compare `flat%` and `slid%`.

## Test

```
python3 -m pytest -q
```

## License

MIT, see [LICENSE](LICENSE).
