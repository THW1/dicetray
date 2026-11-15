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

## Simulate

Run it without `-I`, so the script directory is on the import path.

```
python3 dicesim.py --diameter 150 --throws 12 --mu 0.2,0.4,0.6 --by-type --plot sim.png
```

The friction coefficient on printed plastic is not measured, so `--mu` is a sweep. Dice sizes, mass and restitution are illustrative. Use the numbers to compare trays and see trends, not as a prediction for a specific dice set.

## Test

```
python3 -m pytest -q
```

## License

MIT, see [LICENSE](LICENSE).
