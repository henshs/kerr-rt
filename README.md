# Relativistic Ray Tracing of Thin Accretion Disks Around Kerr Black Holes

A backward ray-tracing code for computing images of a Novikov–Thorne thin Keplerian accretion disk around a Kerr black hole.

Photons are traced backward from an observer's image plane through the Kerr spacetime until they intersect the disk. The code accounts for relativistic effects such as gravitational lensing, Doppler boosting, gravitational redshift, and frame dragging.

This code is developed by Jan Schee (jan.schee@physics.slu.cz) and Sudipta Hensh(sudipta.hensh@aei.mpg.de).

## Components

- **kdfirewall** – backward Kerr ray tracer
- **kdimage** – image construction utility
- **kerrrms** – computes the Kerr ISCO radius

## Requirements

- C++ compiler
- GNU Make
- GNU Scientific Library (GSL)
- GNUPlot
  
## Build

```bash
bash install.sh
```

This creates the `simulations/` directory structure and copies all executables into `simulations/bin/`

## Run

A sample run script `run.sh` is automatically generated:

The main simulation parameters can be modified directly in `run.sh` according to user's needs:

| Parameter | Description |
|------------|------------|
| `a` | Black-hole spin parameter |
| `tho` | Observer inclination angle in degrees (`0°` = face-on, `90°` = edge-on) |
| `rmax` | Outer radius of the emitting disk |
| `alphamin`, `alphamax` | Horizontal field of view on the observer sky |
| `betamin`, `betamax` | Vertical field of view on the observer sky |
| `w`, `h` | Output image resolution |
| `uo` | Observer location parameter (`u = 1/r`) |

After making changes, run:

```bash
bash run.sh
```

to generate a image.

The image is written to:

```text
png/disk_a_<spin>_tho_<inclination>.png
```

## Output Data

The generated image data file contains:

```text
alpha    beta    intensity(cgs)    frequency_shift
```

where

- **alpha** : horizontal observer sky coordinate
- **beta** : vertical observer sky coordinate
- **intensity** : observed specific intensity in CGS units
- **frequency_shift** : redshift/blueshift factor (g-factor)

By default, the supplied GNUPlot script visualizes the intensity column.

To visualize the frequency shift instead, modify the plotting command in `run.sh` from:

```gnuplot
plot "$fimagedata" u 1:2:3 w image
```

to

```gnuplot
plot "$fimagedata" u 1:2:4 w image
```

## References

- Kerr, R. P., Phys. Rev. Lett. 11, 237 (1963).
- Bardeen, J. M., Press, W. H. & Teukolsky, S. A. (1972), ApJ 178, 347–370.
- Novikov, I. D. & Thorne, K. S. (1973), *Astrophysics and black holes*, in *Les Houches Summer School: Black Holes*, 343–550.

## License

This project is open-source. You are free to use, modify, and distribute it with attribution.
