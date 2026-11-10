# Dark-sky measurements: a fixed camera and a Unihedron SQM-L

Night-sky brightness at dark sites in Crete, south-east Poland and central France, measured with one Nikon Z5, one lens, fixed settings and a Unihedron SQM-L.
The aim is to test how well the light-pollution atlas matches the real sky at the dark end.

Write-up with charts: [How dark is a dark sky? Measurements at seven sites](https://borama.co/dark-sky/)

This is a hobby project. If you measure with a Unihedron SQM anywhere in Europe, I'd be glad to compare readings, even without a camera. Please use GitHub Discussions.

**How this was made:** the photos and the meter readings are real and were taken by hand in the field. The scripts, much of the analysis and interpretation, and this text were written with Claude Code, an AI coding assistant. Errors were found and corrected along the way, and some may remain. The raw data is here so you can check it.

## What is here

| Path | Contents |
|---|---|
| `scripts/solve_sky.py` | Plate-solves a NEF against bright stars, then measures sky brightness |
| `scripts/sky_background.py` | Sky background, colour ratios and star sharpness per frame |
| `scripts/horizon_glow.py` | Maps a solved frame to altitude and azimuth, for reading horizon glow |
| `data/readings.csv` | Meter readings with time, site, position and atlas zone |
| `data/z5_frames.csv` | One row per measured Z5 frame, with the paired meter reading |
| `data/stars.csv` | Star positions and magnitudes used for plate solving |

**Raw files:** the NEFs are attached to [release v1](https://github.com/standinga/dark-sky-measurements/releases/tag/v1) as one zip per session, each with an `.xmp` sidecar per frame holding its GPS position. The NEFs themselves are unmodified.

| Zip | Frames |
|---|---|
| `2026-08-08-roztocze.zip` | 6 |
| `2026-09-11-agios-ioannis.zip` | 15 |
| `2026-09-19-rodakino.zip` | 5 |
| `2026-10-03-roztocze.zip` | 10 |
| `2026-10-06-roztocze.zip` | 12 |
| `2026-11-limousin.zip` | 27 |

Requirements: Python 3.11, `numpy`, `scipy`, `rawpy`, `astropy`, `matplotlib`, and `exiftool` on the path.

    python scripts/solve_sky.py --site 45.67281,2.09480 DSC_0523.NEF

## Field method

- **Camera:** Nikon Z5 with a Viltrox AF 16 mm f/1.8 Z, at f/1.8 (a few early frames at f/2.2), 30 s, ISO 800. The Z5 has no GPS, so position comes from a phone app, logged when each sequence starts.
- **Meter:** Unihedron SQM-L, whose lens narrows its field of view to about 20° (FWHM), read by hand and typed into a phone app with GPS and time.
- **Older meter:** some readings used an astronocnik.pl SQM-L, with or without its tube. They are labelled in the `meter` column and should never be pooled with Unihedron readings.

### How the readings were logged

Anyone reusing `readings.csv` needs to know this:

- Several readings were taken at every stop. When they nearly agreed, one value was logged, the repeated one (21.33, 21.33, 21.35 is logged as 21.33). So a stop with one row is not a single measurement, and row counts do not indicate quality.
- The zenith was always measured, but some stops also have extra readings deliberately aimed away from it. These are not marked. Use one median per stop, not pooled readings or the best value.
- Rows with a value in `excluded` were dropped from the results: twilight, or an outlier more than 0.3 mag from its stop's median.

## Sky brightness from the camera

Each frame is plate-solved, which gives its pointing and a per-frame zero point from the stars. Zero points scatter too much to compare sites, so the site comparison uses a fixed scale for this camera and these settings:

    mu = 21.62 + 2.5 log10(296.7 / ADU)

ADU is the raw green sky background in a 15° cone (half-angle), corrected for cos⁴ vignetting and taken at the 25th percentile of blocks to reject stars. Frames at f/2.2 are multiplied by 1.494 to match f/1.8. The constants come from one reference frame, so this is a relative scale for this camera, not a true SQM value. The per-frame zero point is used only as a check on transparency.

**Using another camera:** the 21.62 and 296.7 will not transfer. Derive your own from one well-solved, clear, dark frame at your fixed settings.

### Columns in `z5_frames.csv`

- `like_for_like`: altitude at least 60°, at least 14° from the galactic plane, sun below −18° and moon below the horizon. Only these frames should be compared between sites.
- `mu_z5`: sky brightness on the fixed camera scale above.
- `meter_stop_median`: median Unihedron reading at the same stop, within 150 m and 60 minutes.
- `meter_values_logged`: how many logged values that median uses.
- `meter_gap_min`: minutes between the frame and the nearest of those readings.

## Results so far

| Site | Z5 fixed scale, like-for-like | Unihedron SQM-L |
|---|---|---|
| Saint-Merd-de-Lapleau (Corrèze, France) | 21.99 to 22.06 | 21.33 to 21.40 |
| Agios Ioannis, Sfakia (Crete, Greece) | 22.00 | astronocnik only |
| Coiral (Corrèze, France) | 21.91 | 21.25 to 21.36 |
| Millevaches plateau (Limousin, France) | 21.87 | 21.16 to 21.27 per stop |
| Rodakino, Rethymno (Crete, Greece) | 21.79 | astronocnik only |
| Near Horyniec-Zdrój (Roztocze, Poland) | 21.78 | 21.26 to 21.29 around midnight |

The Millevaches plateau is the Millevaches en Limousin International Dark Sky Reserve. Its readings come from nine stops across the plateau on one night.

These results come from few frames and few nights. Sites were measured at different times of night, in different seasons and at different levels of solar activity, and Z5 and Unihedron numbers are on different scales.

## Open questions

Help is very welcome on any of these:

1. **How well does the atlas match at the dark end?** A difference of 2.9 vs 3.9 on [lightpollutionmap.info](https://www.lightpollutionmap.info), from the Lorenz atlas, predicts about 0.15 mag, close to what was measured. Is that general?
2. **Can artificial light be separated from airglow at the zenith?** Horizon glow is clearly artificial because it stays at a fixed bearing. A single zenith number cannot be split. Timing lamp switch-offs, many nights at one site, same-night comparisons and spectra are all candidates.
3. **Does sky colour follow solar activity?** The red/green ratio rose with F10.7 across these sessions, but site, season and solar flux are fully confounded here.
4. **Same protocol elsewhere.** If you have a Unihedron and a camera you can lock to fixed settings, comparable data from other sites would be the most useful thing of all.

## Location data

The Roztocze site is given as near Horyniec-Zdrój, with its atlas zone. To plate-solve its frames, use `--site 50.191,23.363`: the difference is too small to affect the solve.

## Licence

TODO: licence for the scripts and data. `data/stars.csv` is CC BY-SA 4.0 (see Credits).

## Credits

- Star data derives from the [HYG Database](https://github.com/astronexus/HYG-Database) by David Nash, CC BY-SA 4.0. `data/stars.csv` is shared under the same licence.
- Light-pollution zones: David Lorenz's light-pollution atlas (2025).

Written with [Claude Code](https://claude.com/claude-code) by Anthropic, from my own measurements and photos.
