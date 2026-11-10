#!/usr/bin/env python3
"""
Measure night-sky surface brightness (mag/arcsec^2) from a Nikon NEF.

Self-calibrating: plate-solves the frame against a bright-star catalog, derives the
photometric zero point from the identified stars, then converts the sky background
in the same ADU units into mag/arcsec^2.

Usage: solve_sky.py --site LAT,LON [--cone DEG] FRAME.NEF [FRAME2.NEF ...]
The Z5 has no GPS, so --site is required.
"""

import csv
import math
import os
import subprocess
import sys

import numpy as np
import rawpy
from scipy.spatial import cKDTree

LAT, LON = None, None              # always set from --site
SENSOR_WIDTH_MM = 35.9
SENSOR_WIDTH_PX = 6016
EXTINCTION_K = 0.20                # mag per airmass, green, clear dark site
SKY_CONE_DEG = 15.0                # half-angle of the measured sky patch; match the
                                   # meter you compare against (SQM-L ~10, SQM ~42)
CATALOG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "stars.csv")

D2R = math.pi / 180.0


def load_catalog():
    with open(CATALOG_PATH, newline="") as handle:
        rows = list(csv.DictReader(handle))
    ra_hours = np.array([float(row["ra_hours"]) for row in rows])
    dec_deg = np.array([float(row["dec_deg"]) for row in rows])
    mag = np.array([float(row["magnitude"]) for row in rows])
    return ra_hours * 15.0, dec_deg, mag


def exif(path):
    fields = ["DateTimeOriginal", "OffsetTimeOriginal", "PitchAngle", "RollAngle", "YawAngle",
              "ExposureTime", "ISO", "FNumber", "FocalLength"]
    out = subprocess.run(
        ["exiftool", "-T", "-n"] + [f"-{f}" for f in fields] + [path],
        capture_output=True, text=True, check=True).stdout.strip().split("\t")
    stamp = out[0]
    year, month, day = (int(v) for v in stamp[:10].split(":"))
    hour, minute, second = (int(v) for v in stamp[11:].split(":"))
    # The camera writes its own zone (+02:00 in Poland, +03:00 on Crete); a missing tag must stop the run, a wrong hour silently rotates the sky.
    offset_sign = -1.0 if out[1].startswith("-") else 1.0
    offset_hours, offset_minutes = (int(v) for v in out[1].lstrip("+-").split(":"))
    offset = offset_sign * (offset_hours + offset_minutes / 60)
    return {
        "utc": (year, month, day, hour - offset + minute / 60 + second / 3600),
        "pitch": float(out[2]), "roll": float(out[3]), "yaw": float(out[4]),
        "exposure": float(out[5]), "iso": float(out[6]),
        "fnumber": float(out[7]), "focal_mm": float(out[8]),
    }


def julian_day(year, month, day, hour):
    if month <= 2:
        year -= 1
        month += 12
    a = year // 100
    b = 2 - a + a // 4
    return (math.floor(365.25 * (year + 4716)) + math.floor(30.6001 * (month + 1))
            + day + hour / 24.0 + b - 1524.5)


def star_altaz(ra_deg, dec_deg, jd):
    # J2000 -> epoch of date. 24 years of precession is ~0.33 deg, well above the
    # matching tolerance, so it cannot be skipped.
    t = (jd - 2451545.0) / 36525.0
    zeta = (2306.2181 * t + 0.30188 * t * t) / 3600.0 * D2R
    z = (2306.2181 * t + 1.09468 * t * t) / 3600.0 * D2R
    theta = (2004.3109 * t - 0.42665 * t * t) / 3600.0 * D2R
    ra, dec = ra_deg * D2R, dec_deg * D2R
    a = np.cos(dec) * np.sin(ra + zeta)
    b = np.cos(theta) * np.cos(dec) * np.cos(ra + zeta) - np.sin(theta) * np.sin(dec)
    c = np.sin(theta) * np.cos(dec) * np.cos(ra + zeta) + np.cos(theta) * np.sin(dec)
    ra = np.arctan2(a, b) + z
    dec = np.arcsin(c)

    gmst = (280.46061837 + 360.98564736629 * (jd - 2451545.0)) % 360.0
    hour_angle = ((gmst + LON) * D2R - ra)
    lat = LAT * D2R
    alt = np.arcsin(np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(hour_angle))
    az = np.arctan2(-np.sin(hour_angle) * np.cos(dec),
                    np.cos(lat) * np.sin(dec) - np.sin(lat) * np.cos(dec) * np.cos(hour_angle))
    alt_deg = alt / D2R
    visible = alt_deg > 1.0
    refraction = np.zeros_like(alt_deg)
    refraction[visible] = 1.02 / np.tan(
        (alt_deg[visible] + 10.3 / (alt_deg[visible] + 5.11)) * D2R) / 60.0
    return alt_deg + refraction, (az / D2R) % 360.0


def unit_vectors(alt_deg, az_deg):
    alt, az = alt_deg * D2R, az_deg * D2R
    return np.stack([np.cos(alt) * np.sin(az),   # east
                     np.cos(alt) * np.cos(az),   # north
                     np.sin(alt)], axis=-1)      # up


def camera_basis(alt_deg, az_deg, roll_deg, handedness):
    view = unit_vectors(np.array(alt_deg), np.array(az_deg))
    zenith = np.array([0.0, 0.0, 1.0])
    up = zenith - np.dot(zenith, view) * view
    up /= np.linalg.norm(up)
    right = np.cross(view, up) * handedness
    roll = roll_deg * D2R
    return (right * math.cos(roll) + up * math.sin(roll),
            -right * math.sin(roll) + up * math.cos(roll),
            view)


def project(star_vectors, alt_deg, az_deg, roll_deg, focal_cells, handedness):
    right, up, view = camera_basis(alt_deg, az_deg, roll_deg, handedness)
    depth = star_vectors @ view
    with np.errstate(divide="ignore", invalid="ignore"):
        x = focal_cells * (star_vectors @ right) / depth
        y = -focal_cells * (star_vectors @ up) / depth
    return x, y, depth


def green_cells(path):
    """Mean of the two green photosites per Bayer cell, black-level subtracted."""
    raw = rawpy.imread(path)
    image = raw.raw_image_visible.astype(np.float32) - float(raw.black_level_per_channel[0])
    pattern = raw.raw_pattern.tolist()
    saturated = raw.raw_image_visible.max() >= raw.white_level
    green_offsets = [(r, c) for r in (0, 1) for c in (0, 1)
                     if raw.color_desc[pattern[r][c]:pattern[r][c] + 1] == b"G"]
    planes = [image[r::2, c::2] for r, c in green_offsets]
    return (0.5 * (planes[0] + planes[1]), saturated,
            float(raw.white_level - raw.black_level_per_channel[0]))


def detect_sources(cells, block=48, threshold_sigma=8.0, limit=400):
    height, width = cells.shape
    rows, cols = height // block, width // block
    tiles = cells[:rows * block, :cols * block].reshape(rows, block, cols, block)
    coarse = np.median(tiles, axis=(1, 3))
    background = np.repeat(np.repeat(coarse, block, axis=0), block, axis=1)
    residual = cells[:rows * block, :cols * block] - background
    sigma = 1.4826 * np.median(np.abs(residual - np.median(residual)))

    peak = residual.copy()
    shifted = [np.roll(residual, (dy, dx), axis=(0, 1))
               for dy in (-1, 0, 1) for dx in (-1, 0, 1) if (dy, dx) != (0, 0)]
    is_peak = np.all([peak > s for s in shifted], axis=0) & (peak > threshold_sigma * sigma)
    ys, xs = np.nonzero(is_peak)
    strength = residual[ys, xs]
    order = np.argsort(strength)[::-1][:limit]
    return np.stack([xs[order], ys[order]], axis=1).astype(np.float64), residual, sigma


def solve_pointing(sources, star_vectors, focal_cells, centre, exif_data):
    tree = cKDTree(sources)
    best = None
    for handedness in (1.0, -1.0):
        for pitch_step, roll_step, az_step, tol in ((2.0, 2.0, 2.0, 45.0), (0.4, 0.4, 0.4, 14.0)):
            if best is None or best[1] != handedness:
                pitch_grid = exif_data["pitch"] + np.arange(-6, 6.01, pitch_step)
                roll_grid = exif_data["roll"] + np.arange(-6, 6.01, roll_step)
                az_grid = np.arange(0, 360, az_step)
            else:
                pitch_grid = best[2] + np.arange(-3, 3.01, pitch_step)
                roll_grid = best[3] + np.arange(-3, 3.01, roll_step)
                az_grid = best[4] + np.arange(-3, 3.01, az_step)
            for pitch in pitch_grid:
                for roll in roll_grid:
                    for az in az_grid:
                        x, y, depth = project(star_vectors, pitch, az, roll, focal_cells, handedness)
                        inside = depth > 0.3
                        if inside.sum() < 5:
                            continue
                        points = np.stack([x[inside] + centre[0], y[inside] + centre[1]], axis=1)
                        distance, _ = tree.query(points, distance_upper_bound=tol)
                        score = np.sum(np.isfinite(distance))
                        if best is None or score > best[0]:
                            best = (score, handedness, pitch, roll, az)
    return best


def measure(path, cone_deg=SKY_CONE_DEG):
    exif_data = exif(path)
    jd = julian_day(*exif_data["utc"])
    cells, saturated, saturation_level = green_cells(path)
    height, width = cells.shape
    centre = (width / 2.0, height / 2.0)

    pixel_mm = SENSOR_WIDTH_MM / SENSOR_WIDTH_PX
    arcsec_per_cell = 206265.0 * (2 * pixel_mm) / exif_data["focal_mm"]
    focal_cells = exif_data["focal_mm"] / (2 * pixel_mm)

    ra, dec, mag = load_catalog()
    alt, az = star_altaz(ra, dec, jd)
    bright = (mag < 5.6) & (alt > 2.0)
    vectors = unit_vectors(alt[bright], az[bright])

    sources, residual, sigma = detect_sources(cells)
    solution = solve_pointing(sources, vectors, focal_cells, centre, exif_data)
    score, handedness, pitch, roll, az_solved = solution

    # Photometry on identified stars, central field only: off-axis cos^4 vignetting
    # and the gnomonic scale change both stay under ~0.1 mag inside 15 degrees.
    x, y, depth = project(vectors, pitch, az_solved, roll, focal_cells, handedness)
    px, py = x + centre[0], y + centre[1]
    off_axis = np.arctan(np.hypot(x, y) / focal_cells)
    tree = cKDTree(sources)
    catalog_mag = mag[bright]
    catalog_alt = alt[bright]

    zero_points, saturated_stars = [], 0
    for i in range(len(px)):
        if depth[i] <= 0 or off_axis[i] > 25 * D2R:
            continue
        distance, index = tree.query([px[i], py[i]], distance_upper_bound=14.0)
        if not np.isfinite(distance):
            continue
        sx, sy = sources[index]
        radius = 7
        y0, y1 = int(sy - radius), int(sy + radius + 1)
        x0, x1 = int(sx - radius), int(sx + radius + 1)
        if y0 < 0 or x0 < 0 or y1 > height or x1 > width:
            continue
        stamp = cells[y0:y1, x0:x1]
        # A clipped core loses flux, which biases the zero point low and the sky
        # bright. The green cell is a 2-photosite mean, so clipping shows at half
        # the raw ceiling; stay well under it.
        if stamp.max() >= 0.75 * saturation_level:
            saturated_stars += 1
            continue
        local_background = np.median(cells[max(0, y0 - 20):y1 + 20, max(0, x0 - 20):x1 + 20])
        flux = float(np.sum(stamp - local_background))
        if flux <= 0:
            continue
        airmass = 1.0 / math.sin(max(catalog_alt[i], 5.0) * D2R)
        vignette = math.cos(off_axis[i]) ** 4
        zero_points.append(catalog_mag[i] + EXTINCTION_K * airmass
                           + 2.5 * math.log10(flux / vignette))

    zero_points = np.array(zero_points)
    keep = np.abs(zero_points - np.median(zero_points)) < 3 * 1.4826 * np.median(
        np.abs(zero_points - np.median(zero_points)))
    zero_point = float(np.median(zero_points[keep]))

    # Sky background: darkest quartile of block medians in the central field rejects
    # stars, the Milky Way, and any ground at the frame edge.
    block = 64
    rows, cols = height // block, width // block
    tiles = cells[:rows * block, :cols * block].reshape(rows, block, cols, block)
    block_median = np.median(tiles, axis=(1, 3))
    yy, xx = np.mgrid[0:rows, 0:cols]
    block_x = (xx + 0.5) * block - centre[0]
    block_y = (yy + 0.5) * block - centre[1]
    block_theta = np.arctan(np.hypot(block_x, block_y) / focal_cells)
    central = block_theta < cone_deg * D2R
    values = block_median[central] / np.cos(block_theta[central]) ** 4
    sky_adu = float(np.percentile(values, 25))
    sky_median = float(np.percentile(values, 50))

    surface_brightness = zero_point - 2.5 * math.log10(sky_adu / arcsec_per_cell ** 2)
    median_brightness = zero_point - 2.5 * math.log10(sky_median / arcsec_per_cell ** 2)

    return {
        "frame": path.split("/")[-1], "score": score, "stars_used": int(keep.sum()),
        "saturated_stars": saturated_stars,
        "zero_point": zero_point, "zp_scatter": float(np.std(zero_points[keep])),
        "alt": pitch, "az": az_solved, "roll": roll, "handedness": handedness,
        "exif_pitch": exif_data["pitch"], "saturated": saturated,
        "sky_adu": sky_adu, "mu_dark": surface_brightness, "mu_median": median_brightness,
        "arcsec_per_cell": arcsec_per_cell, "cone_deg": cone_deg,
    }


if __name__ == "__main__":
    arguments = sys.argv[1:]
    cone = SKY_CONE_DEG
    if "--site" not in arguments:
        raise SystemExit("--site LAT,LON is required: the Z5 has no GPS")
    if "--site" in arguments:
        at = arguments.index("--site")
        LAT, LON = (float(v) for v in arguments[at + 1].split(","))
        arguments = arguments[:at] + arguments[at + 2:]
    if "--cone" in arguments:
        at = arguments.index("--cone")
        cone = float(arguments[at + 1])
        arguments = arguments[:at] + arguments[at + 2:]
    for frame in arguments:
        result = measure(frame, cone)
        print("{frame}  alt {alt:5.1f} az {az:6.1f} roll {roll:6.1f}  "
              "matched {score:3d}  ZP {zero_point:6.2f} +/- {zp_scatter:.2f} ({stars_used} stars)  "
              "sat_skip {saturated_stars:2d}  cone {cone_deg:4.1f}  sky {sky_adu:7.1f} ADU  mu_dark {mu_dark:5.2f}  mu_median {mu_median:5.2f}".format(**result))
