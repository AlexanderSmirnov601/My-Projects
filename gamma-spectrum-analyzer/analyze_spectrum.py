"""Gamma-spectrum peak finder, isotope identifier and activity calculator.

Reads an Ortec/Maestro ``.Spe`` file, smooths and convolution-filters the
spectrum, locates peaks, matches them against a gamma library, fits Gaussians to
the peaks of interest and derives the sample activity and production yield.

Everything run- and sample-specific lives in ``settings.ini`` next to this file.
"""

import configparser
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.text import Text
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
DATE_FORMAT = '%m/%d/%Y %H:%M:%S'

# A plain "<mass number><element symbol>" library entry, optionally metastable.
ISOTOPE_NAME = re.compile(r'^(\d{1,3})([A-Za-z]{1,2})m?$')

# How far either side of a detected position to look for the true peak channel,
# and how far a fitting window may extend before it is cut off.
PEAK_ANCHOR_HALF_WIDTH = 3
PEAK_MAX_HALF_WINDOW = 20

# Most x-axis ticks to place on the summary spectrum plot.
MAX_XTICKS = 20


# --------------------------------------------------------------------------- #
# configuration
# --------------------------------------------------------------------------- #
def load_settings(path=None):
    """Read settings.ini. A ``settings.local.ini`` beside it takes precedence."""
    path = Path(path) if path else HERE / 'settings.ini'
    settings = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
    if not settings.read([path, path.with_suffix('.local.ini')],
                         encoding='utf-8'):
        raise FileNotFoundError(f'no settings file found at {path}')
    return settings


def resolve(path):
    """Interpret a configured path relative to this script's directory."""
    path = Path(path)
    return path if path.is_absolute() else HERE / path


def parse_bands(raw):
    """``"0:82:0  82:383:2.0"`` -> ``[(0.0, 82.0, 0.0), (82.0, 383.0, 2.0)]``.

    An empty upper bound means "to the end of the spectrum".
    """
    bands = []
    for chunk in raw.replace(',', ' ').split():
        low, high, factor = chunk.split(':')
        bands.append((float(low),
                      float(high) if high.strip() else np.inf,
                      float(factor)))
    return bands


# --------------------------------------------------------------------------- #
# file readers
# --------------------------------------------------------------------------- #
def _floats(line):
    """Leading numeric tokens of a line, stopping at the first non-number."""
    values = []
    for token in line.split():
        try:
            values.append(float(token))
        except ValueError:
            break                       # trailing unit such as "keV"
    return values


def read_spe(path, encoding):
    """Parse an Ortec ``.Spe`` file.

    Returns ``(counts, meas_time, dead_time, acquired_at, calibration)``, where
    ``calibration`` holds the $MCA_CAL (or $ENER_FIT) polynomial coefficients,
    lowest order first, or ``None`` when the file carries no calibration.
    $MEAS_TIM holds "<live time> <real time>"; the second value is used as the
    dead time.
    """
    with open(path, 'r', encoding=encoding, errors='replace') as handle:
        lines = handle.read().splitlines()

    counts, meas_time, dead_time = [], None, None
    acquired_at, calibration, ener_fit = None, None, None

    index = 0
    while index < len(lines):
        marker = lines[index].strip()
        if marker == '$MEAS_TIM:':
            meas_time, dead_time = _floats(lines[index + 1])[:2]
            index += 2
        elif marker == '$DATE_MEA:':
            acquired_at = datetime.strptime(lines[index + 1].strip(),
                                            DATE_FORMAT)
            index += 2
        elif marker == '$DATA:':
            first, last = (int(value) for value in lines[index + 1].split()[:2])
            channels = last - first + 1
            counts = [float(value)
                      for value in lines[index + 2:index + 2 + channels]]
            index += 2 + channels
        elif marker == '$MCA_CAL:':
            calibration = _floats(lines[index + 2])
            index += 3
        elif marker == '$ENER_FIT:':
            ener_fit = _floats(lines[index + 1])
            index += 2
        else:
            index += 1

    if not counts:
        raise ValueError(f'no $DATA: block found in {path}')
    if acquired_at is None or meas_time is None:
        raise ValueError(f'{path} is missing $DATE_MEA: or $MEAS_TIM:')

    return counts, meas_time, dead_time, acquired_at, calibration or ener_fit


def load_gamma_library(path, max_mass_number):
    """Read ``energy;intensity;isotope;half_life;unit`` records into a list."""
    library = []
    with open(path, 'r', encoding='utf-8', errors='replace') as handle:
        for line in handle:
            fields = line.strip().split(';')
            if len(fields) != 5:
                continue
            name = ISOTOPE_NAME.match(fields[2])
            if not name or int(name.group(1)) > max_mass_number:
                continue
            library.append([float(fields[0]), float(fields[1]), fields[2],
                            float(fields[3]), fields[4].strip()])
    if not library:
        raise ValueError(f'no usable records in {path}')
    return library


# --------------------------------------------------------------------------- #
# curve fitting
# --------------------------------------------------------------------------- #
def gauss(x, H, A, x0, sigma):
    return H + np.abs(A) * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def gauss2(x, H, H1, A, A1, x0, x01, sigma, sigma1):
    return (H + np.abs(A) * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
            + H1 + np.abs(A1) * np.exp(-(x - x01) ** 2 / (2 * sigma1 ** 2)))


def gauss_fit(x, y):
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, _ = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma],
                        maxfev=50000)
    return popt


def gauss_fit2(x, y, x1, y1, x0, y0):
    mean = sum(x * y) / sum(y)
    mean1 = sum(x1 * y1) / sum(y1)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    sigma1 = np.sqrt(sum(y1 * (x1 - mean) ** 2) / sum(y1))
    popt, _ = curve_fit(gauss2, x0, y0,
                        p0=[min(y), min(y1), max(y), max(y1),
                            mean, mean1, sigma, sigma1],
                        maxfev=50000)
    return popt


# --------------------------------------------------------------------------- #
# spectrum processing
# --------------------------------------------------------------------------- #
def moving_average(values):
    """Symmetric 3-point moving average; the two endpoints are left untouched."""
    values = np.asarray(values, dtype=float)
    smoothed = values.copy()
    smoothed[1:-1] = (values[:-2] + values[1:-1] + values[2:]) / 3
    return smoothed.tolist()


def convolution_filter(counts):
    """Background suppression: narrow positive core, wide negative wings."""
    c = np.asarray(counts, dtype=float)
    filtered = c.copy()
    lo, hi = 5, len(c) - 5
    if hi <= lo:
        return filtered
    filtered[lo:hi] = (
        -(c[lo - 3:hi - 3] + c[lo + 3:hi + 3]
          + c[lo - 4:hi - 4] + c[lo + 4:hi + 4]) / 4
        + (c[lo - 1:hi - 1] + c[lo:hi] + c[lo - 2:hi - 2]
           + c[lo + 2:hi + 2] + c[lo + 1:hi + 1]) / 5
    )
    return filtered


def apply_thresholds(filtered, energies, bands, sensitivity):
    """Zero everything that is not a credible peak, band by band."""
    search = np.asarray(filtered, dtype=float).copy()

    keep = np.zeros(len(search), dtype=bool)
    for low, high, factor in bands:
        if factor:
            keep |= (energies >= low) & (energies < high)
    search[~keep] = 0.0

    reference_max = search.max() if search.size else 0.0
    for low, high, factor in bands:
        if not factor:
            continue
        in_band = (energies >= low) & (energies < high)
        search[in_band & (search < reference_max / sensitivity * factor)] = 0.0
    return search


def detect_peaks(search, energies):
    """Energy and height of the tallest channel in each nonzero run."""
    peak_energies, peak_heights = [], []
    best_height, best_energy, inside_peak = 0.0, 0.0, False

    for index, value in enumerate(search):
        if value > best_height:
            inside_peak = True
            best_height, best_energy = value, energies[index]
        if value == 0 and inside_peak:
            peak_energies.append(best_energy)
            peak_heights.append(best_height)
            best_height, best_energy, inside_peak = 0.0, 0.0, False

    return peak_energies, peak_heights


def energy_axis(counts, calibration, settings):
    """Build the keV axis, preferring the file's own calibration."""
    channels = np.arange(len(counts), dtype=float)
    mode = settings['spectrum'].get('calibration', 'file').strip().lower()

    if mode == 'file' and calibration:
        return sum(coefficient * channels ** order
                   for order, coefficient in enumerate(calibration))
    if mode == 'file':
        print('warning: no $MCA_CAL/$ENER_FIT in the spectrum file, '
              'falling back to the linear energy_max_kev axis')
    return np.linspace(0.0, settings['spectrum'].getfloat('energy_max_kev'),
                       len(counts))


def index_of_energy(energies, energy):
    """Channel index whose energy equals ``energy`` exactly."""
    matches = np.flatnonzero(np.asarray(energies) == energy)
    if not matches.size:
        raise ValueError(f'energy {energy} is not on the spectrum axis')
    return int(matches[-1])


# --------------------------------------------------------------------------- #
# library matching
# --------------------------------------------------------------------------- #
HOURS_PER_UNIT = {'h': 1.0, 'd': 24.0, 'y': 24.0 * 365.0}


def decay_criterion(record, hours_after_irradiation):
    """Rank a candidate line by how much of it should still be present.

    Emission intensity per hour of half-life, decayed forward to the moment of
    acquisition. Minute- and second-lived lines are pushed to the bottom.
    """
    intensity, half_life, unit = record[1], record[3], record[4]
    if intensity > 100:
        intensity /= 10
    if unit in HOURS_PER_UNIT:
        half_life_hours = half_life * HOURS_PER_UNIT[unit]
        return (intensity / half_life_hours
                * 2 ** (-hours_after_irradiation / half_life_hours))
    if unit == 'm':
        return -1.0
    if unit == 's':
        return -10.0
    return float('-inf')


def match_library(peak_energies, library, tolerance, hours_after_irradiation):
    """For each detected peak, the one or two most plausible library lines."""
    matches = []
    for energy in peak_energies:
        candidates = [record for record in library
                      if abs(record[0] - energy) <= tolerance]
        if not candidates:
            matches.append(['null'])
        elif len(candidates) == 1:
            matches.append([candidates[0]])
        else:
            ranked = sorted(
                candidates,
                key=lambda record: -decay_criterion(record,
                                                    hours_after_irradiation))
            matches.append([ranked[0], ranked[1]])
    return matches


def deduplicate(peak_matches):
    """Trim ambiguous matches without overriding the ranking.

    Within one peak, the second candidate is dropped when it repeats the
    first's isotope or line energy. Across peaks, when two peaks share a
    primary isotope -- or a primary mass number, since members of one mass
    chain coexist in a sample -- the shared primary already explains both
    peaks and their secondary candidates are dropped. A peak's own
    best-ranked candidate is never displaced.
    """
    for k in range(len(peak_matches)):
        if len(peak_matches[k]) == 2:
            if (peak_matches[k][0][2] == peak_matches[k][1][2]
                    or peak_matches[k][0][0] == peak_matches[k][1][0]):
                del peak_matches[k][1]

        for i in range(k + 1, len(peak_matches)):
            if (len(peak_matches[i]) == 2
                    and peak_matches[k][0][2] == peak_matches[i][0][2]):
                del peak_matches[i][1]
                if len(peak_matches[k]) == 2:
                    del peak_matches[k][1]
            if (len(peak_matches[i]) == 2
                    and peak_matches[k][0][2][0:3] == peak_matches[i][0][2][0:3]):
                del peak_matches[i][1]
                if len(peak_matches[k]) == 2:
                    del peak_matches[k][1]


# --------------------------------------------------------------------------- #
# peak area extraction
# --------------------------------------------------------------------------- #
def interpolate(energies, counts):
    """Straight-line background through the three channels at each edge."""
    nodes, node_energies, slope = [], [], 0.0
    for i in range(3):
        nodes.append(counts[i])
        node_energies.append(energies[i])
    for i in range(3):
        nodes.append(counts[-3 + i])
        node_energies.append(energies[-3 + i])
    for i in range(3):
        slope += ((nodes[i + 3] - nodes[i])
                  / (node_energies[i + 3] - node_energies[i]) / 3)
    interpolation = [slope * energy + (nodes[0] - slope * node_energies[0])
                     for energy in energies]
    return interpolation, sum(interpolation[2:-3])


def _finish_cut(energies, counts, filtered_spectrum):
    """Package a windowed peak and subtract its interpolated background."""
    if len(counts) < 6:
        raise ValueError(
            f'peak window is only {len(counts)} channels wide; interpolate() '
            'needs at least 6 to fit a background. Check the peak_detection '
            'thresholds in settings.ini.')
    result = [energies, counts, filtered_spectrum]
    interpolation, background_area = interpolate(result[0], result[1])
    result.append(interpolation)

    net = [result[1][i] - result[3][i] for i in range(len(result[3]))]
    minimum = min(net)
    if minimum < 0:
        net = [value - minimum for value in net]
    result[3] = net
    return result, background_area


def area_cutting(peak_energies, filtered_spectrum, counts, energies, k,
                 sensitivity, alignment):
    """Cut a window around peak ``k-1``, absorbing neighbours within ``sensitivity``."""
    flag1 = flag2 = 0
    boundary1 = boundary2 = 0
    number = 1
    left_or_right = 0

    if (k - 2 >= 0
            and np.abs(peak_energies[k - 1] - peak_energies[k - 2]) <= sensitivity):
        x1 = index_of_energy(energies, peak_energies[k - 2])
        left_or_right = -1
        number += 1
    else:
        x1 = index_of_energy(energies, peak_energies[k - 1])

    if (k + 1 <= len(peak_energies)
            and np.abs(peak_energies[k - 1] - peak_energies[k]) <= sensitivity):
        x2 = index_of_energy(energies, peak_energies[k])
        number += 1
        left_or_right = 1
    else:
        x2 = index_of_energy(energies, peak_energies[k - 1])

    step = min(min(filtered_spectrum[x1 - 10:x2 + 10]) / 15, -3)

    for i in range(100):
        if filtered_spectrum[x2 + i] >= step and flag1 == 1:
            if number == 2:
                boundary2 = i + 8
            elif alignment == 1:
                x3 = index_of_energy(energies, peak_energies[k])
                boundary2 = np.abs((x3 - x1) // 2) - 2
            else:
                boundary2 = i
            flag1 = 2
            break
        elif filtered_spectrum[x2 + i] < step and flag1 == 0:
            flag1 = 1

    for i in range(100):
        if filtered_spectrum[x1 - i] >= step and flag2 == 1:
            if number == 2:
                boundary1 = i + 5
            elif alignment == -1:
                x3 = index_of_energy(energies, peak_energies[k - 2])
                boundary1 = np.abs((x3 - x1) // 2) - 2
            else:
                boundary1 = i - 3
            flag2 = 2
            break
        elif filtered_spectrum[x1 - i] < step and flag2 == 0:
            flag2 = 1

    window = slice(x1 - boundary1, x2 + boundary2)
    result, background_area = _finish_cut(energies[window], counts[window],
                                          filtered_spectrum[window])
    return result, number, left_or_right, background_area


def cut_single_peak(peak_energies, filtered_spectrum, counts, energies, k):
    """Cut a window around an isolated peak, stopping where the slope turns."""
    x = index_of_energy(energies, peak_energies[k - 1])

    # Anchor on the tallest channel near the detected position; the slice is
    # clamped so a peak at the spectrum edge cannot wrap around.
    low = max(x - PEAK_ANCHOR_HALF_WIDTH, 0)
    high = min(x + PEAK_ANCHOR_HALF_WIDTH, len(counts))
    neighbourhood = counts[low:high]
    x = low + neighbourhood.index(max(neighbourhood))
    x = min(max(x, 1), len(counts) - 2)     # keep x-1 and x+1 addressable

    boundary1 = boundary2 = PEAK_MAX_HALF_WINDOW
    slope = (counts[x - 1] - counts[x]) / (energies[x - 1] - energies[x])
    offset = counts[x] - slope * energies[x]
    for i in range(2, min(PEAK_MAX_HALF_WINDOW, x + 1)):
        if counts[x - i] > (slope * energies[x - i] + offset):
            boundary1 = i + 1
            break

    slope = (counts[x] - counts[x + 1]) / (energies[x] - energies[x + 1])
    offset = counts[x] - slope * energies[x]
    for i in range(2, min(PEAK_MAX_HALF_WINDOW, len(counts) - x)):
        if counts[x + i] > (slope * energies[x + i] + offset):
            boundary2 = i + 1
            break

    window = slice(max(x - boundary1, 0), min(x + boundary2, len(counts)))
    result, _ = _finish_cut(energies[window], counts[window],
                            filtered_spectrum[window])
    return result


# --------------------------------------------------------------------------- #
# physics
# --------------------------------------------------------------------------- #
def activity_calculation(half_life_time, mass_attenuation_coeff,
                         sample_half_height, gamma_yield, count_time,
                         dead_time, cps, count_efficiency,
                         time_after_irradiation, irradiation_time,
                         thousands_of_particles):
    """Saturation activity and production yield from a measured count rate."""
    decay_constant = np.log(2) / half_life_time
    absorbtion_coeff = np.exp(-mass_attenuation_coeff * sample_half_height)
    decay_coeff = float(
        decay_constant * count_time
        / (1 - np.exp(-decay_constant * count_time))
        * (np.exp(decay_constant * dead_time) - 1)
        / (decay_constant * dead_time))
    activity = float(
        cps / (count_efficiency * absorbtion_coeff * gamma_yield)
        * decay_coeff * np.exp(decay_constant * time_after_irradiation))
    element_yield = float(
        activity / (1 - np.exp(-decay_constant * irradiation_time))
        * (1 - np.exp(-decay_constant * 3600)) * 360 / thousands_of_particles)
    return activity, element_yield


def _report_peak(net_counts, constant, meas_time, sample, dead_time,
                 elapsed_seconds, report):
    """Print and log the count rate, activity and yield for one fitted peak."""
    count_rate = sum(net_counts) / meas_time
    error = np.sqrt(sum(net_counts) * (1 + constant)) / meas_time

    activity, element_yield = activity_calculation(
        sample['half_life_time'], sample['mass_attenuation_coeff'],
        sample['sample_half_height'], sample['gamma_yield'], meas_time,
        dead_time, count_rate, sample['count_efficiency'], elapsed_seconds,
        sample['irradiation_time'], sample['thousands_of_particles'])

    print('Count rate w/ bg noise = ', round(count_rate, 2), 'counts/s')
    print('standart error = ', round(error, 2), 'counts/s')
    print('activity = ', round(activity, 2), 'decay/s')
    print('yield = ', round(element_yield, 2), 'Bq/mkA*h')
    report.write(f'{round(count_rate, 2)}    {round(error, 2)}\n')


def analyse_peak(target_energy, peak_energies, peak_matches, counts, energies,
                 filtered_spectrum, meas_time, dead_time, elapsed_seconds,
                 sample, report):
    """Fit the peak nearest ``target_energy`` and report activity from it."""
    distances = [np.abs(energy - target_energy) for energy in peak_energies]
    peak_index = distances.index(min(distances))
    print(f'--- {target_energy} keV -> peak at '
          f'{peak_energies[peak_index]:.2f} keV: {peak_matches[peak_index]}')

    result, peaks_crossed, side, _ = area_cutting(
        peak_energies, filtered_spectrum, counts, energies, peak_index + 1, 5, 0)

    if peaks_crossed == 1:
        result = cut_single_peak(peak_energies, filtered_spectrum, counts,
                                 energies, peak_index + 1)
        _, amplitude, centre, sigma = gauss_fit(result[0], result[1])
        net = gauss(result[0], 0, amplitude, centre, sigma)
        constant = (4 * (1 + (2 * len(result[0]) + 1) / 6)
                    * (max(result[1]) - max(result[3])) / max(result[1]))
        _report_peak(net, constant, meas_time, sample, dead_time,
                     elapsed_seconds, report)

        plt.figure(figsize=(14, 12))
        plt.plot(result[0], result[1], 'red')
        plt.plot(result[0], result[2], 'purple')
        plt.plot(result[0], gauss(result[0], *gauss_fit(result[0], result[1])),
                 '--b', label='fit')
        plt.grid()
        plt.show()

    elif peaks_crossed == 2:
        result0, _, _, _ = area_cutting(peak_energies, filtered_spectrum,
                                        counts, energies, peak_index + 1, 1,
                                        side)
        if side == -1:
            result1, _, _, background_area = area_cutting(
                peak_energies, filtered_spectrum, counts, energies,
                peak_index, 1, 1)
        else:
            result1, _, _, background_area = area_cutting(
                peak_energies, filtered_spectrum, counts, energies,
                peak_index + 2, 1, -1)

        parameters = gauss_fit2(result0[0], result0[1], result1[0], result1[1],
                                result[0], result[1])
        net = gauss(result[0], 0, np.abs(parameters[2]), parameters[4],
                    parameters[6])
        constant = (4 * (1 + (2 * len(result0[0]) + 1) / 6)
                    * background_area / sum(net))
        _report_peak(net, constant, meas_time, sample, dead_time,
                     elapsed_seconds, report)

        plt.figure(figsize=(10, 10))
        plt.plot(result[0], result[1], 'green')
        plt.plot(result[0], gauss2(result[0], *parameters), '--b', label='fit')
        plt.plot(result1[0], result1[1], 'orange')
        plt.plot(result0[0], result0[1], 'red')
        plt.plot(result[0], result[2], 'purple')
        plt.grid()
        plt.show()

    else:
        print(f'{peaks_crossed} overlapping peaks is not supported, skipping')


# --------------------------------------------------------------------------- #
# summary plot
# --------------------------------------------------------------------------- #
def tick_step(span):
    """A round tick spacing that keeps at most MAX_XTICKS ticks in the window."""
    for step in (0.5, 1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000):
        if span / step <= MAX_XTICKS:
            return step
    return 2000


def _overlaps(a, b, pad):
    return (a.x0 < b.x1 + pad and a.x1 > b.x0 - pad
            and a.y0 < b.y1 + pad and a.y1 > b.y0 - pad)


def _text_extent(annotation, renderer):
    """Pixel box of the annotation's text alone.

    Annotation.get_window_extent() unions the text with the connector arrow,
    which stays anchored at the peak, so collision tests must measure the text
    by itself.
    """
    annotation.update_positions(renderer)
    return Text.get_window_extent(annotation)


def spread_labels(ax, annotations, pad=2.0):
    """Lift peak labels until none of them overlap.

    Each label is measured on the renderer and pushed up until it clears every
    label already placed, working left to right; the connector line back to its
    peak grows to match. Labels near the window edge are first shifted inward
    so they stay inside the axes.
    """
    if not annotations:
        return

    figure = ax.get_figure()
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()

    pixels_per_unit = (ax.transData.transform((0, 1))[1]
                       - ax.transData.transform((0, 0))[1])
    pixels_per_x_unit = (ax.transData.transform((1, 0))[0]
                         - ax.transData.transform((0, 0))[0])
    if pixels_per_unit <= 0 or pixels_per_x_unit <= 0:
        return

    # Keep every label inside the axes horizontally: a peak at the very edge
    # of the window otherwise draws its label over the y-axis. The connector
    # still points at the true peak position, only the text shifts inward.
    axes_box = ax.get_window_extent(renderer)
    for annotation in annotations:
        box = _text_extent(annotation, renderer)
        shift_px = 0.0
        if box.x0 < axes_box.x0 + pad:
            shift_px = (axes_box.x0 + pad) - box.x0
        elif box.x1 > axes_box.x1 - pad:
            shift_px = (axes_box.x1 - pad) - box.x1
        if shift_px:
            x, y = annotation.xyann
            annotation.xyann = (x + shift_px / pixels_per_x_unit, y)

    placed = []
    for annotation in sorted(annotations, key=lambda a: a.xyann[0]):
        box = _text_extent(annotation, renderer)
        for _ in range(len(annotations) + 5):
            clash = next((other for other in placed
                          if _overlaps(box, other, pad)), None)
            if clash is None:
                break
            # An Annotation positions its text from xyann on every draw, so
            # xyann is the value to move.
            x, y = annotation.xyann
            annotation.xyann = (
                x, y + ((clash.y1 + pad) - box.y0) / pixels_per_unit)
            box = _text_extent(annotation, renderer)
        placed.append(box)

    # Tighten the headroom back down onto the tallest label.
    highest = max(box.y1 for box in placed) + 4 * pad
    ax.set_ylim(top=ax.transData.inverted().transform((0, highest))[1])


def plot_spectrum(counts, energies, peak_energies, peak_matches, settings):
    """Bar chart of the spectrum window with identified peaks annotated."""
    low = settings['plot'].getfloat('energy_min_kev')
    high = settings['plot'].getfloat('energy_max_kev')
    label_threshold = settings['plot'].getfloat('label_threshold')

    inside = np.flatnonzero((energies >= low) & (energies <= high))
    if not inside.size:
        print('warning: the plot energy window contains no channels')
        return
    window = slice(int(inside[0]), int(inside[-1]) + 1)
    counts, energies = counts[window], energies[window]

    figure, ax = plt.subplots(figsize=(16, 9))
    width = (high - low) / len(energies) * 1.6
    ax.bar(energies, counts, width=width)

    step = tick_step(high - low)
    ax.set_xticks(np.arange(low, high + step / 2, step))
    ax.set_xlim(low, high)
    # Headroom for the labels; spread_labels() trims it back afterwards.
    ax.set_ylim(bottom=0, top=max(counts) * 2.5)

    peak_lookup = {energy: i for i, energy in enumerate(peak_energies)}
    annotations = []
    for k, energy in enumerate(energies):
        i = peak_lookup.get(energy)
        if i is None or counts[k] <= max(counts) / label_threshold:
            continue
        match = peak_matches[i]
        if len(match) == 2:
            label = f'{match[0][2]}/{match[1][2]}?'
        elif match != ['null']:
            label = match[0][2]
        else:
            continue
        annotations.append(ax.annotate(
            label, xy=(energy, counts[k] + max(counts) / 45),
            xytext=(energy, counts[k] + max(counts) / 20),
            arrowprops=dict(arrowstyle='-', connectionstyle='arc3',
                            linewidth=0.8, color='0.4'),
            fontsize=12, rotation='vertical', ha='center',
            va='bottom', fontweight='bold'))

    spread_labels(ax, annotations)

    ax.set_xlabel('Energy of Gamma quants, keV', fontsize=12)
    ax.set_ylabel('counts', fontsize=12)
    ax.tick_params(axis='x', which='both', bottom=True, top=False,
                   labelbottom=True, labelsize=12)
    ax.tick_params(axis='y', which='both', right=False, left=True,
                   labelleft=True, labelsize=12)
    for position in ('right', 'top'):
        ax.spines[position].set_visible(False)
    plt.show()


# --------------------------------------------------------------------------- #
# entry point
# --------------------------------------------------------------------------- #
def main():
    settings = load_settings()
    paths, detection = settings['paths'], settings['peak_detection']

    counts, meas_time, dead_time, acquired_at, calibration = read_spe(
        resolve(paths['spectrum_file']), paths.get('encoding', 'cp1251'))
    energies = energy_axis(counts, calibration, settings)
    print(f'{len(counts)} channels spanning '
          f'{energies[0]:.2f}-{energies[-1]:.2f} keV')

    irradiated_at = datetime.strptime(settings['sample']['irradiation_date'],
                                      DATE_FORMAT)
    elapsed = acquired_at - irradiated_at
    hours_after_irradiation = elapsed.total_seconds() / 3600
    print(f'acquired {acquired_at}, {round(hours_after_irradiation, 2)} hours '
          f'/ {round(hours_after_irradiation / 24, 2)} days after irradiation')

    counts = moving_average(counts)
    filtered_spectrum = convolution_filter(counts)
    search = apply_thresholds(filtered_spectrum, energies,
                              parse_bands(detection['bands']),
                              detection.getfloat('sensitivity'))
    peak_energies, _ = detect_peaks(search, energies)
    print(f'{len(peak_energies)} peaks detected')

    library = load_gamma_library(resolve(paths['gamma_library']),
                                detection.getint('max_mass_number'))
    peak_matches = match_library(peak_energies, library,
                                 detection.getfloat('match_tolerance_kev'),
                                 hours_after_irradiation)
    deduplicate(peak_matches)

    sample = {
        'half_life_time': settings['sample'].getfloat('half_life_days') * 24 * 3600,
        'irradiation_time': settings['sample'].getfloat('irradiation_time_min') * 60,
        'gamma_yield': settings['sample'].getfloat('gamma_yield'),
        'count_efficiency': settings['sample'].getfloat('count_efficiency'),
        'mass_attenuation_coeff': settings['sample'].getfloat('mass_attenuation_coeff'),
        'sample_half_height': settings['sample'].getfloat('sample_half_height'),
        'thousands_of_particles': settings['sample'].getfloat('thousands_of_particles'),
    }
    peaks_of_interest = [float(value) for value
                         in settings['analysis']['peaks_of_interest']
                         .replace(',', ' ').split()]

    with open(resolve(paths['report_file']), 'a', encoding='utf-8') as report:
        report.write('START OF REPORT\n')
        report.write(f'irradiation date: {irradiated_at}\n')
        report.write(f'aquisition date: {acquired_at}\n')
        for target_energy in peaks_of_interest:
            analyse_peak(target_energy, peak_energies, peak_matches, counts,
                         energies, filtered_spectrum, meas_time, dead_time,
                         elapsed.total_seconds(), sample, report)
        report.write('END OF REPORT\n\n\n')

    plot_spectrum(counts, energies, peak_energies, peak_matches, settings)


if __name__ == '__main__':
    main()
