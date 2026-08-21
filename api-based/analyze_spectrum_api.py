"""Gamma-spectrum analyzer backed by the IAEA Live Chart API.

Runs the exact pipeline of ``gamma-spectrum-analyzer/analyze_spectrum.py`` --
same peak detection, fitting and activity calculation -- but identifies peaks
against a gamma library fetched from the public IAEA nuclear data service
instead of the curated ``gamma_library.txt``.

The library is built once (about a thousand small API requests) and cached in
``data/api_gamma_library.csv``; later runs read the cache.
"""

import importlib.util
import sys
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent

sys.path.insert(0, str(HERE))
import iaea_client


def load_analyzer():
    """Import the file-based analyzer module and reuse its pipeline."""
    path = HERE.parent / 'gamma-spectrum-analyzer' / 'analyze_spectrum.py'
    spec = importlib.util.spec_from_file_location('analyze_spectrum', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules['analyze_spectrum'] = module
    spec.loader.exec_module(module)
    return module


def main():
    sr = load_analyzer()
    settings = sr.load_settings()
    paths, detection = settings['paths'], settings['peak_detection']

    counts, meas_time, dead_time, acquired_at, calibration = sr.read_spe(
        sr.resolve(paths['spectrum_file']), paths.get('encoding', 'cp1251'))
    energies = sr.energy_axis(counts, calibration, settings)
    print(f'{len(counts)} channels spanning '
          f'{energies[0]:.2f}-{energies[-1]:.2f} keV')

    irradiated_at = datetime.strptime(settings['sample']['irradiation_date'],
                                      sr.DATE_FORMAT)
    elapsed = acquired_at - irradiated_at
    hours_after_irradiation = elapsed.total_seconds() / 3600

    counts = sr.moving_average(counts)
    filtered_spectrum = sr.convolution_filter(counts)
    search = sr.apply_thresholds(filtered_spectrum, energies,
                                 sr.parse_bands(detection['bands']),
                                 detection.getfloat('sensitivity'))
    peak_energies, _ = sr.detect_peaks(search, energies)
    print(f'{len(peak_energies)} peaks detected')

    # The one difference from the file-based analyzer: the library source.
    library = iaea_client.library_records(iaea_client.build_gamma_library(
        progress=lambda i, n: print(f'  building library {i}/{n} nuclides')))
    print(f'{len(library)} gamma lines from the IAEA API')
    peak_matches = sr.match_library(peak_energies, library,
                                    detection.getfloat('match_tolerance_kev'),
                                    hours_after_irradiation)
    sr.deduplicate(peak_matches)

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

    report_path = HERE / 'report_api.txt'
    with open(report_path, 'a', encoding='utf-8') as report:
        report.write('START OF REPORT (IAEA API library)\n')
        report.write(f'irradiation date: {irradiated_at}\n')
        report.write(f'aquisition date: {acquired_at}\n')
        for target_energy in peaks_of_interest:
            sr.analyse_peak(target_energy, peak_energies, peak_matches, counts,
                            energies, filtered_spectrum, meas_time, dead_time,
                            elapsed.total_seconds(), sample, report)
        report.write('END OF REPORT\n\n\n')

    sr.plot_spectrum(counts, energies, peak_energies, peak_matches, settings)


if __name__ == '__main__':
    main()
