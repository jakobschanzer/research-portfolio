# Validation and quality controls

## Automated checks

`make check` performs three compilations with warnings enabled:

- JANAL plus the six-plane volume backend
- hydrogen-bond analysis
- β-CD ring-geometry analysis

It also runs a publication-boundary check that permits only the reviewed processed/derived input files and rejects local private paths or raw trajectory data.

`make test` exercises the Python hydrogen-bond compiler with synthetic fixtures. The test verifies deterministic bootstrap output under a fixed seed, expected summary dimensions, a known mean, and a clear error for missing temperature files. Synthetic values exist only inside the test and never replace scientific results.

## Scientific checks

- Each plot reads data from a named, versioned workbook rather than hard-coded curves.
- The spatial-map bulk hydrogen-bond references are three documented scalars extracted from the corresponding 25 °C production workbooks.
- Occupancy-transition widths are written to a derived audit workbook with definitions and uncertainty assumptions.
- OA occupancy alignment is exported as an audit sheet so plotted coordinates can be compared with tabulated means.
- Figure 6 uses one explicit coordinate convention: negative `z` is left and positive `z` is right. CB7 and OA use the positive-`z` approach; β-CD shows its negative- and positive-portal branches on their corresponding sides.
- Wet-basin quadratic fits and enhanced low-occupancy points in Figure 6 are exported to `data/derived/pn_fattail_fits.xlsx` for numerical audit.
- Generated plots use redundant markers/line styles and grayscale-compatible encodings.

## Known limitations

- Raw trajectories and production run directories are not distributed, so this repository reproduces analysis from compact processed inputs rather than rerunning molecular simulations.
- The checked-in volume backend has no matching production `vol` input in the archive used to curate this repository; it is compiled and documented, but a host assignment is not guessed.
- Cross-bin covariance was unavailable for drying-transition uncertainty propagation.
- XLSX containers record creation timestamps, and the GR plotting backend assigns generated clip-path IDs to SVG exports. A repeat run can therefore change those file bytes even when workbook cell contents and rendered figures are identical; scientific validation compares the numerical and visual content rather than container hashes alone.
- The manuscript is under peer review; numerical conclusions remain subject to author review.
