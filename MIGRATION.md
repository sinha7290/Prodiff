# Migration notes

This repository previously carried copies of lab-internal scripts (`bone.py`,
`HegemonUtil.py`, `StepMiner.py`, `explore.conf`, `hegemonutils.pl`). Those were
roughly a megabyte of code that only ran against one filesystem, and they have
been removed.

The helper functions the notebooks actually used now come from
[bioutils](https://github.com/sinha7290/bioutils):

| was | now |
|---|---|
| `bone.readList` | `bu.read_list` |
| `bone.saveList` | `bu.save_list` |
| `bone.printOLS` | `bu.ols_table` |
| `bone.getCode` | `bu.significance_code` |
| `bone.getPDF` / `closePDF` | `bu.open_pdf` / `bu.close_pdf` |
| `hu.plotCoef` | `bu.plot_coefficients` |
| `hu.censor` | `bu.censor` |
| `hu.uniq` | `bu.unique` |
| `StepMiner.fitstep` | `bu.step_threshold` |

## Calls that still need attention

These have no drop-in equivalent, either because the signature differs or
because they reach into the Hegemon expression database, which is lab-internal
and not public:

- `hu.Multivariate(df)` - `bu.univariate_then_multivariate(data, outcome, predictors)`
  takes the outcome and predictor columns rather than a prefit frame.
- `hu.survival(time, status, groups)` - `bu.kaplan_meier(...)` returns
  `(ax, stats)` rather than just an axis; `stats` carries the log-rank p-value.
- `hu.getThrData(x)` - `bu.threshold_bounds(x)` returns `(threshold, low, high)`.
- `bone.adj_light(c, f, n)` - `bu.adjust_lightness(c, f)`; the third argument is gone.
- `bone.getEntries`, `bone.processGeneGroups`, `bone.processGeneGroupsDf`,
  `bone.BINetwork`, `bone.MacAnalysis` - all require the Hegemon database.

Note that `bu.step_threshold` is an independent implementation of single-step
StepMiner thresholding, so thresholds can differ marginally from the original
and figures will not be bit-identical.
