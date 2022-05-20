WY-Wind-Model
=============

_by Wooyoung Jeon and Ray Zimmerman_

[WY-Wind-Model][1] provides a time-series wind model for use by the [MATPOWER
Optimal Scheduling Tool][2] (MOST), which is part of [MATPOWER][3].


System Requirements
-------------------
*   [MATLAB][4] version 9.1 (R2016b) or later, or
*   [GNU Octave][5] version 5.x or later
*   [MP-Test][6] version 7.1 or later


Installation
------------

Installation and use of WY-Wind-Model requires familiarity with the basic operation
of MATLAB or Octave, including setting up your MATLAB/Octave path.

1.  Clone the repository or download and extract the zip file of the WY-Wind-Model
    distribution from the [WY-Wind-Model project page][1] to the location of your
    choice.

2. Add the following directories to your MATLAB/Octave path:
    * `<WYWIND>/lib`
    * `<WYWIND>/lib/t`

    where `<WYWIND>` is used to denote the path to the `mp-element`
    directory you cloned or unzipped in step 1.

3.  At the MATLAB/Octave prompt, type `test_wy_wind` to run the test
    suite and verify that WY-Wind-Model is properly installed and functioning.
    The result should resemble the following:
```
  >> test_wy_wind
  t_wy_wind_date2pidx....ok
  t_wy_wind_pidx2date....ok
  t_wy_wind..............ok
  t_wy_wind_model........ok
  All tests successful (126 of 126)
  Elapsed time 0.09 seconds.
```


Usage
-----

*   Create a WY-Wind-Model object, with models for sites 2, 6, and 15 from
    the 16-site NPCC model:
    ```matlab
    widx = [2;6;15];
    wm = wy_wind_model('model_npcc', widx);
    ```

*   Generate transition probabilities for [MOST][2] for a 24 period horizon with
    4 wind scenarios in each period.
    ```matlab
    tp = wm.transition_probs(24, 4);
    ```

*   Load the original NPCC historical wind speed data (raw wind speeds in m/s)
    and convert it to `log(raw_wind_speed+1)`, since that is what this model
    is based on.
    ```matlab
    s = load('winddata_npcc');
    wind_data = s.wind_data;
    log_wind_data = log10(wind_data + 1);
    ```

*   Extract a 48-hour wind speed realization for the 3 sites of interest from
    the original NPCC historical data, beginning at noon on June 1, 2005 (data
    starts at 1am on Jan 1, 2004).
    ```matlab
    pidx0 = wy_wind_date2pidx(24, [2004 1 1 1 0 0], [2005 6 1 12 0 0]);
    wsr = wm.realizations(pidx0, 48, log_wind_data);
    ```

*   Generate a 24-hour wind speed realization from the model for the 3 sites
    of interest, beginning at the same hour.
    ```matlab
    wsr = wm.realizations(pidx0, 24);
    ```

*   Convert the wind speed realizations to wind power realizations, expressed
    as fractions of installed capacity, using the default (multi-turbine)
    power curve.
    ```matlab
    wpr = wm.speed2power(wsr);
    ```

*   Generate a 24-hour wind speed forecast from the model for the 3 sites
    of interest, beginning at the same hour, using 4 bins to represent the
    forecast distribution.
    ```matlab
    ws0 = log_wind_data(pidx0-1, widx);
    wsf = wm.forecasts(pidx0, ws0, 24, 4);
    ```

*   Convert the wind speed forecasts to wind power realizations, expressed
    as fractions of installed capacity, using the 4th power curve (off-shore)
    defined in `'WindPowerCurveIEC.txt'`.
    ```matlab
    s2p = wy_wind_power_curve(4, 'WindPowerCurveIEC.txt');
    wm = wy_wind_model('model_npcc', [2;6;15], s2p);
    wpf = wm.speed2power(wsf);
    ```

*   See `ex_wy_wind_simulation.m` for an example of using WY-Wind-Model
    to generate inputs for a receding horizon simulation.


Documentation
-------------

The primary sources of documentation for WY-Wind-Model are this section of
this README file and the built-in `help` command. As with the built-in
functions and toolbox routines in MATLAB and Octave, you can type `help`
followed by the name of a command or M-file to get help on that
particular function.

### Notation
- `np_all` — total number of periods in raw wind data
- `nw_all` — total number of wind sites in raw data (16 for `winddata_npcc.mat`)
- `np` — number of periods of interest (e.g. for planning horizon)
- `nw` — number of wind sites of interest
- `npd` — number of periods per day (typically 24, for hourly data)
- `nb` — number of bins used for wind model inputs for MOST
- `widx` — (`nw x 1`) vector of indices of wind sites of interest
- `pidx` — scalar period index for raw historical wind data
- `dt` — 1 x 6 standard Matlab date vector [_yr, mo, day, hr, min, sec_],
  specifying a specific period in the raw historical wind data,
  _yr_ = 4-digit year,
  _mo_ = month-of-year (1-12),
  _day_ = day-of-month (1-31),
  _hr_ = hour-of-day (0-23),
  _min_ = minute (0-59),
  _sec_ = second (0-59)
- `wind_data` - (`np_all x nw_all`) matrix of wind speeds (in m/s),
  corresponding to `np_all` periods, `nw_all` sites, a particular `npd`, and a
  starting date/time (`dt0`)
- `log_wind_data` — (`np_all x nw_all`) matrix equal to `log10(wind_data + 1)`
- `bins` —  bin specification, supplied as either:
  1. number of bins (`nb`), or
  2. (`1 x nb-1`) vector of bin boundaries (standard deviation coefficients),
     where initial `-Inf` and final `+Inf` are assumed, but not included
- `model` — struct with fields:
    - `type` — type of wind speed model
      - 0 = based on *raw_wind_speed* in m/s
      - 1 = based on log(*raw_wind_speed* + 1)
    - `ar1` — (`nw_all x 1`) vector of AR[1] coefficients for individual sites
    - `ols` — (`nw_all x 9`) matrix of OLS estimation parameters for individual
        sites:  
        [_C CY1 SY1 CY2 SY2 CD1 SD1 CD2 SD2_]
    - `var_wnr` — (`nw_all x nw_all`) covariance matrix for individual sites
    - `ar1_total` — scalar AR[1] coefficient for total wind
    - `ols_total` — `1 x 9` vector of OLS estimation parameters for total wind
    - `var_wnr_total` — scalar variance for total wind
- `ws` — wind speed quantity, units depend on the type of the model

### Data Files

- **winddata_npcc.mat** — raw historical wind speed data from NPCC
  - `wind_data` - (`26303 x 16`) matrix of wind speeds in m/s corresponding to:
    - 16 sites (ny1–ny9, ne1–ne7)
    - 3 years (2004, 2005, 2006) _(8760 * 3 + 24 -1 = 26303)_
      - beginning @ 2004-01-01 1:00 (1am)
      - ending at 2006-12-31 23:00 (11pm)

- **model_npcc.mat** — type 1 model for the 16 sites in NPCC
  - `model` - struct with model parameters for 16 sites based on
    **winddata_npcc.mat**

- **WindPowerCurveIEC.txt** — power curve lookup tables
  - col 1 is wind speeds in m/s from 0 to 30
  - cols 2–6 contain 5 power curves as follows, expressed as fraction of
    installed capacity:  
        1 = IEC1, 2 = IEC2, 3 = IEC3, 4 = offshore, 5 = multi-turbine


### Object Interface

- __wy_wind_model__ — constructor, creates a WY-Wind-Model object from a
  MAT-file containing a `model` struct, optionally selecting a subset of the
  sites indexed by `widx`, and/or providing a specific power curve for
  conversion from speed to power
    ```matlab
    wm = wy_wind_model(model_fname)
    wm = wy_wind_model(model_fname, widx)
    wm = wy_wind_model(model_fname, widx, s2p)
    ```

- __transition_probs__ — generates transition probabilities to use with
  [MOST][2] from the model
    ```matlab
    tp = wm.transition_probs(np, bins)
    ```

- __realizations__ — generates realizations from the model, or extracts them
  from the historical data, given the starting period and the number of
  periods desired
    ```matlab
    wsr = wm.realizations(pidx0, np)
    wsr = wm.realizations(pidx0, np, wind_data)
    ```

- __forecasts__ — generates forecasts of wind speed bin means from the model,
  for use with [MOST][2], given the starting period, initial wind speed, number
  of periods desired, and the bin specification
    ```matlab
     wsf = wm.forecasts(pidx0, ws0, np, bins)
    ```

- __speed2power__ — converts wind speed data to wind power, expressed as
  fraction of installed capacity
    ```matlab
    wp = wm.speed2power(ws)
    ```

- __display__ — called to display object on command line
    ```matlab
    wm
    ```

### Function Interface

- __wy_wind_trans_probs__ — generates transition probabilities for [MOST][2]
  from the model
    ```matlab
    tp = wy_wind_trans_probs(model, np, bins)
    ```

- __wy_wind_realizations__ — generates realizations from the model, or extracts
  them from the historical data, given the starting period and the number of
  periods desired
    ```matlab
    wsr = wy_wind_realizations(model, widx, pidx0, np)
    wsr = wy_wind_realizations(wind_data, widx, pidx0, np)
    wsr = wy_wind_realizations(log_wind_data, widx, pidx0, np)
    ```

- __wy_wind_forecasts__ — generates forecasts of wind speed bin means from the
  model, for use with [MOST][2], given the starting period, initial wind speed,
  number of periods desired, and the bin specification
    ```matlab
    wsf = wy_wind_forecasts(model, widx, pidx0, ws0, np, bins)
    ```

- __wy_wind_speed2power__ — converts wind speed data to wind power, expressed
  as fraction of installed capacity, based on provided power curve
    ```matlab
    wp = wy_wind_speed2power(ws, s2p)
    wp = wy_wind_speed2power(ws, idx)
    wp = wy_wind_speed2power(ws, s2p, ws_type)
    wp = wy_wind_speed2power(ws, idx, ws_type)
    ```

### Other Functions

#### Utility Functions

- __wy_wind_date2pidx__ — converts a date vector to a raw period index
    ```matlab
    pidx = wy_wind_date2pidx(npd, dt0, dt)
    ```

- __wy_wind_pidx2date__ — converts a raw period index to date vector
    ```matlab
    dt = wy_wind_pidx2date(npd, dt0, pidx)
    ```

- __wy_wind_power_curve_data__ — loads a power curve table from a file
    ```matlab
    s2p = wy_wind_power_curve_data()
    s2p = wy_wind_power_curve_data(idx)
    s2p = wy_wind_power_curve_data(idx, fname)
    ```

#### Statistics Toolbox Replacement Functions

The following functions are included in order to eliminate any dependence
on the Statistics Toolbox.

- __mvnrnd_nst__ — replacement for `mvnrnd()` based on `randn()`.
    ```matlab
    r = mvnrnd_nst(mu, sigma, n)
    ```

- __normcdf_nst__ — replacement for `normcdf()` based on `erfc()`.
    ```matlab
    pout = normcdf_nst(xin, mu, sigma)
    ```

- __norminv_nst__ — replacement for `norminv()` based on `erfcinv()`.
    ```matlab
    xout = norminv_nst(pin)
    ```

#### Other Function _(used internally)_

- __wy_wind_bins__ — returns number of bins and bin boundaries given a bin
  specification `bins`
    ```matlab
    [nb, bin_bounds] = wy_wind_bins(bins)
    ```


Publications
------------

1.  R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas,
    ["MATPOWER: Steady-State Operations, Planning and Analysis Tools
    for Power Systems Research and Education,"][7] *Power Systems, IEEE
    Transactions on*, vol. 26, no. 1, pp. 12–19, Feb. 2011.  
    doi: [10.1109/TPWRS.2010.2051168][7].


[Citing WY-Wind-Model][8]
--------------------------

We request that publications derived from the use of the WY-Wind-Model
explicitly acknowledge that fact by citing the following paper.

>   R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas, "MATPOWER:
    Steady-State Operations, Planning and Analysis Tools for Power Systems
    Research and Education," *Power Systems, IEEE Transactions on*, vol. 26,
    no. 1, pp. 12-19, Feb. 2011.  
    doi: [10.1109/TPWRS.2010.2051168][7]

License
-------

WY-Wind-Model is distributed under the [3-clause BSD license][9].


Acknowledgments
---------------

This material is based upon work supported in part by the National Science
Foundation under Grant No. ???????. Any opinions, findings, and
conclusions or recommendations expressed in this material are those of the
author(s) and do not necessarily reflect the views of the funding agencies.

----
[1]: https://github.com/MATPOWER/wy-wind-model
[2]: https://github.com/MATPOWER/most
[3]: https://github.com/MATPOWER/matpower
[4]: https://www.mathworks.com/
[5]: https://www.gnu.org/software/octave/
[6]: https://github.com/MATPOWER/mptest
[7]: https://doi.org/10.1109/TPWRS.2010.2051168
[8]: CITATION
[9]: LICENSE
