function success = test_wy_wind(verbose, exit_on_fail)
%TEST_WY_WIND  Run all WY-Wind-Model tests.
%   TEST_WY_WIND
%   TEST_WY_WIND(VERBOSE)
%   TEST_WY_WIND(VERBOSE, EXIT_ON_FAIL)
%   SUCCESS = TEST_WY_WIND(...)
%
%   Runs all of the WY-Wind-Model tests. If VERBOSE is true (false by default),
%   it prints the details of the individual tests. If EXIT_ON_FAIL is true
%   (false by default), it will exit MATLAB or Octave with a status of 1
%   unless T_RUN_TESTS returns ALL_OK.
%
%   See also T_RUN_TESTS.

%   WY-Wind-Model
%   Copyright (c) 2022, Wooyoung Jeon, Ray Zimmerman
%   by Ray Zimmerman
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

if nargin < 2
    exit_on_fail = 0;
    if nargin < 1
        verbose = 0;
    end
end

tests = {};

tests{end+1} = 't_wy_wind_date2pidx';
tests{end+1} = 't_wy_wind_pidx2date';
tests{end+1} = 't_wy_wind';
tests{end+1} = 't_wy_wind_model';

%% run the tests
all_ok = t_run_tests( tests, verbose );

%% handle success/failure
if exit_on_fail && ~all_ok
    exit(1);
end
if nargout
    success = all_ok;
end
