function r = mvnrnd_nst(mu, sigma, n)
%MVNRND_NST  Replacement for MVNRND based only on RANDN.
%   R = MVNRND_NST(MU, SIGMA, N)
%   Only handles the case where MU is a 1-by-D vector and SIGMA is a D-by-D
%   co-variance matrix.
%
%   See https://stackoverflow.com/a/14517624 for description of the math.

%   WY-Wind-Model
%   Copyright (c) 2022, Ray Zimmerman
%   by Ray Zimmerman
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

d = size(mu, 2);
r = repmat(mu, n, 1) + randn(n, d) * chol(sigma);
