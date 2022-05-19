function pout = normcdf_nst(xin, mu, sigma)
%NORMCDF_NST  Replacement for NORMCDF based on ERFC.
%   POUT = NORMCDF_NST(XIN, MU, SIGMA)
%   Assumes normalized function, use of standard function ERFC.
%   Removes dependency on Statistics Toolbox.
%
%   See https://www.mathworks.com/help/stats/normcdf.html

%   WY-Wind-Model
%   Copyright (c) 2022, Alberto J. Lamadrid L.
%   by Alberto J. Lamadrid L., 2022.05.02
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

if nargin<3
    sigma = 1;
    if nargin<2
        mu = 0;
    end
end

znorm = (xin-mu) ./ sigma;

pout = 0.5 * erfc(-znorm ./ sqrt(2));
