function xout = norminv_nst(pin)
%NORMINV_NST  Replacement for NORMINV based on ERFCINV.
%   XOUT = NORMINV_NST(PIN)
%   Assumes normalized function, use of standard function ERFCINV.
%   Removes dependency on Statistics Toolbox.
%
%   See https://www.mathworks.com/help/stats/norminv.html

%   WY-Wind-Model
%   Copyright (c) 2022, Alberto J. Lamadrid L.
%   by Alberto J. Lamadrid L., 2022.04.27
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

xout = -sqrt(2).*erfcinv(2*pin);
