function pout = normcdf_nst(xin, mu, sigma)
% assumes normalized function, use of standard function erfc
% See
% https://www.mathworks.com/help/stats/normcdf.html
% probcumb = normcdf_nst(bin_bound);

% 2022.05.02
% Alberto J. Lamadrid L.

if nargin<3
    sigma = 1;
    if nargin<2
        mu = 0;
    end
end

znorm = (xin-mu) ./ sigma;

pout = 0.5 * erfc(-znorm ./ sqrt(2));
