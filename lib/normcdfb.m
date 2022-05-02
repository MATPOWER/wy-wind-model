function pout = normcdfb(xin, mu, sigma)
% assumes normalized function, use of standard function erfc
% See
% https://www.mathworks.com/help/stats/normcdf.html
% probcumb = normcdfb(bin_bound);

% 2022.04.27
% Alberto J. Lamadrid L.

if nargin<3
	sigma = 1;
	if nargin<2
		mu = 0;
	end
end

znorm = (xin-mu) ./ sigma;

pout = 0.5 * erfc(-znorm ./ sqrt(2));