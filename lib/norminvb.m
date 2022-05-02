function xout = norminvb(pin)
% assumes normalized function, use of standard function erfcinv
% See 
% https://www.mathworks.com/help/stats/norminv.html
% bin_meanb = norminvb(probref);

% 2022.04.27
% Alberto J. Lamadrid L.

xout = -sqrt(2).*erfcinv(2*pin);

