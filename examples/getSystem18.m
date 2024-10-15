function [f, g, h] = getSystem18()
%getSystem18 Returns the 3DOF linear model from Bass & Webber 1966.
%
%   Usage:  [f,g,h] = getSystem18()
%
%   References: Bass & Webber 1966
%
%%

f = [0 1 0; 
    0 0 1; 
    -6 -11 -6];
g = [0;0;1];
h = g.';
end

