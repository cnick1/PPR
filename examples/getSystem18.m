function [f, g, h] = getSystem18()
%getSystem18 Returns the 3D linear model from [1].
%
%   Usage:  [f,g,h] = getSystem18()
%
%       \dot{x_1} = x2
%       \dot{x_2} = x3
%       \dot{x_3} = -6x1 - 11 x2 - 6 x3 + u
%
%   References: [1] R. Bass and R. Webber, "Optimal nonlinear feedback
%               control derived from quartic and higher-order performance
%               criteria,” IEEE Transactions on Automatic Control, vol. 11,
%               no. 3, pp. 448–454, Jul. 1966, doi: 10.1109/tac.1966.1098381.
%
%%

f = [0 1 0; 
    0 0 1; 
    -6 -11 -6];
g = [0;0;1];
h = g.';
end

