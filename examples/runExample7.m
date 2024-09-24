function runExample7(exportData)
%runExample7 Runs 3D aircraft stall model
%   Usage:  [v, w] = runExample7()
%
%   runExample7() runs the aircraft stall stabilization model from
%   Garrard 1977 [1].
%
%   Outputs:
%       v,w - coefficients of the past and future energy function
%             approximations, respectively.
%
%   Reference: [1] W. L. Garrard and J. M. Jordan, “Design of nonlinear
%               automatic flight control systems,” Automatica, vol. 13,
%               no. 5, pp. 497–505, Sep. 1977,
%               doi: 10.1016/0005-1098(77)90070-x
%
%   Part of the PPR repository.
%%
fprintf('Running Example 7\n')
if nargin < 1, exportData = false; end
options = odeset('Events', @myEvent); % Kill sim if the solution blows up

%% Define original dynamics and controllers
[f, g, ~] = getSystem7();
F = @(x) kronPolyEval(f, x);
G = @(x) (g{1} + kronPolyEval(g(2:end), x));

Q = 0.25; R = 1; 

%% Plot x1 and control u
%  Compute the polynomial approximations to the future energy function
degree = 8;
[~,Gains] = ppr(f, g, Q, R, degree);

tspan = [0, 12];

for alpha0 = [25, 27, 30, 35]
    X0 = [pi / 180 * alpha0; 0; 0]; % Initial condition, convert to radians
    figure('Position', [415 47 793.3333 347.3333]);
    subplot(1, 2, 1); hold on; subplot(1, 2, 2); hold on;

    legendEntries = {}; Ts = {}; X1s = {}; 
    fprintf('Controller costs to t=T for alpha=%i\n    d  &   costPPR    \n',alpha0)
    for d = 2:2:degree
        % uPPR = @(x) (-R\g{1}.'*kronPolyDerivEval(wPPR(1:d), x).'/2);
        uPPR = @(x) (kronPolyEval(Gains(1:d-1), x));

        % Simulate closed-loop dynamics
        [tPPR, XPPR] = ode45(@(t, x) F(x) + G(x) * uPPR(x), tspan, X0, options);

        % Plot PPR results
        subplot(1, 2, 1)
        plot(tPPR, XPPR(:, 1) / (pi / 180))
        legendEntries{end + 1} = sprintf('order-%i PPR controller ', d - 1);

        % Compute controller costs
        usPPR = zeros(length(tPPR),1); JsPPR = zeros(length(tPPR),1);
        for i = 1:length(tPPR)
            usPPR(i) = uPPR(XPPR(i, :).');
            JsPPR(i) = 0.25 * XPPR(i, :) * XPPR(i, :).' + uPPR(XPPR(i, :).') .^ 2;
        end
        cost_PPR = trapz(tPPR, JsPPR) / 2; % Approximate integral with trapezoidal rule
        fprintf('    %i  &  %f   \n', d, cost_PPR)

        subplot(1, 2, 2)
        plot(tPPR, usPPR / (pi / 180))

        % Add legends, plot range, etc.
        subplot(1, 2, 1)
        legend(legendEntries)
        ylim([0 alpha0 + 5])
        xlim(tspan)

        subplot(1, 2, 2)
        legend(legendEntries)
        ylim([-20 40])
        xlim(tspan)

        % Save for printing to text file 
        Ts{end + 1} = tPPR; X1s{end + 1} = XPPR(:, 1) / (pi / 180); 
    end

    if exportData
        fprintf('Writing data to plots/example7_alpha%i_x1.dat \n', alpha0)
        fileID = fopen(sprintf('plots/example7_alpha%i_x1.dat', alpha0), 'w');
        fprintf(fileID, '# Figure X-a Data\n');
        fprintf(fileID, '# aircraft stall angle, angle of attack data\n');

        % Calculate the maximum number of points in any line
        max_points = max(cellfun(@numel, Ts));

        % Determine the number of lines (sets of points)
        num_lines = length(Ts);

        % Write the header
        fprintf(fileID, '       t01     &      x01      & ');

        % Write the rest of the header
        for i = 2:num_lines - 1
            fprintf(fileID, '      t%02d     &      x%02d      & ', i, i);
        end
        fprintf(fileID, '      t%02d     &      x%02d      \n ', i + 1, i + 1);

        % Iterate over the number of points
        for j = 1:max_points
            % Iterate over each line
            for i = 1:num_lines
                t_line = Ts{i};
                x_line = X1s{i};
                if j <= numel(t_line) && abs(x_line(j)) < 90 % If this line has j or more points
                    fprintf(fileID, '%+1.6e & %+1.6e', t_line(j), x_line(j));
                    % Add '&' delimiter unless it's the last line
                    if i < num_lines
                        fprintf(fileID, ' & ');
                    else
                        fprintf(fileID, ' \n ');
                    end
                else % this line doesn't have a jth point
                    fprintf(fileID, '              &              ');
                    % Add '&' delimiter unless it's the last line
                    if i < num_lines
                        fprintf(fileID, ' & ');
                    else
                        fprintf(fileID, ' \n ');
                    end
                end
            end
        end
        % Close the data file
        fclose(fileID);
    end

end

end

function [value, isterminal, direction] = myEvent(T, Y)
value = (abs(Y(1)) > 90 * pi / 180);
isterminal = 1; % Stop the integration
direction = 0;
end
