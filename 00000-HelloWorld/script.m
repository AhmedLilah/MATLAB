%% Pre-Scirpt Routine
% clearing the command window
clc;

format shortEng;

%% Preparing the simulation paramters
% frequency scale
scale = 4e-9;

% time resolution
res = 5e-4;

% time data vecotr for analysis
t = -1 * scale : res*scale : 4*scale;

% impulse response data vector
stepMagintude = 1;
u = (t >= 0) * stepMagintude;

% system's natural frequency
wn = (2 * 2 * 3.14) / scale;

%damping factor
damping_ratio = 0.94429;

% Bounds for settling time
error_tolerance = 1 / (2^(12+1));
upper_bound = stepMagintude * (1 + error_tolerance);
lower_bound = stepMagintude * (1 - error_tolerance);

% Intermediate parameters for the equations
wd = wn*sqrt(1-damping_ratio^2);
theta = asin(sqrt(1-damping_ratio^2));
thetaDeg = theta * 180 / pi;
tanTheta = tan(theta);

% Responses vectors
y1 = [length(t)];
y2 = [length(t)];
y3 = [length(t)];


%% Graphing different system equations

% calculating time response for each system equation range
for indx = 1:length(t)
    y1(indx) = (1 - cos(wn*t(indx))) * u(indx);                                                                             % Undamped system equation
    y2(indx) = (1 - ((exp(-damping_ratio * wn*t(indx))/sqrt(1-damping_ratio^2)) * sin(wd*t(indx)+theta))) * u(indx);        % Under damped system
    y3(indx) = (1 - exp(-wn*t(indx))- wn*t(indx)*exp(-wn*t(indx))) * u(indx);                                               % Critcally damped system
end


figure;
hold on; % Hold the plot for overlaying multiple vectors
plot(t, u,  'linewidth', 1, 'DisplayName', 'Unit step');
% plot(t, y1, 'linewidth', 1, 'DisplayName', 'Undamped response');
% plot(t, y2, 'linewidth', 1, 'DisplayName', 'Under damped response');
plot(t, y3, 'linewidth', 1, 'DisplayName', 'Criticaly damped response');


% Customize plot
xlabel('Time');
ylabel('Magnitude');
title('Second Order System Response');
grid on;

% Add legend
legend('Location', 'best'); % Automatically position the legend
hold off;

%% Finding the settling time for the 
% Find the settling time
within_bounds = (y2 >= lower_bound) & (y2 <= upper_bound);
settling_index = find(~within_bounds, 1, 'last') + 1;

if settling_index <= length(t)
    settling_time = t(settling_index);
else
    % If the response never settles within bounds
    settling_time = t(lenght(t));
end

%% Calculating and plotting the sttling time for different damping factors
segma = 0.5 : res : 0.999;
settling_time_vs_damping_ratio = [length(segma)];

% figure;
% hold on; % Hold the plot for overlaying multiple vectors
upper_bound_vec = upper_bound * ones(1, length(t));
lower_bound_vec = lower_bound * ones(1, length(t));

% plot(t, upper_bound_vec, 'DisplayName', 'upper bound');
% plot(t, lower_bound_vec, 'DisplayName', 'lower bound');

y = [length(t)];

for i = 1 : length(segma)
    for indx = 1:length(t)
        wd = wn*sqrt(1-segma(i)^2);
        y(indx) = (1 - ((exp(-segma(i) * wn*t(indx))/sqrt(1-segma(i)^2)) * sin(wd*t(indx)+theta))) * u(indx);
    end

    % Plotting the Curve
    % plot(t, y,  'linewidth', 1, 'DisplayName', string(segma(i)));

    % Find the settling time
    within_bounds = (y >= lower_bound) & (y <= upper_bound);
    settling_index = find(~within_bounds, 1, 'last') + 1;

    if settling_index <= length(t)
        % If the response settles within bounds
        settling_time_vs_damping_ratio(i) = t(settling_index);
    else
        % If the response never settles within bounds
        settling_time_vs_damping_ratio(i) = t(length(t));
    end

end

% % Customize plot
% xlabel('Time');
% ylabel('System Response');
% title('Second Order System Response Vs Damping Factor');
% grid on;
% 
% % Add legend
% legend('Location', 'best'); % Automatically position the legend
% hold off;


%% Ploting the settling time vs damping factor

figure;
hold on; % Hold the plot for overlaying multiple vectors
plot(segma, settling_time_vs_damping_ratio,  'linewidth', 1, 'DisplayName', 'Settling Time');

% Customize plot
xlabel('Damping Ratio');
ylabel('Settling Time');
title('Second Order System Settling Time Vs Damping Factor');
grid on;

% Add legend
legend('Location', 'best'); % Automatically position the legend
hold off;

disp("Step Size: ");
disp(stepMagintude);
disp("Wn: ");
disp(wn);
disp("Damping Factor: ");
disp(damping_ratio);
disp("Error Tollerance: ");
disp(error_tolerance);
disp("Settling Time: ");
disp(settling_time);
disp("ts: ");
disp(wn);