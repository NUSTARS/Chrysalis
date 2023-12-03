% Flight Simulator

clc; clear; close all;

W = 5;
rho = 0.0023769;
g = 32.17;

CLift = 1.61; %%% CHECK
CDrag_0 = 0.16; %%% CHECK

chord = 4.5/12;
span = 22/12; %%% CHECK
S = span * chord;
AR = span/chord;
e = 0.7; %%% CHECK

m = W/g;

dt = 0.01;
t_0 = 0;
t_f = 1000;
N = (t_f-t_0)/dt;

time = zeros(N, 1);
time(1,:) = 0;
states = zeros(N, 5);

%%% STARTING CONDITIONS %%%
x_0 = 0;
y_0 = 400;
dx_0 = 0;
dy_0 = -20;

%%% SET UP %%%
x_axis_vec = [1,0];
velocity_vec = [dx_0, dy_0];
angle_0 = vectorAngle(x_axis_vec,velocity_vec);
states(1,:) = [x_0, y_0, dx_0, dy_0, angle_0];

vmag_vec = zeros(N,1);
impact = 500;

for ii = 2:length(time)

    time(ii) = dt*ii;
    
    vmag = sqrt(states(ii-1, 3)^2 + states(ii-1, 4)^2);

    vmag_vec(ii) = vmag;

    x_axis_vec = [1,0];
    velocity_vec = [states(ii-1, 3), states(ii-1, 4)] / vmag;

    glide_angle = vectorAngle(x_axis_vec,velocity_vec);

    L = 1/2 * rho * vmag^2 * CLift * S;
    CDrag = (CDrag_0) + (CLift^2 / (pi * AR * e));
    D = 1/2 * rho * vmag^2 * CDrag * S;

    % Transform aerodynamic forces from body frame to world frame
    L_world_xaxis = L * sin(glide_angle);
    L_world_yaxis = L * cos(glide_angle);
    D_world_xaxis = D * cos(glide_angle);
    D_world_yaxis = D * sin(glide_angle);

    world_y_forces = L_world_yaxis + D_world_yaxis - W;
    world_x_forces = L_world_xaxis - D_world_xaxis;

    %%% IN WORLD REFERENCE FRAME %%%
    ax = (1/m) * world_x_forces;
    ay = (1/m) * world_y_forces;

    vx = states(ii-1,3) + ax * dt;
    vy = states(ii-1, 4) + ay * dt;

    x = states(ii-1, 1) + vx * dt;
    y = states(ii-1, 2) + vy * dt;

    if isnan(vx) || isnan(vy) || isnan(x) || isnan(y)
        impact = ii - 1;
        disp("Error with signs — NAN");
    elseif vx > 10e3 || vy > 10e3
        impact = ii - 2;
        disp("Error with signs — EXPLODES");
    elseif y <= 0
        impact = ii - 1;
        break;
    else
        states(ii,:) = [x, y, vx, vy, glide_angle];
    end

end

impact_x = states(impact,1);
starting_y = states(1,2);
glide_ratio = impact_x / starting_y;
disp("Glide ratio: " + glide_ratio);

impact_velocity = vmag_vec(impact);
disp("Impact velocity: " + impact_velocity + " ft/s");

impact_y_velocity = states(impact,4);
disp("Y-direction: " + impact_y_velocity + " ft/s");

impact_x_velocity = states(impact,3);
disp("X-direction: " + impact_x_velocity + " ft/s");


% Create a single figure window
figure('Position', [0, 0, 1000, 1500]);

% Plotting the first subplot
subplot(2, 2, 1);
plot(states(1:impact, 1), states(1:impact, 2))
title("Vertical Height vs. Horizontal Distance")
ylabel("Height (ft)")
xlabel("Distance (ft)")

% Plotting the second subplot
subplot(2, 2, 2);
plot(time(1:impact),states(1:impact,2))
title("Vertical Height vs. Time")
ylabel("Height (ft)")
xlabel("Time (s)")

% Plotting the third subplot
subplot(2, 2, 3);
plot(time(1:impact),states(1:impact,2))
ylabel("Height (ft)")
hold on;
yyaxis right
ylabel("Velocity (ft/s)")
plot(time(1:impact),states(1:impact,3));
plot(time(1:impact),states(1:impact,4));
plot(time(1:impact),vmag_vec(1:impact));
legend("Height", "x-dot", "y-dot", "v-mag")
title("Vertical Height vs. Time with Velocities")
hold off

% Plotting the fourth subplot
subplot(2, 2, 4);
plot(time(1:impact),states(1:impact,3));
ylabel("Velocity (ft/s)")
hold on;
plot(time(1:impact),states(1:impact,4));
yyaxis right
ylabel("Glide Angle (º)")
plot(time(1:impact),states(1:impact, 5)*180/pi)
legend("x-dot", "y-dot", "Glide Angle")
title("Velocities and Glide Angle over Time")
hold off;

% hold on;
% figure;
% yyaxis left
% plot(time,states(:,3))
% plot(time,states(:,4))
% plot(time,vmag_vec);
% yyaxis right
% plot(time,psi_vec)
% title("Vmag and Psi")
% hold off;

function angle = vectorAngle(a, b)
    % Check if the vectors have the same size
    if numel(a) ~= numel(b)
        error('Vectors must have the same size');
    end
    
    % Calculate the dot product of the vectors
    dotProduct = dot(a, b);
    
    % Calculate the magnitudes of the vectors
    magnitudeA = norm(a);
    magnitudeB = norm(b);
    
    % Check if the magnitudes are non-zero
    if magnitudeA == 0 || magnitudeB == 0
        %error('One or both vectors have zero magnitude, angle is undefined');
        angle = 0;
        return;
    end
    
    % Calculate the cosine of the angle
    cosineTheta = dotProduct / (magnitudeA * magnitudeB);
    
    % Ensure the cosine value is within the valid range [-1, 1]
    cosineTheta = max(min(cosineTheta, 1), -1);
    
    % Calculate the angle in radians
    angle = acos(cosineTheta);
end
