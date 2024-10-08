clear all
close all
clc

%% Defining the mesh
n_points = 101;
dom_size = 1;
h = dom_size/(n_points - 1);
dt = 0.000001;
alpha = dt/(h*h);

%% Initialising the problem
y(n_points, n_points) = 0;
y(1,:) = 1;

y_new(n_points, n_points) = 0;
y_new(1,:) = 1;

error_mag = 1;
error_req = 1e-6;
iterations = 0;

%Tracking the error magnitude
error_track = 0;

%% Calculations
while error_mag > error_req
    for i = 2:(n_points-1)
        for j = 2:(n_points-1)
            y_new(i,j) = y(i,j) + alpha.*(y(i-1,j) + y(i+1,j) + y(i,j-1) + y(i,j+1) - 4*y(i,j));
        end
    end
    iterations = iterations + 1;
    % Calculation of error magnitude
    error_mag = 0;
    for i = 2:(n_points-1)
        for j = 2:(n_points-1)
            error_mag = error_mag + abs(y(i,j) - y_new(i,j));
            error_track(iterations) = error_mag;
        end
    end
    if rem(iterations, 1000) == 0
        iterations
        error_mag
    end
    % Assigning new to be old
    y = y_new;
end

%% Plotting
x_dom = ((1:n_points)-1).*h;
y_dom = 1-((1:n_points)-1).*h;
[X,Y] = meshgrid(x_dom,y_dom);
contourf(X,Y,y, 12)
colorbar

% Plot the error with time
figure;
time = dt.*(1:iterations);
plot(time, error_track)

%% Concept of subplots
figure;
subplot(2,1,1)
plot(time, error_track)
subplot(2,1,2)
contourf(X,Y,y, 12)
colorbar