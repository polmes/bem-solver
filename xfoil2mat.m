%% XFOIL DATA: alpha CL CD ...
filename = 'naca64A006_Re2M';
data = csvread([filename '.csv'], 11);

data_alpha = deg2rad([data(1, 1) - 75; data(1, 1) - 60; data(1, 1) - 20; data(:, 1); data(end, 1) + 20; data(end, 1) + 60; data(end, 1) + 75]);
data_Cl = [0; 0; 0; data(:, 2); 0; 0; 0];
data_Cd = [0.7; 0.7; 0.4; data(:, 3); 0.4; 0.7; 0.7];

%% PCHIP
pp_Cl = pchip(data_alpha, data_Cl);
pp_Cd = pchip(data_alpha, data_Cd);

%% EFFICIENCY
[mv, mi] = max(data(:, 2) ./ data(:, 3));
alpha_opt = deg2rad(data(mi, 1));

%% PLOTS

% Range
alphas = linspace(-30, 50, 100);

% Cl = f(alpha)
figure;
hold('on');
grid('on');

plot(data(:, 1), data(:, 2), 'x');
plot(alphas, ppval(pp_Cl, deg2rad(alphas)));

xlabel('$$\alpha$$ [$^{\circ}$]', 'Interpreter', 'latex');
ylabel('$$C_l$$', 'Interpreter', 'latex');
legend('Xfoil', 'PCHIP');

% Cd = f(alpha)
figure;
hold('on');
grid('on');

plot(data(:, 1), data(:, 3), 'x');
plot(alphas, ppval(pp_Cd, deg2rad(alphas)));

xlabel('$$\alpha$$ [$^{\circ}$]', 'Interpreter', 'latex');
ylabel('$$C_d$$', 'Interpreter', 'latex');
legend('Xfoil', 'PCHIP');

% E = Cl/Cd = f(alpha)
figure;
hold('on');
grid('on');

plot(data(:, 1), data(:, 2) ./ data(:, 3), 'x-');
plot(rad2deg(alpha_opt), data(mi, 2) ./ data(mi, 3), '*');

xlabel('$$\alpha$$ [$^{\circ}$]', 'Interpreter', 'latex');
ylabel('$$E = \frac{C_l}{C_d}$$', 'Interpreter', 'latex');

% Cd = f(Cl)
figure;
hold('on');
grid('on');

plot(data(:, 2), data(:, 3), 'x-');

xlabel('$$C_l$$', 'Interpreter', 'latex');
ylabel('$$C_d$$', 'Interpreter', 'latex');

%% SAVE
save([filename '.mat'], 'pp_Cl', 'pp_Cd', 'alpha_opt');
