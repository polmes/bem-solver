%% CONSTANTS
disp('INITIALIZING');
disp(' ');

% Global variables
global pp_Cl pp_Cd lambda_c rho Omega R n_blades theta_1 sigma_0 sigma_1 M_tip;

% Helicopter
M = 3300; % mass [kg]
DL = 350; % disk loading [N/m^2]
h = 1500; % 300; design altitude [m]
M_tip = 0.5; % tip Mach number
v_c = 0; % climbing speed [m/s]
r_root = 0.1; % non-dimensional root start position

% Other parameters
c_roots = [0.5 0.75];
n_blades_lim = 7;

% Numeric
r_val = 0.7;
N = 100;
DELTA = 1e-6;
EPSILON = 1e-6;
theta_0_k = deg2rad(15); % first initial guess for theta_0

% Debug
options = optimset('Display', 'off');

% Universal
g = 9.80665; % m/s^2

% Air
gamma = 1.4;
r_air = 287;

% ISA
T_0 = 288.15; % K
rho_0 = 1.225; % kg/m^3
T = T_0 - 6.5e-3 * h;
rho = rho_0 * (T/T_0)^(g/(6.5e-3 * r_air) - 1);

% Speed of sound
a = sqrt(gamma * r_air * T); % m/s

% Vertical flight at constant climbing speed (or hovering)
W = M * g;
T = W;

% Rotor
A = T/DL;
R = sqrt(A / (pi * (1 - r_root^2))); % radius
v_tip = a * M_tip;
Omega = v_tip / R; % rotation speed

% Airfoil
filename = 'naca64A006_Re2M';
load([filename '.mat']); % pp_Cl pp_Cd alpha_opt

%% IDEAL ROTOR
disp('IDEAL')
disp(' ');

tic;

% MTH
syms v_i;
mth = T == 2 * rho * A * v_i * (v_c + v_i);
assume(in(v_i, 'real') & v_i > 0);
sol_mth = solve(mth, v_i);
v_i = double(sol_mth);

% Non-dimensional velocities
lambda_i = v_i / (Omega * R);
lambda_c = v_c / (Omega * R);

% The variables
syms r sigma;
phi = atan((lambda_c + lambda_i) / r);
alpha = alpha_opt;
Cl_opt = ppval(pp_Cl, alpha_opt);
Cd_opt = ppval(pp_Cd, alpha_opt);

% BEM (ideal rotor): dT = dFz
bem = 8 * lambda_i * (lambda_c + lambda_i) * r == (r^2 + (lambda_c + lambda_i)^2) * (Cl_opt * cos(phi) - Cd_opt * sin(phi)) * sigma;
sol_bem = solve(bem, sigma);

% BEM (simplified ideal rotor)
bem_sim = 8 * lambda_i * (lambda_c + lambda_i) * r == (r^2) * (Cl_opt) * sigma;
sol_sim = solve(bem_sim, sigma);

% Main results: f(r)
theta = alpha + phi;
sigma = sol_bem;
sigma_sim = sol_sim;

toc;

disp(' ');

%% REAL ROTOR
disp('REAL');
disp(' ');

tic;

% The variables
syms r theta_0 lambda_i_real;

% Solidity
sigma_1 = double(subs(diff(sigma, r), r_val));
sigma_0 = double(subs(sigma, r_val)) - sigma_1 * r_val;
sigma_real = sigma_0 + sigma_1 * r;

% Chord
sigma_root = double(subs(sigma_real, r, r_root));
c_root = c_roots(1);
n_blades = ceil(sigma_root * pi * R / c_root);
if n_blades >= n_blades_lim
    c_root = c_roots(2);
    n_blades = ceil(sigma_root * pi * R / c_root);
end
c = sigma_real * pi * R / n_blades; % f(r)

% Mechanical angle
theta_1 = double(subs(diff(theta, r), r_val));
theta_real = theta_0 + theta_1 * r;

% More angles
phi_real = atan((lambda_c + lambda_i_real) / r);
alpha_real = theta_real - phi_real;

% BEM
rs = linspace(r_root, 1, N);

lambda_is = zeros(1, N);
% dFzs = sym(zeros(1, N));

k = 0;
delta = 1;
while delta > DELTA
    for e = 1:N
        r_e = rs(e);
        
        % f(lambda_i_real)
        sigma_e = double(subs(sigma_real, r, r_e));
        theta_e = double(subs(theta_real, [r theta_0], [r_e theta_0_k]));
        bem_loss = @(lambda_i_real) f_bem(lambda_i_real, r_e, sigma_e, theta_e);
                
        try
            lambda_is(e) = fsolve(bem_loss, lambda_is(e-1), options);
        catch
            lambda_is(e) = fsolve(bem_loss, lambda_i, options);
        end
    end
    
    % f(theta_0)
    int_bem = @(theta_0) integral(@(r) f_int_bem(r, theta_0, rs, lambda_is), r_root, 1);
    
    theta_0_prev = theta_0_k;
    theta_0_k = fsolve(@(theta_0) int_bem(theta_0) - T, theta_0_k, options);

    delta = abs(theta_0_k - theta_0_prev);
    % theta_0_k = ... % relaxation factor
    
    k = k + 1;
    disp(['k = ' num2str(k) ' -> ' num2str(delta)]);
end

% Redefine with found parameter
theta_real = subs(theta_real, theta_0, theta_0_k);

toc;

disp(' ');

%% PRANDTL + COMPRESSIBILITY LOSSES
disp('LOSSES');
disp(' ');

tic;

% The variables
syms r theta_0 lambda_i_loss;
theta_loss = theta_0 + theta_1 * r;


lambda_is_loss = zeros(1, N);
Fs = zeros(1, N);

k = 0;
delta = 1;
theta_0_k_loss = theta_0_k; % real value as initial guess
while delta > DELTA
    for e = 1:N
        r_e = rs(e);
        
        % Prandtl
        phi_trial = atan((lambda_c + lambda_is(e)) / r_e);
        epsilon = 1;
        while epsilon > EPSILON
            f_tip = (n_blades/2)*((1-r_e)/(r_e*sin(phi_trial)));
            f_root = (n_blades/2)*((r_e-r_root)/(r_e*sin(phi_trial)));
            F = (4/pi^2)*acos(exp(-f_tip))*acos(exp(-f_root));
            phi_prev = phi_trial;
            phi_trial = atan((lambda_c + (lambda_is(e)/F)) / r_e);
            
            epsilon = abs(phi_trial - phi_prev);
        end
        Fs(e) = F;
        
        % f(lambda_i_loss)
        sigma_e = double(subs(sigma_real, r, r_e));
        theta_e = double(subs(theta_loss, [r theta_0], [r_e theta_0_k_loss]));
        bem_loss = @(lambda_i_loss) f_bem_loss(lambda_i_loss, r_e, sigma_e, theta_e, F);
        
        try
            lambda_is_loss(e) = fsolve(bem_loss, lambda_is(e-1), options);
        catch
            lambda_is_loss(e) = fsolve(bem_loss, lambda_i, options);
        end
    end
    
    % f(theta_0)
    int_bem_loss = @(theta_0) integral(@(r) f_int_bem_loss(r, theta_0, rs, lambda_is_loss, Fs), r_root, 1);
    
    theta_0_prev = theta_0_k_loss;
    theta_0_k_loss = fsolve(@(theta_0) int_bem_loss(theta_0) - T, theta_0_k_loss, options);

    delta = abs(theta_0_k_loss - theta_0_prev);
    % theta_0_k = ... % relaxation factor
    
    k = k + 1;
    disp(['k = ' num2str(k) ' -> ' num2str(delta)]);
end

% Redefine with found parameter
theta_loss = subs(theta_loss, theta_0, theta_0_k_loss);

toc;

disp(' ');

%% FINAL CALCULATIONS
disp('CALCULATIONS');
disp('');

% Cl/d(r)
thetas = double(subs(theta_real, r, rs));
phis_real = atan((lambda_c + lambda_is) ./ rs);
alphas = thetas - phis_real;
Cls = ppval(pp_Cl, alphas);
Cds = ppval(pp_Cd, alphas);

thetas_loss = double(subs(theta_loss, r, rs));
phis_loss = atan((lambda_c + lambda_is_loss./Fs) ./ rs);
alphas_loss = thetas_loss - phis_loss;
Cls_loss = ppval(pp_Cl, alphas_loss);
Cds_loss = ppval(pp_Cd, alphas_loss);

% dFz/x(r)
phis = atan((lambda_c + lambda_i) ./ rs);
cs = double(subs(c, r, rs));
dFzs = 1/2 * rho * (Omega * R)^2 * (rs.^2 + (lambda_c + lambda_i)^2) * n_blades .* cs .* (Cl_opt * cos(phis) - Cd_opt * sin(phis)) * R;
dFxs = 1/2 * rho * (Omega * R)^2 * (rs.^2 + (lambda_c + lambda_i)^2) * n_blades .* cs .* (Cl_opt * sin(phis) + Cd_opt * cos(phis)) * R;
dFzs_real = 1/2 * rho * (Omega * R)^2 * (rs.^2 + (lambda_c + lambda_is).^2) * n_blades .* cs .* (Cls .* cos(phis_real) - Cds .* sin(phis_real)) * R;
dFxs_real = 1/2 * rho * (Omega * R)^2 * (rs.^2 + (lambda_c + lambda_is).^2) * n_blades .* cs .* (Cls .* sin(phis_real) + Cds .* cos(phis_real)) * R;
dFzs_loss = 1/2 * rho * (Omega * R)^2 * (rs.^2 + (lambda_c + lambda_is_loss).^2) * n_blades .* cs .* (Cls_loss .* cos(phis_loss) - Cds_loss .* sin(phis_loss)) * R;
dFxs_loss = 1/2 * rho * (Omega * R)^2 * (rs.^2 + (lambda_c + lambda_is_loss).^2) * n_blades .* cs .* (Cls_loss .* sin(phis_loss) + Cds_loss .* cos(phis_loss)) * R;

% Induced Power
Pi = T * (v_c + v_i); % MTH "ideal"
% dPis = 4 * pi * rho * (Omega * R)^3 * R^2 * lambda_is .* (lambda_c + lambda_is).^2 .* rs;
dPis = (1/2 * rho * (Omega * R)^2 * (rs.^2 + (lambda_c + lambda_is).^2) * n_blades .* cs .* (Cls .* sin(phis_real)) * R) .* (Omega * R .* rs);
Pi_real = trapz(rs, dPis);
% dPis_loss = 4*pi * (Omega * R)^3 * rho * R^2 * lambda_is_loss .* (lambda_c + lambda_is_loss).^2 .* rs;
dPis_loss = (1/2 * rho * (Omega * R)^2 * (rs.^2 + (lambda_c + lambda_is_loss).^2) * n_blades .* cs .* (Cls_loss .* sin(phis_loss)) * R) .* (Omega * R .* rs);
Pi_loss = trapz(rs, dPis_loss);

% Parasite Power
dPos = (1/2 * rho * (Omega * R)^2 * (rs.^2 + (lambda_c + lambda_is).^2) * n_blades .* cs .* (Cds .* cos(phis_real)) * R) .* (Omega * R .* rs);
Po = trapz(rs, dPos);
dPos_loss = (1/2 * rho * (Omega * R)^2 * (rs.^2 + (lambda_c + lambda_is_loss).^2) * n_blades .* cs .* (Cds_loss .* cos(phis_loss)) * R) .* (Omega * R .* rs);
Po_loss = trapz(rs, dPos_loss);

% Total Power
Pt = Pi_real + Po;
Pt_loss = Pi_loss + Po_loss;

disp(' ');

%% PLOTS
disp('POST-PROCESSING');
disp(' ');

% Solidity
figure;
hold('on');

fplot(sigma);
fplot(sigma_sim);
fplot(sigma_real);

xlabel('$$r = \frac{y}{R}$$', 'Interpreter', 'latex');
ylabel('$$\sigma$$', 'Interpreter', 'latex');
legend('Ideal', 'Simplified', 'Real');

grid('on');
xlim([0 1]);
ylim([0 inf]);

% Chord
figure;
hold('on');

fplot(c);
line([r_root r_root], [0 1], 'Color', 'k');

xlabel('$$r = \frac{y}{R}$$', 'Interpreter', 'latex');
ylabel('$$c$$ [m]', 'Interpreter', 'latex');

grid('on');
xlim([0 1]);
ylim([0 1]);

% Induced velocity
figure;
hold('on');

fplot(lambda_i);
plot(rs, lambda_is);
plot(rs, lambda_is_loss);

xlabel('$$r = \frac{y}{R}$$', 'Interpreter', 'latex');
ylabel('$$\lambda_i = \frac{v_i}{\Omega R}$$', 'Interpreter', 'latex');
legend('Ideal', 'Real', 'Prandtl + Compressibility', 'Location', 'southeast');

grid('on');
xlim([0 1]);
ylim([0 inf]);

% Mechanical angle
figure;
hold('on');

fplot(180/pi * theta);
fplot(180/pi * theta_real);
fplot(180/pi * theta_loss);

xlabel('$$r = \frac{y}{R}$$', 'Interpreter', 'latex');
ylabel('$$\theta$$ [$^{\circ}$]', 'Interpreter', 'latex');
legend('Ideal', 'Real', 'Prandtl + Compressibility');

grid('on');
xlim([0 1]);
ylim([0 inf]);

% Lift coefficient
figure;
hold('on');

plot(rs, ones(1, length(rs)) * Cl_opt);
plot(rs, Cls);
plot(rs, Cls_loss);

xlabel('$$r = \frac{y}{R}$$', 'Interpreter', 'latex');
ylabel('$$C_l$$', 'Interpreter', 'latex');
legend('Ideal', 'Real', 'Prandtl + Compressibility', 'Location', 'southeast');

grid('on');

% dFx/dr
figure;
hold('on');

plot(rs, dFxs_real);
plot(rs, dFxs_loss);

xlabel('$$r = \frac{y}{R}$$', 'Interpreter', 'latex');
ylabel('$$\frac{dF_x}{dr}$$ [N/m]', 'Interpreter', 'latex');
legend('Real', 'Prandtl + Compressibility', 'Location', 'southeast');

grid('on');

% F
figure;
hold('on');

plot(rs, Fs);

xlabel('$$r = \frac{y}{R}$$', 'Interpreter', 'latex');
ylabel('$$F$$', 'Interpreter', 'latex');

grid('on');

%% UTIL
function s = f_bem(lambda_i_real, r_e, sigma_e, theta_e) % f(lambda_i_real)
    global pp_Cl pp_Cd lambda_c;
    
    phi_e = atan((lambda_c + lambda_i_real) / r_e);
    alpha_e = theta_e - phi_e;
    Cl_e = ppval(pp_Cl, alpha_e);
    Cd_e = ppval(pp_Cd, alpha_e);
    
    s = 8 * lambda_i_real * (lambda_c + lambda_i_real) * r_e - (r_e^2 + (lambda_c + lambda_i_real)^2) * (Cl_e * cos(phi_e) - Cd_e * sin(phi_e)) * sigma_e;
end

function s = f_int_bem(r, theta_0, rs, lambda_is) % f(r)
    global pp_Cl pp_Cd lambda_c rho Omega R n_blades theta_1 sigma_0 sigma_1;
    
    lambda_i_e = interp1(rs, lambda_is, r);
    phi_e = atan((lambda_c + lambda_i_e) ./ r);
    theta_e = theta_0 + theta_1 * r;
    c_e = (sigma_0 + sigma_1 * r) * pi * R / n_blades;
    alpha_e = theta_e - phi_e;
    Cl_e = ppval(pp_Cl, alpha_e);
    Cd_e = ppval(pp_Cd, alpha_e);
    
    s = 1/2 * rho * (Omega * R)^2 * (r.^2 + (lambda_c + lambda_i_e).^2) * n_blades .* c_e .* (Cl_e .* cos(phi_e) - Cd_e .* sin(phi_e)) * R;
end

function s = f_bem_loss(lambda_i_loss, r_e, sigma_e, theta_e, F) % f(lambda_i_real)
    global pp_Cl pp_Cd lambda_c M_tip;
    
    phi_e = atan((lambda_c + lambda_i_loss / F) / r_e);
    alpha_e = theta_e - phi_e;
    Cl_e = ppval(pp_Cl, alpha_e);
    Cd_e = ppval(pp_Cd, alpha_e);
    
    s = 8 * lambda_i_loss * (lambda_c + lambda_i_loss) * r_e - (r_e^2 + (lambda_c + lambda_i_loss)^2) * (Cl_e / (sqrt(1-(M_tip^2)*(r_e^2 +...
        (lambda_c + lambda_i_loss)^2))) * cos(phi_e) - Cd_e * sin(phi_e)) * sigma_e;
end

function s = f_int_bem_loss(r, theta_0, rs, lambda_is, Fs) % f(r)
    global pp_Cl pp_Cd lambda_c rho Omega R n_blades theta_1 sigma_0 sigma_1;
    
    lambda_i_e = interp1(rs, lambda_is, r);
    F_e = interp1(rs, Fs, r);
    phi_e = atan((lambda_c + lambda_i_e./F_e) ./ r);
    theta_e = theta_0 + theta_1 * r;
    c_e = (sigma_0 + sigma_1 * r) * pi * R / n_blades;
    alpha_e = theta_e - phi_e;
    Cl_e = ppval(pp_Cl, alpha_e);
    Cd_e = ppval(pp_Cd, alpha_e);
    
    s = 1/2 * rho * (Omega * R)^2 * (r.^2 + (lambda_c + lambda_i_e).^2) * n_blades .* c_e .* (Cl_e .* cos(phi_e) - Cd_e .* sin(phi_e)) * R;
end
