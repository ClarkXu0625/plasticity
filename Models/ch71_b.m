%% parameter
theta_E = -5;
theta_I = 0;
r_max = 100; % Hz
alpha_E = 0.05;
alpha_I = 1;
W_EE = 2;
W_EI = 2.5;
W_IE = -2.5;
W_II = -2;

W_EI_X = 1.75;

%% Time vector
dt = 0.1e-3; % 0.1 ms
t = 0:dt:3;

%% Create vector for extra current
Iapp_E_base = 25 * ones(1, length(t));   % Iapp excitatory base
Iapp_I_base = 20 * ones(1, length(t));   % Iapp inhibitory base

Istim_X = zeros(1, length(t));
Istim_X(1/dt:1.1/dt) = 10;
Istim_Y = zeros(1, length(t));
Istim_Y(2/dt:2.1/dt) = 10;

Iapp_E_X = Iapp_E_base;
Iapp_I_X = Iapp_I_base + Istim_X;
Iapp_E_Y = Iapp_E_base;
Iapp_I_Y = Iapp_I_base + Istim_Y;

tau_E = 5e-3;
tau_I = 5e-3;

%% firing rate r vector
r_E_X = zeros(1, length(t));
r_I_X = zeros(1, length(t));
I_E_X = zeros(1, length(t));
I_I_X = zeros(1, length(t));

r_E_Y = zeros(1, length(t));
r_I_Y = zeros(1, length(t));
I_E_Y = zeros(1, length(t));
I_I_Y = zeros(1, length(t));


for i = 2: length(t)
    %% I for cell X
    I_E_X(1, i-1) = W_EE*r_E_X(1, i-1) + W_IE*r_I_X(1, i-1) + Iapp_E_X(1, i-1);
    I_I_X(1, i-1) = W_EI_X * r_E_Y(1, i-1) + W_EI*r_E_X(1, i-1) + W_II*r_I_X(1, i-1) + Iapp_I_X(1, i-1);


    dr_E_X = (-r_E_X(1, i-1) + alpha_E * ((I_E_X(1, i-1)-theta_E)^2) * sign(I_E_X(1, i-1)-theta_E))*dt/tau_E;
    r_E_X(1, i) = r_E_X(1, i-1)+dr_E_X;

    % adjust the firing rate if it is out of r_max range
    if r_E_X(1, i) < 0
        r_E_X(1, i) = 0;
    elseif r_E_X(1, i) > r_max
        r_E_X(1, i) = r_max;
    end

    dr_I_X = (-r_I_X(1, i-1) + alpha_I * (I_I_X(1, i-1) - theta_I))* dt/tau_I;
    r_I_X(1, i) = r_I_X(1, i-1) +dr_I_X;

    % adjust the firing rate if it is out of r_max range
    if r_I_X(1, i) < 0
        r_I_X(1, i) = 0;
    elseif r_I_X(1, i) > r_max
        r_I_X(1, i) = r_max;
    end


    %% Cell Y
    I_E_Y(1, i-1) = W_EE*r_E_Y(1, i-1) + W_IE*r_I_Y(1, i-1) + Iapp_E_Y(1, i-1);
    I_I_Y(1, i-1) = W_EI*r_E_Y(1, i-1) + W_II*r_I_Y(1, i-1) + Iapp_I_Y(1, i-1);

    dr_E_Y = (-r_E_Y(1, i-1) + alpha_E * ((I_E_Y(1, i-1)-theta_E)^2) * sign(I_E_Y(1, i-1)-theta_E))*dt/tau_E;
    r_E_Y(1, i) = r_E_Y(1, i-1)+dr_E_Y;

    % adjust the firing rate if it is out of r_max range
    if r_E_Y(1, i) < 0
        r_E_Y(1, i) = 0;
    elseif r_E_Y(1, i) > r_max
        r_E_Y(1, i) = r_max;
    end

    dr_I_Y = (-r_I_Y(1, i-1) + alpha_I * (I_I_Y(1, i-1) - theta_I))* dt/tau_I;
    r_I_Y(1, i) = r_I_Y(1, i-1) +dr_I_Y;

    % adjust the firing rate if it is out of r_max range
    if r_I_Y(1, i) < 0
        r_I_Y(1, i) = 0;
    elseif r_I_Y(1, i) > r_max
        r_I_Y(1, i) = r_max;
    end
end



set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

 
    

plot(t, r_I_X(1, :));
legend_vec{1} = "r_I_X";

hold on
plot(t, r_E_X(1, :));
legend_vec{2} = "r_E_X";

hold on
plot(t, r_I_Y(1, :));
legend_vec{3} = "r_I_Y";

hold on
plot(t, r_E_Y(1, :));
legend_vec{4} = "r_E_Y";

xlabel("time (s)")
ylabel("firing rate (Hz)")
legend(legend_vec);

ylim([0 25])


