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

%% Time vector
dt = 1e-4; % 0.1 ms
t = 0:dt:3;

%% Create vector for extra current
Iapp_E = zeros(4, length(t));   % Iapp excitatory base
Iapp_I = zeros(4, length(t));   % Iapp inhibitory base
Istim_I = zeros(1, length(t));
Istim_I(1/dt:2/dt) = 20;

Iapp_E(1, :) = 5*ones(1, length(t));
Iapp_I(1, :) = 5*ones(1, length(t)) + Istim_I;
Iapp_E(2, :) = 5*ones(1, length(t));
Iapp_I(2, :) = 10*ones(1, length(t)) + Istim_I;
Iapp_E(3, :) = 5*ones(1, length(t));
Iapp_I(3, :) = 15*ones(1, length(t)) + Istim_I;
Iapp_E(4, :) = 5*ones(1, length(t));
Iapp_I(4, :) = 40*ones(1, length(t)) + Istim_I;


%% firing rate r vector
r_E = zeros(4, length(t));
r_I = zeros(4, length(t));
I_E = zeros(4, length(t));
I_I = zeros(4, length(t));

%% initialize parameters
tau_E = 5e-3 * ones(1, 4);
tau_I = 5e-3 * ones(1, 4);

tau_E(1, 3) = 2e-3;
tau_E(1, 4) = 2e-3;
tau_I(1, 3) = 10e-3;
tau_I(1, 4) = 10e-3;


for question_number = 1: 4

    for i = 2: length(t)
        %% get I value for last time slice
        I_E(question_number, i-1) = W_EE*r_E(question_number, i-1) + ...
            W_IE*r_I(question_number, i-1) + Iapp_E(question_number, i-1);
        
        I_I(question_number, i-1) = W_EI*r_E(question_number, i-1) + ...
            W_II*r_I(question_number, i-1) + Iapp_I(question_number, i-1);

        dr_E = (-r_E(question_number, i-1) + alpha_E * ((I_E(question_number, i-1)-theta_E)^2) * ...
            sign(I_E(question_number, i-1)-theta_E))*dt/tau_E(1, question_number);


        r_E(question_number, i) = r_E(question_number, i-1) + dr_E;

        % adjust the firing rate if it is out of r_max range
        if r_E(question_number, i) < 0
            r_E(question_number, i) = 0;
        elseif r_E(question_number, i) > r_max
            r_E(question_number, i) = r_max;
        end


        dr_I = (-r_I(question_number, i-1) + alpha_I * (I_I(question_number, i-1) - theta_I))* ...
            dt/tau_I(1, question_number);

        r_I(question_number, i) = r_I(question_number, i-1) + dr_I;

        % adjust the firing rate if it is out of r_max range
        if r_I(question_number, i) < 0
            r_I(question_number, i) = 0;
        elseif r_I(question_number, i) > r_max
            r_I(question_number, i) = r_max;
        end
    end
end


set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');



for i = 1: 4
    figure(i);
    
    

    plot(t, r_I(i, :));
    legend_vec{1} = "r_I";
    
    hold on
    plot(t, r_E(i, :));
    legend_vec{2} = "r_E";

    xlabel("time (s)")
    ylabel("firing rate (Hz)")
    legend(legend_vec);
    
    title("Question "+ int2str(i))
end

