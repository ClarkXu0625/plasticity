%% Parameters
tau = 0.020;
r_max = 100;
I_threshold = 50;
I_delta = 5;
W_initial = 0.2;
W_max=2*W_initial;
sigma_I = 1; % noise term
threshold = 40; % Hz
trials = 800;
plasticity_rate=0.04;
plasticity_rule = 2;

%% set up probability. 
% (i, j) refers the probability of P(Unit i | condition j)
prob_matrix = zeros(5, 2);
prob_matrix(1, 2) = 0.95;
prob_matrix(2, 2) = 0.75;
prob_matrix(3, 2) = 0.5;
prob_matrix(4, 2) = 0.25;
prob_matrix(5, 2) = 0.05;

for i = 1:5
    prob_matrix(i, 1) = 1-prob_matrix(i, 2);
end

%% set up time
tmax = 0.500;
dt = 0.001;
t = 0:dt:tmax;
Nt = length(t);


% set up W_input and W_recurrent
W_input = W_initial*ones(5, 2);
W_recurrent = 0.5*ones(2, 2);
W_recurrent(1, 2) = -0.5;
W_recurrent(2, 1) = -0.5;

% noise current
I_noise = randn(Nt, 2)*(sigma_I/sqrt(dt));

% count number of correct response
correct_response = zeros(1, 2);

reward=zeros(5,2);


synaptic_strength = zeros(trials, 5, 2);    % record the synaptic strength
correct = zeros(1,trials); % if the system make correct predictions in each trial
prediction = ones(1, trials); % create an array of predictions
cumulative_correctness=zeros(1,trials);



for trial=1:trials
    % for every trial, record the last trial's synaptic strength
    synaptic_strength(trial,:,:) = W_input;


    % determine the weather randomly
    % 1 is rainy and 2 is sunny
    weather = 2-(rand()>0.5);
    
    current = zeros(1, 2);
    rate = zeros(Nt, 5);
    rate_circuit = zeros(Nt, 2);
    
    
    % determine the fired input units in this trial by the given
    % probability matrix
    firing = find(rand(5, 2) < prob_matrix );
    fired_inputs = zeros(5, 2);
    fired_inputs(firing) = 1;    

    % set up input current corresponding to cues
    input_current = zeros(Nt, 5);
    input_current(0.1/dt+1:end,:) = 50;
    active_input=fired_inputs(:,weather);
    
    % alter the input current
    for i = 1:5
        if fired_inputs(i, weather)==0
            input_current(:,i) = 0;
        end
    end
        
    %% do the 500 ms simulation
    for i = 2:Nt
        % Calculate the input firing rate
        ssr_input = r_max./(1+exp((I_threshold*ones(1,5)-input_current(i, :))/I_delta));        
        rate(i,:) = rate(i-1,:)+(ssr_input-rate(i-1,:))*(dt/tau);

        % Calculate intput current for decision-making units
        current(1) = rate(i-1, :)*W_input(:, 1)+ ...
                       rate_circuit(i-1,:)*W_recurrent(:,1)+...
                       I_noise(i-1, 1);
        current(2) = rate(i-1, :)*W_input(:, 2)+ ...
                       rate_circuit(i-1,:)*W_recurrent(:,2)+...
                       I_noise(i-1, 2);
        
        % Calculate the firing rate for the two units in the decision 
        % making circuit
        ssr = r_max./(1+exp((-current+I_threshold)/I_delta));
        rate_circuit(i,:) = rate_circuit(i-1,:) + (ssr-rate_circuit(i-1,:))*(dt/tau);
        

        
        % determine which unit first reach the threshold
        % do not alter the prediction after it's predicted
        if prediction(1, trial)==0            
            if rate_circuit(i, 1)>=threshold
                prediction(1, trial) = 1;
            elseif rate_circuit(i, 2)>=threshold
                prediction(1, trial) = 2;
            end
        end
        %%%%%%
        %prediction(trial) = 1;                             % find which decision is made, start by assuming 1
        %maxrate = rate_circuit(Nt,prediction(trial));       % initialize maximum rate at end of simulation
        
        %if rate_circuit(Nt,2) > maxrate   % if rate of other cell is higher than max
            %prediction(trial) = 2;                  % that cell is the winner
            %maxrate = rate_circuit(length(t),2); % reset max rate
        %end
        
                    
    end
    %% Simulation ends here
        
    active_unit=prediction(1,trial);
    correct(1,trial)=(active_unit==weather);


    switch plasticity_rule
        case 1
            rew_pred_err = (correct(1,trial) - 0.5);
        otherwise
            if ( trial > 10 )
                rew_pred_err = correct(1,trial) - mean(correct(1,trial-10:trial-1));
            else
                rew_pred_err = correct(1,trial) - 0.5;
            end
    end
    
    switch plasticity_rule
        case { 1,2 }
            for cell = find(active_input)         % potentiate synapses from stimuli that were on
                if active_unit==0
                    W_input(cell,weather) = W_input(cell,weather) - ...
                        rew_pred_err*plasticity_rate;
                else
                    W_input(cell,active_unit) = W_input(cell,active_unit) + ...
                        rew_pred_err*plasticity_rate;
                    W_input(cell,3-active_unit) = W_input(cell,3-active_unit) - ...
                        rew_pred_err*plasticity_rate;
                end
             end
                 
        case 3
            if rew_pred_err > 0                 % with a correct response
                for cell = find(active_input)         % potentiate synapses from stimuli that were on
                    W_input(cell,active_unit) = W_input(cell,active_unit) + ...
                        rew_pred_err*plasticity_rate*(W_max-W_input(cell,active_unit))/W_initial;
                    W_input(cell,3-active_unit) = W_input(cell,3-active_unit) - ...
                        rew_pred_err*plasticity_rate*W_input(cell,3-active_unit)/W_initial;
                end
            else
                for cell = find(active_input)         % potentiate synapses from stimuli that were on
                    W_input(cell,active_unit) = W_input(cell,active_unit) + ...
                        rew_pred_err*plasticity_rate*W_input(cell,active_unit)/W_initial;
                    W_input(cell,3-active_unit) = W_input(cell,3-active_unit) - ...
                        rew_pred_err*plasticity_rate*(W_max-W_input(cell,3-active_unit))/W_initial;
                end
                
            end
        case 4
                for cell = find(active_input)         % potentiate synapses from stimuli that were on
                    W_input(cell,active_unit) = W_input(cell,active_unit) + ...
                        rew_pred_err*plasticity_rate;
                end
        case 5
            for cell = find(active_input)         % potentiate synapses from stimuli that were on
                W_input(cell,active_unit) = W_input(cell,active_unit) + ...
                    rew_pred_err*plasticity_rate;
            end
            for cell = find(~active_input)         % depress synapses from stimuli that were not on
                W_input(cell,active_unit) = W_input(cell,active_unit) - ...
                    rew_pred_err*plasticity_rate;
            end
    end

    decay_trials = 50;
    smooth_correct = zeros(1,trials);
    for i = 1:trials
        smooth_correct(i) = mean(exp(([1:i]-i)/decay_trials).*correct(1:i))/mean(exp(([1:i]-i)/decay_trials));
    end
    if trial<=100
        cumulative_correctness(1,trial)=sum(correct(1:trial))/trial;
    else
        cumulative_correctness(1,trial)=sum(correct(trial-100:trial))/100;
    end
    %W_input = max(W_input,0);               % weights can not be negative
    %W_input = min(W_input,Wmax);            % or greater than some maximum
end

for i = 1:5
    for j = 1:2
        figure(j)
        plot(synaptic_strength(:,i,j))
        
        hold on
        
        ylabel("connection strength")
        xlabel("trial number")
    end
end
figure(1)
title("rainy")
legend('unitA','unitB','unitC','unitD','unitE',...
            'Location','northwest')
figure(2)
title("sunny")
legend('unitA','unitB','unitC','unitD','unitE',...
            'Location','northwest')


x = zeros(1, 5);
y = zeros(1, 5);
for i =1:5
    y(i) = W_input(i,2)-W_input(i,1);
    x(i) = log(prob_matrix(i,2)/prob_matrix(i,1));
end

figure(3)
plot(1:trials,smooth_correct)

figure(99)
plot(x,y)
xlabel('log probability')
ylabel('synaptic strength difference')

