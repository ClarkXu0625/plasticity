%% Parameters
tau = 0.020;
r_max = 100;
I_threshold = 50;
I_delta = 5;
W_initial = 0.2;
sigma_I = 1; % noise term
threshold = 40; % Hz
trials = 800;

rule = 5;

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
prediction = zeros(1, trials); % create an array of predictions

mean_reward = 0.5*ones(5,2); % used for rule B
mean_reward_temp = zeros(5,2);

synaptic_strength = zeros(trials, 5, 2);    % record the synaptic strength

disp(size(synaptic_strength))

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
                    
    end
    %% Simulation ends here
        
    active_unit=prediction(1,trial);
    
    % Alter the connections by rules
    switch rule
        
        case 1 
            %% rule A
            E=0.5;
            plasticity_rate = 0.04;
            dW = E*plasticity_rate;
            
            % without a decision
            if active_unit==0
                for i=1:5
                    if active_input(i)
                        W_input(i,:) = W_input(i,:)+dW;
                    end
                end
            else
               
                for i = 1:5

                    % alter connection from active input unit to active
                    % decision unit only
                    if active_input(i)

                        % correct prediction
                        if active_unit==weather
                            E = 0.5;
                            % to the positive active unit on coorect trials
                            if W_input(i,active_unit)>0
                                W_input(i,active_unit)=W_input(i,active_unit)+E*plasticity_rate;
                            end
                            
                            % to the negative inactive unit on correct trials
                            if W_input(i,3-active_unit)<0
                                W_input(i,3-active_unit)=W_input(i,3-active_unit)-E*plasticity_rate;
                            end
                          
                        else
                            E=-0.5;
                            % to the negative active unit on incorrect trial
                            if W_input(i,active_unit)<0
                                W_input(i,active_unit)=W_input(i,active_unit)+E*plasticity_rate;
                            end 
                               
                            % to the positive inactive unit on incorrect trial
                            if W_input(i,3-active_unit)>0
                                W_input(i,3-active_unit)=W_input(i,3-active_unit)-E*plasticity_rate;
                            end
                        end

                    end
                end
            end
        
        case 2 
            %% rule B
            plasticity_rate = 0.04;
            
            % matrix of Reward in this trial
            if weather==active_unit
                reward(:,weather)=active_input;
            end
            
            % add the reward to temp to calculate mean reward
            mean_reward_temp = mean_reward_temp+reward;
            
            % update mean reward every 10 trials
            if mod(trial,10)==0
                mean_reward = mean_reward_temp/10;
                mean_rewared_temp=0;
            end

            E=reward-mean_reward;

            % without a decision
            if active_unit==0
                for i=1:5
                    if active_input(i)
                        W_input(i,:) = W_input(i,:)+plasticity_rate*E(i,:);
                    end
                end
            else
               
                for i = 1:5

                    % alter connection from active input unit to active
                    % decision unit only
                    if active_input(i)

                        % correct prediction
                        if active_unit==weather
                            % to the positive active unit on coorect trials
                            if W_input(i,active_unit)>0
                                W_input(i,active_unit)=W_input(i,active_unit)+...
                                    E(i,active_unit)*plasticity_rate;
                            end
                            
                            % to the negative inactive unit on correct trials
                            if W_input(i,3-active_unit)<0
                                W_input(i,3-active_unit)=W_input(i,3-active_unit)-...
                                    E(i,3-active_unit)*plasticity_rate;
                            end
                          
                        else
                            % to the negative active unit on incorrect trial
                            if W_input(i,active_unit)<0
                                W_input(i,active_unit)=W_input(i,active_unit)+...
                                    E(i,active_unit)*plasticity_rate;
                            end 
                               
                            % to the positive inactive unit on incorrect trial
                            if W_input(i,3-active_unit)>0
                                W_input(i,3-active_unit)=W_input(i,3-active_unit)-...
                                    E(i,3-active_unit)*plasticity_rate;
                            end
                        end

                    end
                end
            end

        
        case 3
            %% Rule C
            E=0.5;
            W_max=2*W_initial;
            plasticity_rate = 0.04;
            dW = E*plasticity_rate;
            
            % without a decision
            if active_unit==0
                for i=1:5
                    if active_input(i)
                        W_input(i,:) = W_input(i,:)+dW;
                    end
                end
            else
               
                for i = 1:5

                    % alter connection from active input unit to active
                    % decision unit only
                    if active_input(i)

                        % correct prediction
                        if active_unit==weather
                            E = 0.5;
                            % to the positive active unit on coorect trials
                            if W_input(i,active_unit)>0
                                W_input(i,active_unit)=W_input(i,active_unit)+...
                                    E*plasticity_rate*(W_max-W_input(i,active_unit));
                            end
                            
                            % to the negative inactive unit on correct trials
                            if W_input(i,3-active_unit)<0
                                W_input(i,3-active_unit)=W_input(i,3-active_unit)-...
                                    E*plasticity_rate*(W_max-W_input(i,3-active_unit));
                            end
                          
                        else
                            E=-0.5;
                            % to the negative active unit on incorrect trial
                            if W_input(i,active_unit)<0
                                W_input(i,active_unit)=W_input(i,active_unit)+...
                                    E*plasticity_rate*W_input(i,active_unit);
                            end 
                               
                            % to the positive inactive unit on incorrect trial
                            if W_input(i,3-active_unit)>0
                                W_input(i,3-active_unit)=W_input(i,3-active_unit)-...
                                    E*plasticity_rate*W_input(i,3-active_unit);
                            end
                        end

                    end
                end
            end
        
        
        case 4
            %% Case D
            
            plasticity_rate = 0.04;
            
            % matrix of Reward in this trial
            if weather==active_unit
                reward(:,weather)=active_input;
            end
            
            % add the reward to temp to calculate mean reward
            mean_reward_temp = mean_reward_temp+reward;
            
            % update mean reward every 10 trials
            if mod(trial,10)==0
                mean_reward = mean_reward_temp/10;
                mean_rewared_temp=0;
            end

            E=reward-mean_reward;

            % without a decision
            if active_unit==0
                for i=1:5
                    if active_input(i)
                        W_input(i,:) = W_input(i,:)+plasticity_rate*E(i,:);
                    end
                end
            else
               
                for i = 1:5

                    % alter connection from active input unit to active
                    % decision unit only
                    if active_input(i)

                        % correct prediction
                        if active_unit==weather
                            % to the positive active unit on coorect trials
                            if W_input(i,active_unit)>0
                                W_input(i,active_unit)=W_input(i,active_unit)+...
                                    E(i,active_unit)*plasticity_rate;
                            end
                            
                            
                        else
                            % to the negative active unit on incorrect trial
                            if W_input(i,active_unit)<0
                                W_input(i,active_unit)=W_input(i,active_unit)+...
                                    E(i,active_unit)*plasticity_rate;
                            end                                
                            
                        end

                    end
                end
            end

        case 5
            %% Rule E
            plasticity_rate = 0.04;
            
            % matrix of Reward in this trial
            if weather==active_unit
                reward(:,weather)=active_input;
            end
            
            % add the reward to temp to calculate mean reward
            mean_reward_temp = mean_reward_temp+reward;
            
            % update mean reward every 10 trials
            if mod(trial,10)==0
                mean_reward = mean_reward_temp/10;
                mean_rewared_temp=0;
            end

            E=reward-mean_reward;

            % without a decision
            if active_unit==0
                for i=1:5
                    if active_input(i)
                        W_input(i,:) = W_input(i,:)+plasticity_rate*E(i,:);
                    end
                end
            else
               
                for i = 1:5

                    % alter connection from active input unit to active
                    % decision unit
                    if active_input(i)

                        % correct prediction
                        if active_unit==weather
                            % to the positive active unit on coorect trials
                            if W_input(i,active_unit)>0
                                W_input(i,active_unit)=W_input(i,active_unit)+...
                                    E(i,active_unit)*plasticity_rate;
                            end
                            
                            
                        else
                            % to the negative active unit on incorrect trial
                            if W_input(i,active_unit)<0
                                W_input(i,active_unit)=W_input(i,active_unit)+...
                                    E(i,active_unit)*plasticity_rate;
                            end                                
                            
                        end

                    else
                        % Update the connection from inactive input unit to
                        % active decision-making unit
                        W_input(i,active_unit)=W_input(i,active_unit)-...
                                    E(i,active_unit)*plasticity_rate;
                    end
                end
            end




    end

    if trial>700        
        if prediction==weather
            correct_response(1) = correct_response(1)+1;
        end
        
    end
    
    
    %W_input=min(W_input, 2*W_initial);
    %W_input=max(W_input,0);


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
figure(99)
scatter(x,y)
xlabel('log probability')
ylabel('synaptic strength difference')

