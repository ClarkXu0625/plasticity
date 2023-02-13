N1 = 17;
viewmatrix = zeros(N1,N1); % square array for viewing
Ncells = numel(viewmatrix); % total number of cells
Ntrials = 400;  % Number of trials
inputstrength = 50; % input strength during training


%% Time simulation setup
dt = 0.001;         % time step for simulation
tau_r = 0.010;        % time constant for cells
tmax = 1;           % maximum time to wait
t = 0:dt:tmax;      %time vector
Nt = length(t);     % number of time points

%% Parameters
rmax = 50;
gain = 1;
threshold = 10; 
del_I = 1;


%% Initialize patterns
Npatterns = 4; % Number of stimulus patterns to define and learn
pattern1 = zeros(size(viewmatrix));
pattern2 = zeros(size(viewmatrix));
pattern3 = zeros(size(viewmatrix));
pattern4 = zeros(size(viewmatrix));


%% create pattern
% First pattern is letter "X"
for i = 1:N1
    pattern1(i,i) = 1;                  % Top left to bottom right
    pattern1(i,N1+1-i) = 1;             % Top right to bottom left
end

%Second pattern is letter "Y"
for i = 1:(N1)/2-1
    pattern2(i,i+2) = 1;                % Top left to center
    pattern2(i,i+3) = 1;                % Double the thickness
    pattern2(i,N1-1-i) = 1;             % Top right to center
    pattern2(i,N1-2-i) = 1;             % Double the thickness
end
for i = ceil((N1+2)/2)-3:N1
    pattern2(i,floor((N1+2)/2)) = 1;    % vertical line in center
end

% 3rd pattern is letter "Z"
for j = 1:N1
    pattern3(1,j) = 1;                  % Top horizontal line
    pattern3(N1,j) = 1;                 % Bottom horizontal line
    if ( j > 1 )
        pattern3(N1+2-j,j) = 1;          % Diagonal of the "Z"
    end
end

% 4th pattern is letter "O"
for i = 2:N1-1
    pattern4(i,2) = 1;                  % Top horizontal line
    pattern4(i,N1-1) = 1;               % Bottom horizontal line
end
for j = 2:N1-1
    pattern4(2,j) = 1;                  % Left vertical line
    pattern4(N1-1,j) = 1;               % Right vertical line
end

%% Create pattern collection
pattern_collection = zeros(N1, N1, Npatterns);
pattern_collection(:,:,1) = pattern1;
pattern_collection(:,:,2) = pattern2;
pattern_collection(:,:,3) = pattern3;
pattern_collection(:,:,4) = pattern4;

%% set up W
Wmean = -0.3/Ncells;            
W = zeros(Ncells) + Wmean; 
Wmax = 8/Ncells;    
Wmin = -8/Ncells;    
flip_probability = 0.1;


%% define the rules to alter connections
epsilon = 0.1/Ncells;
r_T = 25;


%% Set up loop
for trials = 1: Ntrials
    
    % randomly select one pattern
    pattern = randi(Npatterns);
    input_pattern = pattern_collection(:,:,pattern);
    
    % array of random number from 0 to 1
    flip = rand(Ncells, 1);
    flipinputs = find(flip < flip_probability);
    % flip the number with 0.1 probability
    input_pattern(flipinputs) = 1-input_pattern(flipinputs);  
    
    % firing rate initialized to zero
    rate = zeros(Nt, Ncells);
    current = zeros(Ncells, 1);


    for i = 2: Nt

        % Applied current only exist in the first half of simulation
        applied_current = zeros(N1);
        if i < Nt/2
            applied_current = input_pattern*50;          
        end

        current = rate(i-1,:)*W + reshape(applied_current, [1, Ncells]);
       
        % Calculate rate
        rss = rmax./(1+exp(-del_I*(current-threshold)));         
        rate(i,:) = rate(i-1,:) + (rss-rate(i-1,:))*dt/tau_r;

    end        
    
    dW = epsilon*(rate'>r_T)*((rate>r_T));
    dW = dW - (epsilon/4)* ((rate'<r_T))*((rate>r_T));

    W = W+dW*dt;    % Update connection strengths
    % Ensure W is in proper range
    W = min(W,Wmax);        % maximum connection strength
    W = max(W,Wmin);

    
    figure(pattern);
    subplot(2, 1, 1)
    imagesc(reshape(input_pattern, [N1, N1]))

    subplot(2, 1, 2)
    imagesc(reshape(rate(end, :), [N1, N1]))

end


%% Set up loop
for trials = 1: 4
    
    % randomly select one pattern
    pattern = trials;
    input_pattern = pattern_collection(:,:,pattern);
    
    % array of random number from 0 to 1
    flip = rand(Ncells, 1);
    flip_probability = 0.2;
    flipinputs = find(flip < flip_probability);
    % flip the number with 0.1 probability
    input_pattern(flipinputs) = 1-input_pattern(flipinputs);  
    
    % firing rate initialized to zero
    rate = zeros(Nt, Ncells);
    current = zeros(Ncells, 1);


    for i = 2: Nt

        % Applied current only exist in the first half of simulation
        applied_current = zeros(N1);
        if i < Nt/2
            applied_current = input_pattern*50;          
        end

        current = rate(i-1,:)*W + reshape(applied_current, [1, Ncells]);
       
        % Calculate rate
        rss = rmax./(1+exp(-del_I*(current-threshold)));         
        rate(i,:) = rate(i-1,:) + (rss-rate(i-1,:))*dt/tau_r;

    end
        
    figure(pattern+4);
    subplot(2, 1, 1)
    imagesc(reshape(input_pattern, [N1, N1]))
    
    subplot(2, 1, 2)
    imagesc(reshape(rate(end, :), [N1, N1]))
    disp(rate(end, :))

end

