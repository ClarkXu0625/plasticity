clear;
%% Parameter
V_th = -50e-3;
V_reset = -80e-3;
sigma = 1e-3;%mV
tau = 3e-3;%mS
E_L = -0.070;%mV
E_I = -65e-3;%mV
E_E = 0;%mV
G_L = 50e-12;%pS

W_EE = 25e-9; %nS
W_EI = 4e-9; %nS
W_IE = 800e-9; %nS

%% Question 1 parameters
t_E = 2e-3;%ms
t_I = 5e-3;%ms
alph = 0.2;
dt = 1e-4;
t = 0:dt:2.5;


%% vectors set up

s1_E = zeros(1, length(t));
s2_I = zeros(1, length(t));

% inhibitory and excitatory conductance
G1_E = zeros(1, length(t));
G2_E = zeros(1, length(t));
G1_I = zeros(1, length(t));
G2_I = zeros(1, length(t));

% excitatory input conductance
G1_in = 1e-9; %nS
G2_in = 0; %nS

Vss_1 = zeros(1, length(t));
Vss_2 = zeros(1, length(t));
f1 = zeros(1, length(t));
f2 = zeros(1, length(t));

%firing rate
r1 = zeros(1, length(t));
r2 = zeros(1, length(t));

%% 
for i = 2: length(t)
    %get s
    ds1_E = dt* (-s1_E(i-1)/t_E+alph*r1(i-1)*(1-s1_E(i-1)));
    s1_E(i) = s1_E(i-1)+ds1_E;
    
    ds2_I = dt* (-s2_I(i-1)/t_I+alph*r2(i-1)*(1-s2_I(i-1)));
    s2_I(i) = s2_I(i-1)+ds2_I;

    %get conductance
    G1_E(i) = W_EE*s1_E(i)+G1_in;
    G1_I(i) = W_IE*s2_I(i);
    G2_E(i) = W_EI*s1_E(i)+G2_in;

    %steady state membrane potential
    Vss_1(i) = (G_L*E_L+G1_I(i)*E_I+G1_E(i)*E_E)/ ...
        (G_L+G1_I(i)+G1_E(i));
    Vss_2(i) = (G_L*E_L+G2_I(i)*E_I+G2_E(i)*E_E)/(G_L+G2_I(i)+G2_E(i));
    
    %check if divider is zero, i.e. Vss=Vth±delta && delta=0
    if Vss_1(i) ==V_th
        f1(i) = sigma/(tau*(Vss_1(i)-V_reset)*alph);
    else
        f1(i) = (Vss_1(i)-V_th)/(tau*(V_th-V_reset)*(1-exp(-(Vss_1(i)-V_th)/sigma)));
    end

    if Vss_2(i) ==V_th
        f2(i) = sigma/(tau*(Vss_2(i)-V_reset)*alph);
    else
        f2(i) = (Vss_2(i)-V_th)/(tau*(V_th-V_reset)*(1-exp(-(Vss_2(i)-V_th)/sigma)));
    end

    %get firing rate
    dr1 = (-r1(i-1)+f1(i-1))*dt/tau;
    r1(1,i) = r1(i-1)+dr1;
    dr2 = (-r2(i-1)+f2(i-1))*dt/tau;
    r2(1,i) = r2(i-1)+dr2;
    
    
end

%% plot the firing rate
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

%plot(t, r1);
hold on
%plot(t, r2);



%% Question 2 couting oscillation
maxThreshold_1 = 54.10;
minThreshold_1 = 0.05;
count_1 = 0;
startTime_1 = 0;
endTime_1 = 0;

maxThreshold_2 = 134.30;
minThreshold_2 = 0.080;
count_2 = 0;
startTime_2 = 0;
endTime_2 = 0;

%status is 0 when the firing rate is dropping, is 1 when it is rising
status_1 = 0;  
status_2 = 0;

%% start counting oscillations from 0.5 second
for i = 0.5/dt:length(t)
    % when the firing rate is dropping to the minimum threshold
    if r1(i)<=minThreshold_1 && status_1==0
        status_1 =1;
        %record the start time
        if count_1 ==0
            startTime_1 = i*dt;
        end
    % when the firing rate rises and reached the maximum threshold
    elseif r1(i)>=maxThreshold_1 && status_1==1
        status_1 =0;
        count_1 = count_1+1;
        endTime_1 = i*dt;    
    end
    
    %repeat the counting process of unit 1 to 2
    if r2(i)<=minThreshold_2 && status_2==0
        status_2 =1;
        if count_2 == 0
            startTime_2 = i*dt;
        end
    elseif r2(i)>=maxThreshold_2 && status_2==1
        status_2 = 0;
        count_2 = count_2+1;
        endTime_2 = i*dt;    
    end        
end

period1 = count_1/(endTime_1-startTime_1);
period2 = count_2/(endTime_2-startTime_2);
x1 = ['Frequency of unit 1 is ', num2str(period1), 'Hz'];
x2 = ['Frequency of unit 2 is ', num2str(period2), 'Hz'];
disp(count_1);
disp(x1);
disp(x2);

%% question 3

f1 = r1(0.5/dt+1:length(t));
f2 = r2(0.5/dt+1:length(t));

Vsine = 0.5:dt:2.5;
Vcosine = 0.5:dt:2.5;

f = 0:0.5:100;  % different firing rates from 0 to 100 Hz
disp(size(Vsine));
disp(size(f1));

A = zeros(1, length(f));
B = zeros(1, length(f));
P = zeros(1, length(f));
% loop through every frequency values
for i = 1:length(f)
    
    for j = 1: length(Vsine)
        Vsine(j) = sin(2*pi*f(i)*f1(j));
        Vcosine(j) = cos(2*pi*f(i)*f1(j));        
    end
    A(i) = sum(Vsine.*f1)/length(f1);
    B(i) = sum(Vcosine.*f1)/length(f1);
    P(i) = A(i)^2 + B(i)^2;

end

plot(f(2:length(f)), P(2:length(P))); % plot the P and discard the first value

max = 0;
frequency = 0;
for i = 2:length(P)
    if P(i)>max
        max = P(i);
        frequency = f(i);
    end
end
disp(frequency)

xlabel('f (Hz)')
ylabel('P(f)')