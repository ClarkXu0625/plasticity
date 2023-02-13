%step (a)
dt = 1e-4;
t = 0:dt:4;

%step (b)
firingRate = zeros(1, length(t));
firingRate(1:1/dt) = 20;
firingRate(1/dt+1:2/dt) = 100;
firingRate(2/dt+1:3/dt) = 10;
firingRate(3/dt+1:4/dt) = 20;

% step (c)
Spike = zeros(1, length(t));
randoms = rand(1,length(t));
for i =1 : length(firingRate)
    x = randoms(i);
    randoms(i) = x;
    if x < firingRate(i)*dt
        Spike(i) = 1;
    end
end

% step (d)
Gsyn = zeros(1, length(t));
Tconstant = 0.1; %decay back to zero by 100 ms (0.1s)
decay = 0;
dG = 1e-9; %constant increment
for i = 2 : length(Gsyn)
    if decay<Tconstant
        Gsyn(i) = Gsyn(i-1)+dG;
        decay = decay+1e-4;
    else
        decay = 0;
    end
end

%step (e)
P0 = 0.5;
tD = 0.25;
D = ones(1, length(t));
for i = 2:length(D)
    if Spike(i) == 0
        dD = (1-D(i-1))/tD*dt;
        D(i) = D(i-1) +dD;
    else
        D(i) = D(i-1)*P0;
    end
    
end



%step (f)
dtc = 0.1; %decay time constant
dGmax = 2e-9; %S

G = zeros(1, length(t));
d = 0;
for i = 2: length(G)    
    if d<dtc
        dG = dGmax*P0*D(i-1);
        G(i) = G(i-1) + dG;
        d = d+1e-4;
    else
        G(i) = 0;
        d = 0;
    end    
end

%step (g)
ffac = 0.25;
tF = 0.25;

F = ones(1, length(t));
P = zeros(1, length(t));
P(:) = 0.2;
for i = 2:length(F)
    if Spike(i) == 1
        F(i) = F(i-1)+ ffac * (1/P(i) - F(i-1));
    else
        dF = (1-F(i-1))/tF*dt;
        F(i) = F(i-1)+ dF;
    end
end



%step (h)
D2 = ones(1, length(t));%new synaptic depression vector
for i = 2:length(D2)
    if Spike(i) ==1
        D2(i) = D2(i-1) - P(i) * F(i-1) * D2(i-1);
    else
        dD2 = (1-D2(i-1))/tD*dt;
        D2(i) = D2(i-1) +dD2;
    end
end


%Step (i)
G3 = zeros(1, length(t));
dGmax = 4e-9;
P0 = 0.2;
d = 0;
for i = 2: length(G3)
    if d<dtc        
        dG3 = dGmax *P0 *F(i-1)* D2(i-1);
        G3(i) = G3(i-1) + dG3;
        d = d+dt;
    else
        d = 0;
    end

end

plot(t, D2);