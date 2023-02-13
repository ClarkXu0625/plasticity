% A(i)
C = 1e-9;
R = 10e6;
E = -70e-3;%V
Vth = -54e-3;%V
Vreset = -80e-3; %V
tsyn = 10e-3; %s
pR = 1;
Ibase = 2e-9;%A
G12 = 1e-6; %S
G21 = 1e-6; %S
dt = 0.1e-3; %s
Erev12 = -70e-3; %V
Erev21 = -70e-3; %V

% A(ii)
t = 1:dt:6;

%add the extra applied current to the two vectors
I1app = Ibase * ones(1, length(t));
I1app(1:100e-3/dt) = Ibase+3e-9;

I2app = Ibase*ones(1, length(t));
I2app(round(length(t)/2) : round(length(t)/2+(100e-3/dt))) = Ibase + 3e-9;

%set D1 and D2 as one
D1 = ones(1, length(t));
D2 = ones(1, length(t));


V1 = Vreset*ones(1, length(t));
V2 = Vreset*ones(1, length(t));

s1 = ones(1, length(t));
s2 = ones(1, length(t));

delta = 5e-11; %noise


for i = 2:length(t)

    if V1(i-1)> Vth
        V1(i) = Vreset;
        s1(i) = s1(i-1) + pR*D1(i-1)*(1-s1(i-1));
    else
        ds1 = -s1(i-1)/tsyn*dt;
        s1(i) = s1(i-1)+ds1;

        dV1 = (E-V1(i-1))/R + G21*s2(i-1)*(Erev21 - V1(i-1))+I1app(i-1)+delta*randn;

      
        V1(i) = V1(i-1)+ dV1*dt/C;
    end

    if V2(i-1)> Vth
        V2(i) = Vreset;
        s2(i) = s2(i-1) + pR*D1(i-1)*(1-s2(i-1));
    else
        ds2 = -s2(i-1)/tsyn*dt;
        s2(i) = s2(i-1)+ds2;

        dV2 = (E-V2(i-1))/R + G12*s1(i-1)*(Erev12 - V2(i-1)) + I2app(i-1)+delta*randn;
        V2(i) = V2(i-1)+ dV2*dt/C;
    end

end


plot(t, V2);



