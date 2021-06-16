%clear

%Simulation parameters
Tmin = 0;
Tmax = 0.8;
dt = 0.000005;
Tvec = Tmin:dt:Tmax;

%Pulse parameters
delay = 20e-3;
%Neuron parameters
Gl = 3e-8;
maxGna = 1.2e-5;
maxGk = 2e-6;
maxGh = 200e-9;
maxGa = 4.77e-6;
El = -0.017;
Eh = -0.020;
Ena = 0.055;
Ek = -0.072;
Ea = -0.075;
Cm = 0.1e-9;

Iapp = zeros(1, length(Tvec)) * 850e-12;
Vvec = zeros(1, length(Tvec));
Vvec(1) = El;

m = zeros(1, length(Tvec));
h = zeros(1, length(Tvec));
n = zeros(1, length(Tvec));
m_h = zeros(1, length(Tvec));
a = zeros(1, length(Tvec));
b = zeros(1, length(Tvec));

Il = zeros(1, length(Tvec));
Ina = zeros(1, length(Tvec));
Ik = zeros(1, length(Tvec));
Ia = zeros(1, length(Tvec));
Ih = zeros(1, length(Tvec));

m(1) = 0.0;
h(1) = 0.0;
n(1) = 0.0;
m_h(1) = 0.0;
a(1) = 0.0;
b(1) = 0.0;

%Setup current pulses
for i = 1:length(Tvec)
    if(i*dt >= 100e-3 && i*dt <= 600e-3)
        Iapp(i) = 850e-12;
    end
end

for i = 2:length(Tvec)-1   
    %Compute derivative of membrane potential
    Il(i-1) = (Gl * (El - Vvec(i-1)));
    %Ina(i-1) = (maxGna * m(i-1)^3 * h(i-1) * (Ena - Vvec(i-1)));
    Ina(i-1) = 0;
    %Ik(i-1) = (maxGk * n(i-1)^4 * (Ek - Vvec(i-1)));
    Ik(i-1) = 0;
    Ia(i-1) = (maxGa * a(i-1)^3 * b(i-1) * (Ea - Vvec(i-1)));
    Ih(i-1) = (maxGh * m_h(i-1) * (Eh - Vvec(i-1)));
    dvdt = (Il(i-1) + Ina(i-1) + Ik(i-1) + Ia(i-1) + Ih(i-1) + Iapp(i-1))/Cm;
    
    
    %Compute gating variables
    [alpha_m, beta_m] = compute_gating_var_m(Vvec(i-1));
    [alpha_h, beta_h] = compute_gating_var_h(Vvec(i-1));
    [alpha_n, beta_n] = compute_gating_var_n(Vvec(i-1));
    
    tau_m = 1./(alpha_m+beta_m);      % time constant converted from ms to sec
    m_inf = alpha_m./(alpha_m+beta_m);

    tau_h = 1./(alpha_h+beta_h);      % time constant converted from ms to sec
    h_inf = alpha_h./(alpha_h+beta_h);

    tau_n = 1./(alpha_n+beta_n);      % time constant converted from ms to sec
    n_inf = alpha_n./(alpha_n+beta_n);
    
    [m_h_inf, tau_m_h] = compute_var_H(Vvec(i-1));
    %Compute derivative of gating variables
    dmhdt = (m_h_inf - m_h(i-1))/tau_m_h; 
    [dadt, dbdt] = compute_gating_var_A(Vvec(i-1), a(i-1), b(i-1));
    
    %Compute the value of each variable of the next time step
    Vvec(i) = Vvec(i-1) + dt * dvdt;
    
    m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
    h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
    n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
    m_h(i) = m_h(i-1) + dt * dmhdt;
    a(i) = a(i-1) + dt * dadt;
    b(i) = b(i-1) + dt * dbdt;
    
end

figure(10);
subplot(2, 1, 1);
plot(Tvec, Iapp);
subplot(2, 1, 2);
plot(Tvec, Vvec);

%Compute gating variable for the opening of the sodium channel
function [alpha_m, beta_m] = compute_gating_var_m(Vm)
    alpha_m = compute_alpha_m(Vm);
    beta_m = compute_beta_m(Vm);
end

function alpha_m = compute_alpha_m(Vm)
    alpha_m = 3.80e5*(Vm+0.0297)./(1-exp(-100*(Vm+0.0297)));
end

function beta_m = compute_beta_m(Vm)
    beta_m = 1.52e4*exp(-55.6*(Vm+0.0547));
end

%Compute gating variable for the closing of the sodium channel
function [alpha_h, beta_h] = compute_gating_var_h(Vm)
    alpha_h = compute_alpha_h(Vm);
    beta_h = compute_beta_h(Vm);
end

function alpha_h = compute_alpha_h(Vm)
    alpha_h = 266*exp(-50*(Vm+0.048));
end

function beta_h = compute_beta_h(Vm)
    beta_h = 3800./(1+exp(-100*(Vm+0.018)));
end

%Compute gating variable for potassium
function [alpha_n, beta_n] = compute_gating_var_n(Vm)
    alpha_n = compute_alpha_n(Vm);
    beta_n = compute_beta_n(Vm); 
end

function alpha_n = compute_alpha_n(Vm)
    alpha_n = 2e4*(Vm+0.0457)./(1-exp(-100*(Vm+0.0457)));
end

function beta_n = compute_beta_n(Vm)
    beta_n = 250*exp(-12.5*(Vm+0.0557));
end

%Compute A-type conductance
function [dadt, dbdt] = compute_gating_var_A(Vm, a, b)
    dadt = compute_dadt(Vm, a);
    dbdt = compute_dbdt(Vm, b);
end

function dadt = compute_dadt(Vm, a)
    tau_a = compute_tau_a(Vm);
    a_inf = compute_a_inf(Vm);
    dadt = (a_inf - a) / tau_a; 
end

function tau_a = compute_tau_a(Vm)
    tau_a = 0.3632e-3 + 1.158e-3./(1+exp(49.7*(Vm+0.05596)));
end

function a_inf = compute_a_inf(Vm)
    a_inf = (0.0761*exp(31.4*(Vm+0.09422))./(1+exp(34.6*(Vm+0.00117)))).^(1/3.0);
end

function dbdt = compute_dbdt(Vm, b)
    tau_b = compute_tau_b(Vm);
    b_inf = compute_b_inf(Vm);
    dbdt = (b_inf - b) / tau_b;
end

function tau_b = compute_tau_b(Vm)
    tau_b = 1.24e-3 + 2.678e-3./(1+exp(62.4*(Vm+0.050)));
end

function b_inf = compute_b_inf(Vm)
    b_inf = (1./(1+exp(68.8*(Vm+0.0533)))).^4;
end

%Compute H-type conductance
function [m_h_inf, tau_m_h] = compute_var_H(Vm)
    m_h_inf = compute_m_h_inf(Vm);
    tau_m_h = compute_tau_m_h(Vm);
end

function m_h_inf = compute_m_h_inf(Vm)
    m_h_inf = 1 / (1 + exp((Vm + 0.070) / 0.00873));
end

function tau_m_h = compute_tau_m_h(Vm)
    tau_m_h = 0.272 + (1.499 / (1 + exp(-(Vm + 0.0422)/0.00873)));
end