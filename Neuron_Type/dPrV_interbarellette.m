%clear

%Simulation parameters
Tmin = 0;
Tmax = 0.8;
dt = 0.000005;
Tvec = Tmin:dt:Tmax;

%Neuron parameters
Gl = 3e-8;
%Gna = 1.2e-5;
%Gk = 2e-6;
Gh = 200e-9;
Gt = 0.22e-6;

El = -0.060;
Eh = -0.020;
%Ena = 0.055;
%Ek = -0.072;
Eca = 120e-3;

Cm = 1/(417e6);


Iapp = zeros(1, length(Tvec)) * 850e-12;
Vvec = zeros(1, length(Tvec));
Vvec(1) = El;

%m = zeros(1, length(Tvec));
%h = zeros(1, length(Tvec));
%n = zeros(1, length(Tvec));
m_h = zeros(1, length(Tvec));
m_t = zeros(1, length(Tvec));
h_t = zeros(1, length(Tvec));

Il = zeros(1, length(Tvec));
%Ina = zeros(1, length(Tvec));
%Ik = zeros(1, length(Tvec));
Ih = zeros(1, length(Tvec));
It = zeros(1, length(Tvec));

%m(1) = 0.0;
%h(1) = 0.0;
%n(1) = 0.0;
m_h(1) = 0.0;
m_t(1) = 0.0;
h_t(1) = 0.0;

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
    %Ina(i-1) = 0;
    %Ik(i-1) = (maxGk * n(i-1)^4 * (Ek - Vvec(i-1)));
    %Ik(i-1) = 0;
    It(i-1) = (Gt * (m_t(i-1)^2) * h_t(i-1) * (Eca - Vvec(i-1)));
    Ih(i-1) = (Gh * m_h(i-1) * (Eh - Vvec(i-1)));
    dvdt = (Il(i-1) + It(i-1) + Ih(i-1) + Iapp(i-1))/Cm;
    
    
    %Compute gating variables
    %[alpha_m, beta_m] = compute_gating_var_m(Vvec(i-1));
    %[alpha_h, beta_h] = compute_gating_var_h(Vvec(i-1));
    %[alpha_n, beta_n] = compute_gating_var_n(Vvec(i-1));
    m_h(i) = compute_m_t_inf(Vvec(i-1));
    [h_t_inf, tau_h_t] = compute_h_t_var(Vvec(i-1));
    
    %tau_m = 1./(alpha_m+beta_m);      % time constant converted from ms to sec
    %m_inf = alpha_m./(alpha_m+beta_m);

    %tau_h = 1./(alpha_h+beta_h);      % time constant converted from ms to sec
    %h_inf = alpha_h./(alpha_h+beta_h);

    %tau_n = 1./(alpha_n+beta_n);      % time constant converted from ms to sec
    %n_inf = alpha_n./(alpha_n+beta_n);
    
    [m_h_inf, tau_m_h] = compute_var_H(Vvec(i-1));
    %Compute derivative of gating variables
    dhtdt = (h_t_inf - h_t(i-1)) / tau_h_t;
    dmhdt = (m_h_inf - m_h(i-1)) / tau_m_h; 
    
    %Compute the value of each variable of the next time step
    Vvec(i) = Vvec(i-1) + dt * dvdt;
    
    %m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
    %h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
    %n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
    m_h(i) = m_h(i-1) + dt * dmhdt;
    h_t(i) = h_t(i-1) + dt * dhtdt;
    
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

%Compute T-type conductance
function m_t_inf = compute_m_t_inf(Vm)
    m_t_inf = 1 / (1 + exp(-(Vm + 0.052)/0.0074));
end

function [h_t_inf, tau_h_t] = compute_h_t_var(Vm)
    h_t_inf = compute_h_t_inf(Vm);
    tau_h_t = compute_tau_h_t(Vm);
end

function h_t_inf = compute_h_t_inf(Vm)
    h_t_inf = 1/(1 + exp(500 * (Vm + 0.0076)));
end

function tau_h_t = compute_tau_h_t(Vm)
    if(Vm < - 0.080)
        tau_h_t = 0.001 * exp(15 * (Vm + 0.467)); 
    else
        tau_h_t = 0.028 + 0.001 * exp(-(Vm + 0.022)/0.0105);
    end
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