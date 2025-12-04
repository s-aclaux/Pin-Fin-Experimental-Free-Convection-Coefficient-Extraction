clear;
clc;
figs = findall(0, 'Type', 'figure');
for k = 1:length(figs)
    clf(figs(k));
end

%% DATA %%

%$% Copper (Square) Data %$%
T_copS_avg = [64.77, 61.40, 60.49, 58.16, 54.02, 53.51, 53.75];
TC_copS_dist = [0.00, 3.06, 5.50, 9.82, 17.15, 22.07, 28.07] * 10^(-2);
T_copS_amb = 22.32;

%$% Copper (Circle) Data %$%
T_copC_avg = [75.90, 73.26, 71.88, 69.52, 67.97];
TC_copC_dist = [0.00, 3.09, 5.56, 9.80, 17.15] * 10^(-2);
T_copC_amb = 22.57;

%$% Stainless Steel Data %$%
T_steel_avg = [95.70, 61.15, 46.65, 31.88, 24.04, 22.94, 22.75];
TC_steel_dist = [0.00, 3.10, 5.44, 9.81, 17.22, 22.07, 28.45] * 10^(-2);
T_steel_amb = 21.87;

%$% Aluminium Data %$%
T_alum_avg = [85.77, 79.88, 74.27, 66.52, 59.24, 56.99, 55.05];
TC_alum_dist = [0.00, 4.35, 5.41, 9.72, 17.15, 22.07, 28.07] * 10^(-2);
T_alum_amb = 21.50;

%% PLOTS %%

%$% Copper (Square) Plot %$%
figure(1);
hold on;
plot(TC_copS_dist, T_copS_avg, '-r', 'LineWidth', 1.5);
yline(T_copS_amb, '--k', 'LineWidth', 1.5);
hold off;
xlabel("Distance Along Fin from Base to Tip (m)"); ylabel("Pin Fin Temperature (C)");
title({"Square Copper Pin Fin Steady-State","Temperature Distribution"});
legend("Fin Temp","Air Temp",'Location','northeast');

%$% Copper (Circle) Plot %$%
figure(2);
hold on;
plot(TC_copC_dist, T_copC_avg, '-r', 'LineWidth', 1.5);
yline(T_copC_amb, '--k', 'LineWidth', 1.5);
hold off;
xlabel("Distance Along Fin from Base to Tip (m)"); ylabel("Pin Fin Temperature (C)");
title({"Circle Copper Pin Fin Steady-State","Temperature Distribution"});
legend("Fin Temp","Air Temp",'Location','northeast');

%$% Stainless Steel Plot %$%
figure(3);
hold on;
plot(TC_steel_dist, T_steel_avg, '-r', 'LineWidth', 1.5);
yline(T_steel_amb, '--k', 'LineWidth', 1.5);
hold off;
xlabel("Distance Along Fin from Base to Tip (m)"); ylabel("Pin Fin Temperature (C)");
title({"Stainless Steel Pin Fin Steady-State","Temperature Distribution"});
legend("Fin Temp","Air Temp",'Location','northeast');

%$% Aluminium) Plot %$%
figure(4);
hold on;
plot(TC_alum_dist, T_alum_avg, '-r', 'LineWidth', 1.5);
yline(T_alum_amb, '--k', 'LineWidth', 1.5);
hold off;
xlabel("Distance Along Fin from Base to Tip (m)"); ylabel("Pin Fin Temperature (C)");
title({"Aluminium Pin Fin Steady-State","Temperature Distribution"});
legend("Fin Temp","Air Temp",'Location','northeast');

%% CONVECTION COEFFICIENT CALCS %%

m = @(h, P, k, A) sqrt((h*P) / (k*A));
theta = @(x, h, L, m, k) (cosh(m *(L - x)) + (h / (m*k)) * sinh(m*(L - x))) / (cosh(m*L) + (h / (m*k)) * sinh(m*L));

%$% Copper (Square) h Calculation %$%
h_guess = 11; % W / (m^2 * K)
k_cop = 388; % W / (m * K)
A_copS = (1.27 * 10^(-2))^2; % m^2
L_copS = 28.50 * 10^(-2); % m
P_copS = 4*(1.27*10^(-2)); % m

S_copS = 10;
i = 0;
dS = 1;
while dS > 0.00001 && i < 100
    m_copS = m(h_guess, P_copS, k_cop, A_copS);
    theta_copS = @(n) theta(TC_copS_dist(n), h_guess, L_copS, m_copS, k_cop);
    
    % Thermocouple 1
    T_copS_1 = T_copS_amb + (T_copS_avg(1) - T_copS_amb) * theta_copS(2);
    err_copS_1 = T_copS_1 - T_copS_avg(2);
    % Thermocouple 2
    T_copS_2 = T_copS_amb + (T_copS_avg(1) - T_copS_amb) * theta_copS(3);
    err_copS_2 = T_copS_2 - T_copS_avg(3);    
    % Thermocouple 3
    T_copS_3 = T_copS_amb + (T_copS_avg(1) - T_copS_amb) * theta_copS(4);
    err_copS_3 = T_copS_3 - T_copS_avg(4);    
    % Thermocouple 4
    T_copS_4 = T_copS_amb + (T_copS_avg(1) - T_copS_amb) * theta_copS(5);
    err_copS_4 = T_copS_4 - T_copS_avg(5);    
    % Thermocouple 5
    T_copS_5 = T_copS_amb + (T_copS_avg(1) - T_copS_amb) * theta_copS(6);
    err_copS_5 = T_copS_5 - T_copS_avg(6);    
    % Thermocouple 6
    T_copS_6 = T_copS_amb + (T_copS_avg(1) - T_copS_amb) * theta_copS(7);
    err_copS_6 = T_copS_6 - T_copS_avg(7);
    
    S_copS_new = err_copS_1^2 + err_copS_2^2 + err_copS_3^2 + err_copS_4^2 + err_copS_5^2 + err_copS_6^2;
    i = i + 1;
    dS = abs(S_copS_new - S_copS);
    if dS > 1
        h_guess = h_guess - 0.1*sqrt(S_copS_new / 6);
    elseif dS > 0.01
        h_guess = h_guess - 0.05*sqrt(S_copS_new / 6);
    elseif dS > 0.005
        h_guess = h_guess - 0.005*sqrt(S_copS_new / 6);
    end
    S_copS = S_copS_new;
end
h_copS = h_guess

%$% Copper (Circle) h Calculation %$%
h_guess = 15;
A_copC = pi *(1.27 * 10^(-2) / 2)^2; % m^2
L_copC = 17.20 * 10^(-2); % m
P_copC = 2*pi*(1.27*10^(-2) / 2); % m

S_copC = 10;
i = 0;
dS = 1;
while dS > 0.00001 && i < 100
    m_copC = m(h_guess, P_copC, k_cop, A_copC);
    theta_copC = @(n) theta(TC_copC_dist(n), h_guess, L_copC, m_copC, k_cop);
    
    % Thermocouple 1
    T_copC_1 = T_copC_amb + (T_copC_avg(1) - T_copC_amb) * theta_copC(2);
    err_copC_1 = T_copC_1 - T_copC_avg(2);
    % Thermocouple 2
    T_copC_2 = T_copC_amb + (T_copC_avg(1) - T_copC_amb) * theta_copC(3);
    err_copC_2 = T_copC_2 - T_copC_avg(3);    
    % Thermocouple 3
    T_copC_3 = T_copC_amb + (T_copC_avg(1) - T_copC_amb) * theta_copC(4);
    err_copC_3 = T_copC_3 - T_copC_avg(4);    
    % Thermocouple 4
    T_copC_4 = T_copC_amb + (T_copC_avg(1) - T_copC_amb) * theta_copC(5);
    err_copC_4 = T_copC_4 - T_copC_avg(5);    
    
    S_copC_new = err_copC_1^2 + err_copC_2^2 + err_copC_3^2 + err_copC_4^2;
    i = i + 1;
    dS = abs(S_copC_new - S_copC);
    if dS > 1
        h_guess = h_guess - 0.1*sqrt(S_copC_new / 6);
    elseif dS > 0.05
        h_guess = h_guess - 0.004*sqrt(S_copC_new / 6);
    elseif dS > 0.005
        h_guess = h_guess - 0.00004*sqrt(S_copC_new / 6);
    end
    S_copC = S_copC_new;
end
h_copC = h_guess

%$% Stainless Steel h Calculation %$%
h_guess = 20; % W / (m^2 * K)
k_steel = 16; % W / (m * K)
A_steel = pi * (0.95 * 10^(-2) / 2)^2; % m^2
L_steel = 28.50 * 10^(-2); % m
P_steel = 2*pi*(0.95 * 10^(-2) / 2); % m

S_steel = 5;
i = 0;
dS = 1;
while dS > 0.00001 && i < 100
    m_steel = m(h_guess, P_steel, k_steel, A_steel);
    theta_steel = @(n) theta(TC_steel_dist(n), h_guess, L_steel, m_steel, k_steel);
    
    % Thermocouple 1
    T_steel_1 = T_steel_amb + (T_steel_avg(1) - T_steel_amb) * theta_steel(2);
    err_steel_1 = T_steel_1 - T_steel_avg(2);
    % Thermocouple 2
    T_steel_2 = T_steel_amb + (T_steel_avg(1) - T_steel_amb) * theta_steel(3);
    err_steel_2 = T_steel_2 - T_steel_avg(3);    
    % Thermocouple 3
    T_steel_3 = T_steel_amb + (T_steel_avg(1) - T_steel_amb) * theta_steel(4);
    err_steel_3 = T_steel_3 - T_steel_avg(4);    
    % Thermocouple 4
    T_steel_4 = T_steel_amb + (T_steel_avg(1) - T_steel_amb) * theta_steel(5);
    err_steel_4 = T_steel_4 - T_steel_avg(5);    
    % Thermocouple 5
    T_steel_5 = T_steel_amb + (T_steel_avg(1) - T_steel_amb) * theta_steel(6);
    err_steel_5 = T_steel_5 - T_steel_avg(6);    
    % Thermocouple 6
    T_steel_6 = T_steel_amb + (T_steel_avg(1) - T_steel_amb) * theta_steel(7);
    err_steel_6 = T_steel_6 - T_steel_avg(7);
    
    S_steel_new = err_steel_1^2 + err_steel_2^2 + err_steel_3^2 + err_steel_4^2 + err_steel_5^2 + err_steel_6^2;
    i = i + 1;
    dS = abs(S_steel_new - S_steel);
    if dS > 1
        h_guess = h_guess - 0.5*sqrt(S_steel_new / 6);
    elseif dS > 0.01
        h_guess = h_guess - 0.2*sqrt(S_steel_new / 6);
    elseif dS > 0.005
        h_guess = h_guess - 0.1*sqrt(S_steel_new / 6);
    end
    S_steel = S_steel_new;
end
h_steel = h_guess

%$% Aluminium h Calculation %$%
h_guess = 15; % W / (m^2 * K)
k_alum = 167; % W / (m * K)
A_alum = pi * (1.27 * 10^(-2) / 2)^2; % m^2
L_alum = 28.50 * 10^(-2); % m
P_alum = 2*pi*(1.27 * 10^(-2) / 2); % m

S_alum = 10;
i = 0;
dS = 1;
while dS > 0.00001 && i < 100
    m_alum = m(h_guess, P_alum, k_alum, A_alum);
    theta_alum = @(n) theta(TC_alum_dist(n), h_guess, L_alum, m_alum, k_alum);
    
    % Thermocouple 1
    T_alum_1 = T_alum_amb + (T_alum_avg(1) - T_alum_amb) * theta_alum(2);
    err_alum_1 = T_alum_1 - T_alum_avg(2);
    % Thermocouple 2
    T_alum_2 = T_alum_amb + (T_alum_avg(1) - T_alum_amb) * theta_alum(3);
    err_alum_2 = T_alum_2 - T_alum_avg(3);    
    % Thermocouple 3
    T_alum_3 = T_alum_amb + (T_alum_avg(1) - T_alum_amb) * theta_alum(4);
    err_alum_3 = T_alum_3 - T_alum_avg(4);    
    % Thermocouple 4
    T_alum_4 = T_alum_amb + (T_alum_avg(1) - T_alum_amb) * theta_alum(5);
    err_alum_4 = T_alum_4 - T_alum_avg(5);    
    % Thermocouple 5
    T_alum_5 = T_alum_amb + (T_alum_avg(1) - T_alum_amb) * theta_alum(6);
    err_alum_5 = T_alum_5 - T_alum_avg(6);    
    % Thermocouple 6
    T_alum_6 = T_alum_amb + (T_alum_avg(1) - T_alum_amb) * theta_alum(7);
    err_alum_6 = T_alum_6 - T_alum_avg(7);
    
    S_alum_new = err_alum_1^2 + err_alum_2^2 + err_alum_3^2 + err_alum_4^2 + err_alum_5^2 + err_alum_6^2;
    i = i + 1;
    dS = abs(S_alum_new - S_alum);
    if dS > 1
        h_guess = h_guess - 0.1*sqrt(S_alum_new / 6);
    elseif dS > 0.1
        h_guess = h_guess - 0.01*sqrt(S_alum_new / 6);
    elseif dS > 0.01
        h_guess = h_guess - 0.001*sqrt(S_alum_new / 6);
    end
    S_alum = S_alum_new;
end
h_alum = h_guess

%% HEAT FLOW CALCS %%

q_f = @(h, P, k, A, T_b, T_infty, m, L) sqrt(h*P*k*A)*(T_b - T_infty) * (sinh(m*L) + (h/(m*k))*cosh(m*L)) / (cosh(m*L) + (h/(m*k))*sinh(m*L));

%$% Copper (Square) q_f Calculation %$%
m_copS = m(h_copS, P_copS, k_cop, A_copS);
q_copS = q_f(h_copS, P_copS, k_cop, A_copS, T_copS_avg(1), T_copS_amb, m_copS, L_copS)

%$% Copper (Circle) q_f Calculation %$%
m_copC = m(h_copC, P_copC, k_cop, A_copC);
q_copC = q_f(h_copC, P_copC, k_cop, A_copC, T_copC_avg(1), T_copC_amb, m_copC, L_copC)

%$% Stainless Steel q_f Calculation %$%
m_steel = m(h_steel, P_steel, k_steel, A_steel);
q_steel = q_f(h_steel, P_steel, k_steel, A_steel, T_steel_avg(1), T_steel_amb, m_steel, L_steel)

%$% Aluminium q_f Calculation %$%
m_alum = m(h_alum, P_alum, k_alum, A_alum);
q_alum = q_f(h_alum, P_alum, k_alum, A_alum, T_alum_avg(1), T_alum_amb, m_alum, L_alum)

%% BOUNDARY TEMPERATURE VECTORS %%

T_b_v = [T_copS_avg(1);
         T_copC_avg(1);
         T_steel_avg(1);
         T_alum_avg(1)]

T_infty_v = [T_copS_amb;
             T_copC_amb;
             T_steel_amb;
             T_alum_amb]

T_t_v = [T_copS_avg(7);
         T_copC_avg(5);
         T_steel_avg(7);
         T_alum_avg(7)]
T_b_minus_T_t = T_b_v - T_t_v