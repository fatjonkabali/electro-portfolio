clear variables;

%% Model Parameters
N = 20;                          % Sections
K = 360*2+180;                   % Time horizon (2h)
T = 10/3600;                     % 10 s timestep
tau = 0.02; l = 1.4; m = 4;      % Dynamics parameters
eta = 4; chi = 40;               % Anticipation & smoothing

L = 0.5 * ones(N,1);             % Section length [km]
vf = 120 * ones(N,1);            % Free-flow speed
rho_max = 400 * ones(N,1);       % Max density
rho_cr = 110 * ones(N,1);        % Critical density

Ir = zeros(N,1);
Ir([5,10,12]) = 1;

%% Ramp demand
dem = zeros(N,K+1);
dem(5,:)  = 700;
dem(10,:) = 800;
dem(12,:) = 650;

%% Initial Conditions
rho_iniz = 100 * ones(N,1);       % Higher starting congestion
v_iniz = vf .* (1 - (rho_iniz ./ rho_max).^l).^m;
q_iniz = rho_iniz .* v_iniz;

%% Boundary Conditions - Section 0
rho_sez0 = 90 * ones(K,1);       % Higher upstream density
v_sez0 = vf(1) * ones(K,1);
q_sez0 = rho_sez0 .* v_sez0;

%% Boundary Conditions - Section N
rho_sezfin = 80 * ones(K,1);
v_sezfin = vf(N) .* (1 - (rho_sezfin ./ rho_max(N)).^l).^m;
q_sezfin = rho_sezfin .* v_sezfin;

%% State Variables
rho = zeros(N,K+1);
v = zeros(N,K+1);
q = zeros(N,K+1);
r = dem;                          % Fixed onramp inflow

%% k = 1 initialization
rho(1,1) = rho_iniz(1) + T/L(1)*(q_sez0(1) - q_iniz(1));
v(1,1) = max(0, v_iniz(1) + T/tau*(vf(1)*(1 - (rho_iniz(1)/rho_max(1))^l)^m - v_iniz(1)));

for i=2:N
    rho(i,1) = rho_iniz(i) + T/L(i)*(q_iniz(i-1) + r(i,1)*Ir(i) - q_iniz(i));
    rho(i,1) = min(rho(i,1), rho_max(i));
    v(i,1) = max(0, v_iniz(i) + T/tau*(vf(i)*(1 - (rho_iniz(i)/rho_max(i))^l)^m - v_iniz(i)));
end

%% Simulation loop
for k = 1:K
    q(:,k) = rho(:,k) .* v(:,k);

    % Section 1
    rho(1,k+1) = rho(1,k) + T/L(1)*(q_sez0(k) - q(1,k));
    v(1,k+1) = max(0, v(1,k) + T/tau*(vf(1)*(1 - (rho(1,k)/rho_max(1))^l)^m - v(1,k)));

    % All other sections
    for i = 2:N
        inflow = q(i-1,k) + r(i,k)*Ir(i);
        rho(i,k+1) = rho(i,k) + T/L(i)*(inflow - q(i,k));
        rho(i,k+1) = min(rho(i,k+1), rho_max(i));

        % Anticipation and convection terms
        if i < N
            anticipation = (rho(i+1,k) - rho(i,k)) / (rho(i,k) + chi);
        else
            anticipation = (rho_sezfin(k) - rho(i,k)) / (rho(i,k) + chi);
        end
        convection = v(i,k) * (v(i,k) - v(i-1,k));

        v(i,k+1) = v(i,k) + (T/tau) * ...
            (vf(i)*(1 - (rho(i,k)/rho_max(i))^l)^m - v(i,k) ...
             - (T/L(i)) * convection ...
             - (eta*T/L(i)) * anticipation);
        v(i,k+1) = max(v(i,k+1), 0);
    end
end

%% Metrics
TTT = sum(sum(rho)) * L(1) * T;
avg_rho = mean(rho,2);
avg_v = mean(v,2);
%% Plots
figure; mesh(rho); title('Density (No Control)'); xlabel('Time'); ylabel('Section');
figure; mesh(v); title('Speed (No Control)'); xlabel('Time'); ylabel('Section');
figure; bar(avg_rho); title('Average Density per Section (No Control)'); xlabel('Section');
figure; bar(avg_v); title('Average Speed per Section (No Control)'); xlabel('Section');
fprintf('Total Travel Time (No Control): %.2f veh·h\\n', TTT);

ramp_sections = find(Ir == 1);  % Indices of sections with on-ramps
time = 0:K;                     % Time steps

for idx = 1:length(ramp_sections)
    i_section = ramp_sections(idx);
    
    figure;
    subplot(2,1,1);
    plot(time, rho(i_section, :), 'b-', 'LineWidth', 1.5);
    title(['Density over time at on-ramp section ', num2str(i_section)]);
    xlabel('Time step');
    ylabel('Density (veh/km)');
    grid on;

    subplot(2,1,2);
    plot(time, v(i_section, :), 'r-', 'LineWidth', 1.5);
    title(['Speed over time at on-ramp section ', num2str(i_section)]);
    xlabel('Time step');
    ylabel('Speed (km/h)');
    grid on;
end