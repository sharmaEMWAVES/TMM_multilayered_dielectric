% This calculates the reflection and transmission coefficients of
% stratified dielectric/MTM medium (N-layers) using TMM approach for normal incidence. The
% incident and transmission regions are considered as free-space
% (epsilon_r=1, mu_r=1) and the wave propoagation is in
% +z-direction.

N = 2; % Number of layers
epsilon_r = [-7,1]; % Relative permittivity of each layer
mu_r = [1, 1]; % Relative permeability of each layer
d = [5e-3, 5e-3]; % Thickness of each layer in meters
f_start=9e9; %Minimum frequency
f_stop=11e9; % Maximum frequency
freq = linspace(f_start, f_stop, 1001); % Frequency range with 1001 points

% Constants
c = 3e8; % Speed of light in vacuum


% Transmission and reflection coefficient array
T = zeros(size(freq));
R = zeros(size(freq));

% Loop over frequencies to calculate transmission coefficient
for i = 1:length(freq)
    f = freq(i);
    omega = 2 * pi * f;
    k0 = omega./c; %free-space wave number
    Z0 = 377; %characteristic impedance of free-space
    
    % Initialize overall transfer matrix
    M = eye(2);
    
    % Loop through each layer
    for n = 1:N
        % Propagation constants
        if epsilon_r(n) > 0 && mu_r(n) > 0 %DPS
            k = k0 * sqrt(epsilon_r(n) * mu_r(n)); 
            Z = Z0 * sqrt(mu_r(n) / epsilon_r(n)); % Characteristic impedance of the medium
            
        elseif epsilon_r(n) < 0 && mu_r(n) > 0 %ENG
            k = 1i * k0 * sqrt(abs(epsilon_r(n))) * sqrt(mu_r(n));
            Z = -1i * Z0 * sqrt(mu_r(n) / abs(epsilon_r(n)));
            
        elseif epsilon_r(n) > 0 && mu_r(n) < 0 %MNG
            k = 1i * k0 * sqrt(epsilon_r(n)) * sqrt(abs(mu_r(n)));
            Z = 1i * Z0 * sqrt(abs(mu_r(n)) / epsilon_r(n));
            
        elseif epsilon_r(n) < 0 && mu_r(n) < 0 %DNG
            k = -k0 * sqrt(abs(epsilon_r(n))) * sqrt(abs(mu_r(n)));
            Z = Z0 * sqrt(abs(mu_r(n)) / abs(epsilon_r(n)));
            
        end
        
        % Transfer matrix for current layer
        M_layer = [cos(k * d(n)), (1i / Z) * sin(k * d(n));
                   1i * Z * sin(k * d(n)), cos(k * d(n))];
        
        % Update overall transfer matrix
        M = M_layer * M;
    end
    
    % Transmission coefficient calculation
    T(i) = 2 / (M(1,1) + M(1,2)*Z0 + M(2,1)/Z0 + M(2,2));
    %Reflection coefficient calculation
    %R(i)=(M(1,1)+Z0*M(1,2)-(M(2,1)/Z0)-M(2,2))./(M(1,1)+Z0*M(1,2)+(M(2,1)/Z0)+M(2,2));
end

% Calculate transmission and reflection coefficient magnitude
t = abs(T).^2;
%r =  abs(R).^2;

% Calculate transmission and reflection coefficient in dB
T_dB = round(10 * log10(t),5);

% Plot the transmission coefficient
figure(1);
plot(freq * 1e-9, T_dB,'LineWidth',3);
xlabel('Frequency (GHz)', 'FontSize',22,'FontWeight','bold')
ylabel('Transmission Magnitude (dB)', 'FontSize',22,'FontWeight','bold')
set(gca, 'YGrid', 'on', 'XGrid', 'on','FontSize',22,'FontWeight','bold')
exportgraphics(gcf,'transmission.png','Resolution',300)




