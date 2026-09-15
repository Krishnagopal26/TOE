%% Quartz Stress–Strain for 7 Motions
clc; clear; close all;
disp('AVAILABLE MOTION TYPES:');
disp('  1: Uniaxial Extension');
disp('  2: Biaxial Extension');
disp('  3: Simple Shear');
disp('  4: Double Shear');
disp('  5: Rigid Body Rotation');
disp('  6: Inversion / Compression');
disp('  7: Non-Homogeneous Extension');
% ---------------- Material (Quartz) ----------------
c11 = 86.74;  c33 = 107.20; c44 = 57.94;
c12 = 6.99;   c13 = 11.91;  c14 = -17.91;
c66 = (c11 - c12) / 2;
C = [c11, c12, c13, c14, 0,   0;
     c12, c11, c13, -c14, 0,   0;
     c13, c13, c33, 0,   0,   0;
     c14, -c14, 0,   c44, 0,   0;
     0,   0,   0,   0,   c44, c14;
     0,   0,   0,   0,   c14, c66];

% ---------------- Simulation control ----------------
num_steps = 200;
time = linspace(0,1,num_steps);
I3 = eye(3);
X_loc = 1;

motion_choice = input('enetr motion')
    
    x_small = zeros(num_steps,1);
    y_small = zeros(num_steps,1);
    x_finite = zeros(num_steps,1);
    y_finite = zeros(num_steps,1);
    
    for n = 1:num_steps
        t = time(n);

        % ---------- F depending on motion ----------
        switch motion_choice
            case 1 % Uniaxial
                a = 1 + 0.2*t;
                F = [a 0 0; 0 1 0; 0 0 1];

            case 2 % Biaxial
                a = 1 + 0.1*t;
                b = 1 + 0.15*t;
                F = [a 0 0; 0 b 0; 0 0 1];

            case 3 % Simple shear
                g = 0.2*t;
                F = [1 g 0; 0 1 0; 0 0 1];

            case 4 % Double shear
                g = 0.1*t;
                F = [1 g 0; g 1 0; 0 0 1];

            case 5 % Rotation
                th = (pi/4)*t;
                c = cos(th); s = sin(th);
                F = [c s 0; -s c 0; 0 0 1];

            case 6 % Inversion
                a = 1 + 0.5*t;
                F = [-a 0 0; 0 1 0; 0 0 1];

            case 7 % Non-homogeneous
                a = 0.1*t;
                F11 = 1 + 2*a*X_loc;
                F = [F11 0 0; 0 1 0; 0 0 1];
        end

        % ---------------- SMALL STRAIN ----------------
        eps = 0.5*(F + F.') - I3;

        eps_voigt = [
            eps(1,1);
            eps(2,2);
            eps(3,3);
            2*eps(2,3);
            2*eps(1,3);
            2*eps(1,2)
        ];

        s_voigt = C * eps_voigt;

        sigma = [s_voigt(1), s_voigt(6), s_voigt(5);
                 s_voigt(6), s_voigt(2), s_voigt(4);
                 s_voigt(5), s_voigt(4), s_voigt(3)];

        dev_s = sigma - trace(sigma)/3 * I3;
        y_small(n) = sqrt(1.5 * sum(sum(dev_s.^2)));

        eps_dev = eps - trace(eps)/3 * I3;
        x_small(n) = sqrt((2/3)*sum(sum(eps_dev.^2)));

        % ---------------- FINITE STRAIN ----------------
        E = 0.5*(F.'*F - I3);

        E_voigt = [
            E(1,1);
            E(2,2);
            E(3,3);
            2*E(2,3);
            2*E(1,3);
            2*E(1,2)
        ];

        S_voigt = C * E_voigt;
        S = [S_voigt(1), S_voigt(6), S_voigt(5);
             S_voigt(6), S_voigt(2), S_voigt(4);
             S_voigt(5), S_voigt(4), S_voigt(3)];

        J = det(F);
        if abs(J) < 1e-12
            sigma_f = zeros(3);
        else
            sigma_f = (1/J) * F * S * F.';
        end
        dev_f = sigma_f - trace(sigma_f)/3 * I3;
        y_finite(n) = sqrt(1.5 * sum(sum(dev_f.^2)));

        E_dev = E - trace(E)/3 * I3;
        x_finite(n) = sqrt((2/3)*sum(sum(E_dev.^2)));

    end

    % ---------- PLOT ----------
    figure('Color','w','Position',[100 100 900 600]);
    plot(x_small, y_small, 'b', 'LineWidth', 2); hold on;
    plot(x_finite, y_finite, 'r', 'LineWidth', 2);
    xlabel('Equivalent strain'); ylabel('Equivalent stress (GPa)');
    title(sprintf('Motion %d: Stress–Strain for Quartz', motion_choice));
    legend('Small strain (linear elastic)', 'Finite strain (nonlinear geometric)');
    grid on;

