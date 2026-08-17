clear;
clc;
close all;

% Material parameters
C10 = 2.50;          % kPa
D = 8.95e-7;         % 1/kPa
K1 = [1.40, 1.40];   % kPa
K2 = [1.30, 1.30];
KAPPA = [0.105, 0.105];

% Fiber orientation angles
beta_list = [0, 15, 30, 60, 90];
nBeta = length(beta_list);
% 
% % Initial specimen length
L0 = 10.0;
% 
% Stretch range
lambda = linspace(1.0, 2.0, 201);

% Storage variables
sigma11_all = zeros(nBeta,length(lambda));
lambda_y_all = zeros(nBeta,length(lambda));
lambda_z_all = zeros(nBeta,length(lambda));
J_all = zeros(nBeta,length(lambda));

% Numerical options for solving the lateral traction-free conditions
options = optimoptions('fsolve', ...
    'Display', 'off', ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12);

% Loop over all fiber orientation angles
for b = 1:nBeta

    beta = beta_list(b)*pi/180;

    % Define two symmetric fiber families
    N1 = [cos(beta); sin(beta); 0];
    N2 = [cos(beta); -sin(beta); 0];

    N = [N1, N2];

    % Storage variables for the current fiber angle
    sigma11 = zeros(size(lambda));
    lambda_y = zeros(size(lambda));
    lambda_z = zeros(size(lambda));
    J_array = zeros(size(lambda));

    % Initial guess for the lateral stretches
    lambda_y0 = 1.0;
    lambda_z0 = 1.0;

    % Solve the homogeneous uniaxial deformation problem
    for n = 1:length(lambda)

        lam = lambda(n);

        % Unknown lateral stretches
        x0 = [lambda_y0, lambda_z0];

        % Solve sigma22 = 0 and sigma33 = 0
        x = fsolve(@(x) lateral_stress(x, lam, C10, D, ...
            K1, K2, KAPPA, N), x0, options);

        ly = x(1);
        lz = x(2);

        % Store the converged lateral stretches
        lambda_y(n) = ly;
        lambda_z(n) = lz;

        % Use the current solution as the initial guess for the next point
        lambda_y0 = ly;
        lambda_z0 = lz;

        % Deformation gradient
        F = diag([lam, ly, lz]);

        % Right Cauchy-Green tensor
        C = F.'*F;

        % Determinant
        J = det(F);
        J_array(n) = J;

        % Right Cauchy-Green inverse
        Cinv = inv(C);

        % First invariant
        I1 = trace(C);

        % Modified first invariant
        Jm23 = J^(-2/3);
        I1bar = Jm23*I1;

        % Second Piola-Kirchhoff stress
        S = zeros(3,3);

        % Isotropic matrix contribution
        dI1bar_dC = Jm23*(eye(3) - I1/3*Cinv);

        S = S + 2*C10*dI1bar_dC;

        % Fiber contribution
        for a = 1:2

            Na = N(:,a);

            % Fiber structural tensor
            M = Na*Na.';

            % Fiber invariant
            I4 = Na.'*C*Na;

            % Modified fiber invariant
            I4bar = Jm23*I4;

            % Fiber strain measure
            Ealpha = KAPPA(a)*(I1bar-3) ...
                + (1-3*KAPPA(a))*(I4bar-1);

            % Tension-only fiber contribution
            if Ealpha > 0

                dI4bar_dC = Jm23*(M - I4/3*Cinv);

                dE_dC = KAPPA(a)*dI1bar_dC ...
                    + (1-3*KAPPA(a))*dI4bar_dC;

                fiberCoeff = 2*K1(a)*Ealpha ...
                    *exp(K2(a)*Ealpha^2);

                S = S + fiberCoeff*dE_dC;

            end
        end

        % Volumetric contribution
        Svol = (1/D)*(J^2-1)*Cinv;

        S = S + Svol;

        % Cauchy stress
        sigma = (1/J)*F*S*F.';

        % Store axial stress
        sigma11(n) = sigma(1,1);

    end

    % Store analytical results for the current fiber angle
    sigma11_all(b,:) = sigma11;
    lambda_y_all(b,:) = lambda_y;
    lambda_z_all(b,:) = lambda_z;
    J_all(b,:) = J_array;

end

% Calculate displacement corresponding to the applied stretch
displacement = L0*(lambda-1);

% Plot all analytical solutions
figure;
hold on;

for b = 1:nBeta
    plot(lambda, sigma11_all(b,:), 'LineWidth', 2);
end

xlabel('Stretch, \lambda');
ylabel('\sigma_{11} (kPa)');
legend('\beta=0^\circ', ...
       '\beta=15^\circ', ...
       '\beta=30^\circ', ...
       '\beta=60^\circ', ...
       '\beta=90^\circ', ...
       'Location', 'best');

title('GOH Model: Analytical Solutions for Different Fiber Orientations');
grid on;
box on;

saveas(gcf,'GOH_Analytical_All_Angles.png');
% 
% Read Abaqus numerical results
file_names = { ...
    'GOH_numerical_beta0.DAT', ...
    'GOH_numerical_beta15.DAT', ...
    'GOH_numerical_beta30.DAT', ...
    'GOH_numerical_beta60.DAT', ...
    'GOH_numerical_beta90.DAT'};

% Storage variables for Abaqus results
lambda_abaqus_all = cell(nBeta,1);
sigma11_abaqus_all = cell(nBeta,1);

% Read all Abaqus files
for b = 1:nBeta

    data = readmatrix(file_names{b});

    displacement_abaqus = data(:,1);
    sigma11_abaqus = data(:,2);

    % Convert displacement to stretch
    lambda_abaqus = 1.0 + displacement_abaqus/L0;

    lambda_abaqus_all{b} = lambda_abaqus;
    sigma11_abaqus_all{b} = sigma11_abaqus;

end

% Plot all Abaqus numerical results
figure;
hold on;

for b = 1:nBeta

    plot(lambda_abaqus_all{b}, ...
         sigma11_abaqus_all{b}, ...
         'o', ...
         'MarkerSize', 5, ...
         'LineWidth', 1.2);

end

xlabel('Stretch, \lambda');
ylabel('\sigma_{11} (kPa)');
legend('\beta=0^\circ', ...
       '\beta=15^\circ', ...
       '\beta=30^\circ', ...
       '\beta=60^\circ', ...
       '\beta=90^\circ', ...
       'Location', 'best');

title('GOH Model: ABAQUS Numerical Results');
grid on;
box on;

saveas(gcf,'GOH_Abaqus_All_Angles.png');

% Plot analytical and Abaqus results together
figure;
hold on;

for b = 1:nBeta

    plot(lambda, ...
         sigma11_all(b,:), ...
         'LineWidth', 2);

    plot(lambda_abaqus_all{b}, ...
         sigma11_abaqus_all{b}, ...
         'o', ...
         'MarkerSize', 4, ...
         'LineWidth', 1.0);

end

xlabel('Stretch, \lambda');
ylabel('\sigma_{11} (kPa)');

legend('\beta=0^\circ analytical', ...
       '\beta=0^\circ Abaqus', ...
       '\beta=15^\circ analytical', ...
       '\beta=15^\circ Abaqus', ...
       '\beta=30^\circ analytical', ...
       '\beta=30^\circ Abaqus', ...
       '\beta=60^\circ analytical', ...
       '\beta=60^\circ Abaqus', ...
       '\beta=90^\circ analytical', ...
       '\beta=90^\circ Abaqus', ...
       'Location', 'best');

title('GOH Model: Analytical Solution vs. ABAQUS');

grid on;
box on;

saveas(gcf,'GOH_Abaqus_vs_Analytical_All_Angles.png');

% Calculate errors for each fiber orientation
fprintf('\n');
fprintf('==============================================\n');
fprintf('GOH UMAT VALIDATION\n');
fprintf('==============================================\n');

for b = 1:nBeta

    % Abaqus data
    lambda_abaqus = lambda_abaqus_all{b};
    sigma11_abaqus = sigma11_abaqus_all{b};

    % Interpolate analytical solution at Abaqus data points
    sigma11_analytical = interp1( ...
        lambda, ...
        sigma11_all(b,:), ...
        lambda_abaqus, ...
        'pchip');

    % Calculate absolute error
    absolute_error = abs( ...
        sigma11_abaqus - sigma11_analytical);

    % Calculate relative error
    relative_error = absolute_error ./ ...
        max(abs(sigma11_analytical),1e-12);

    % Remove possible NaN values
    valid = ~isnan(relative_error);

    absolute_error = absolute_error(valid);
    relative_error = relative_error(valid);

    fprintf('\n');
    fprintf('Fiber angle beta = %d degrees\n',beta_list(b));
    fprintf('Maximum absolute error = %.6e kPa\n', ...
        max(absolute_error));
    fprintf('Maximum relative error = %.6e\n', ...
        max(relative_error));
    fprintf('Mean relative error = %.6e\n', ...
        mean(relative_error));

end

fprintf('\n');
fprintf('==============================================\n');


function residual = lateral_stress(x, lambda, C10, D, ...
    K1, K2, KAPPA, N)

    % Lateral stretches
    ly = x(1);
    lz = x(2);

    % Prevent nonphysical negative stretches
    if ly <= 0 || lz <= 0
        residual = [1e10; 1e10];
        return;
    end

    % Deformation gradient
    F = diag([lambda, ly, lz]);

    % Right Cauchy-Green tensor
    C = F.'*F;

    % Determinant
    J = det(F);

    % Inverse of C
    Cinv = inv(C);

    % First invariant
    I1 = trace(C);

    % Modified first invariant
    Jm23 = J^(-2/3);
    I1bar = Jm23*I1;

    % Second Piola-Kirchhoff stress
    S = zeros(3,3);

    % Isotropic matrix contribution
    dI1bar_dC = Jm23*(eye(3) - I1/3*Cinv);

    S = S + 2*C10*dI1bar_dC;

    % Fiber contribution
    for a = 1:2

        Na = N(:,a);

        % Fiber structural tensor
        M = Na*Na.';

        % Fiber invariant
        I4 = Na.'*C*Na;

        % Modified fiber invariant
        I4bar = Jm23*I4;

        % Fiber strain measure
        Ealpha = KAPPA(a)*(I1bar-3) ...
            + (1-3*KAPPA(a))*(I4bar-1);

        % Fiber activation
        if Ealpha > 0

            dI4bar_dC = Jm23*(M - I4/3*Cinv);

            dE_dC = KAPPA(a)*dI1bar_dC ...
                + (1-3*KAPPA(a))*dI4bar_dC;

            fiberCoeff = 2*K1(a)*Ealpha ...
                *exp(K2(a)*Ealpha^2);

            S = S + fiberCoeff*dE_dC;

        end
    end

    % Volumetric contribution
    Svol = (1/D)*(J^2-1)*Cinv;

    S = S + Svol;

    % Cauchy stress
    sigma = (1/J)*F*S*F.';

    % Traction-free lateral surfaces
    residual = [sigma(2,2); sigma(3,3)];

end