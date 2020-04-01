clear
clc
T=0.3;
mu=[2.73 4.10, 4.78, 14.91];
R=[1.00, 0.37, 0.37, 0.00];

% function [R_dyn,mu_dyn] = Extended_SPO2IDA(mu,R,T)
    % Matlab implementation of the extension to SPO2IDA for infilled RC 
    % frames

    % Please cite as:
    % Nafeh AMB, O?Reilly GJ, Monteiro R. Simplified seismic assessment of 
    % infilled RC frame structures. Bulletin of Earthquake Engineering
    % DOI: 10.1007/s10518-019-00758-2.

    % -------------------------------------------------------------------
    % Inputs:
    % mu: Array of ductilities at points in Figure 10 (B-C-D-E)
    % R: Array of strength ratios at points in Figure 10 (B-C-D-E)
    % T: Initial period of system

    % Outputs:
    % R_dyn: Strength ratio (Eq. (7))
    % mu_dyn: Corresponding mu

    %% Define the coefficients {16%, 50%, 84%}
    % Hardening branch  
    a_alpha_1=[ 0.146, 0.8628, 1.024;...
                0.5926, 0.9235, 0.6034;...
                0.07312, 0.9195, 0.2466;...
                0.2965, 0.9632, 0.06141;...
                0.02688, 0.4745, 0.2511;...
                1.063, 0.0654, 0.0001;...
                0.3127, 0.04461, 0.07086];

    b_alpha_1=[ 0.5335, 0.7624, 0.9018;...
                0.4161, 0.5041, 0.1928;...
                0.4495, 0.1785, 0.4758;...
                0.2215, 1.022, 0.6903;...
                0.3699, 0.3253, 0.3254;...
                1.003, 0.4064, 0.939;...
                0.1462, 0.4479,0.3948];

    c_alpha_1=[ 0.03444, 0.1643, 0.6555;...
                0.3194, 0.1701, 0.1072;...
                0.01667, 0.1147, 0.1232;...
                0.1087, 0.1694, 0.05664;...
                0.0158, 0.09403, 0.07067;...
                0.646, 0.02054, 0.00132;...
                0.07181, 0.01584, 0.02287];

    a_beta_1=[  0.2008, -0.1334, 0.7182;...
                0.179, 0.3312, 0.132;...
                0.1425, 0.7985, 0.1233;...
                0.1533, 0.0001, 0.09805;...
                3.623E+12, 0.1543, 0.1429;...
                0.09451, 0.9252, 0.6547;...
                0.1964, 0.2809, 0.0001];

    b_beta_1=[ 1.093, 0.7771, 0.04151;...
                0.7169, 0.7647, 0.6058;...
                0.4876, 0.04284, 0.4904;...
                0.5709, 0.5721, 0.5448;...
                97.61, 0.4788, 0.3652;...
                0.4424, 0.8165, 0.8431;...
                0.3345, 0.3003, 0.7115];


    c_beta_1=[ 0.5405, 0.04907, 0.09018;...
                0.08836, 0.000986, 0.04845;...
                0.04956, 0.09365, 0.04392;...
                0.07256, 0.0001, 0.01778;...
                17.94, 0.105, 0.09815;...
                0.06262, 0.51, 0.7126;...
                0.09522, 0.1216, 0.0001803];
              
    % Softening branch
    a_alpha_2 = [0.03945, 0.01833, 0.009508];
    b_alpha_2 = [-0.03069, -0.01481, -0.007821];
    a_beta_2 = [1.049, 0.8237, 0.4175];
    b_beta_2 = [0.2494, 0.04082, 0.03164];
    a_gamma_2 = [-0.7326, -0.7208, -0.0375];
    b_gamma_2 = [1.116, 1.279, 1.079];
    
    % Residual plateau branch
    a_alpha_3 = [-5.075, -2.099, -0.382];
    b_alpha_3 = [7.112, 3.182, 0.6334];
    c_alpha_3 = [-1.572, -0.6989, -0.051];
    d_alpha_3 = [0.1049, 0.0481, 0.002];
    
    a_beta_3 = [16.16, 8.417, -0.027];
    b_beta_3 = [-26.5, -14.51, -1.80];
    c_beta_3 = [10.92, 6.75, 2.036];
    d_beta_3 = [1.055, 0.9061, 1.067];
    
    % Strength degradation branch
    a_alpha_4 = [-1.564, -0.5954, -0.06693];
    b_alpha_4 = [2.193, 0.817, 0.1418];
    c_alpha_4 = [-0.352, -0.09191, 0.0124];
    d_alpha_4 = [0.0149, 0.001819, -0.002012];
    a_beta_4 = [1.756, 0.7315, -0.408];
    b_beta_4 = [-8.719, -3.703, -1.333];
    c_beta_4 = [8.285, 4.391, 2.521];
    d_beta_4 = [1.198, 1.116, 1.058];
    
    %% Compute the parameters
    % For each fractile i in 16%, 50% and 84% 
    for i = 1:3
        % Hardening branch
        alpha_1(i)=sum(a_alpha_1(:,i).*exp(-((T-b_alpha_1(:,i))./c_alpha_1(:,i)).^2));
        beta_1(i)=sum(a_beta_1(:,i).*exp(-((T-b_beta_1(:,i))./c_beta_1(:,i)).^2));
        
        % Softening branch
        alpha_2(i)=a_alpha_2(i)*T+b_alpha_2(i);
        beta_2(i)=a_beta_2(i)*T+b_beta_2(i);
        gamma_2(i)=a_gamma_2(i)*T+b_gamma_2(i);

        % Residual branch
        alpha_3(i)=a_alpha_3(i)*T^3+b_alpha_3(i)*T^2+c_alpha_3(i)*T+d_alpha_3(i);
        beta_3(i)=a_beta_3(i)*T^3+b_beta_3(i)*T^2+c_beta_3(i)*T+d_beta_3(i);

        % Strength degradation branch
        alpha_4(i)=a_alpha_4(i)*T^3+b_alpha_4(i)*T^2+c_alpha_4(i)*T+d_alpha_4(i);
        beta_4(i)=a_beta_4(i)*T^3+b_beta_4(i)*T^2+c_beta_4(i)*T+d_beta_4(i);
    end
    
    %% Fit the branches
    % Initialise some arrays
    mu_dyn=linspace(1,max(mu),500);
    
    for i=1:3
        for j=1:length(mu_dyn)
            if mu_dyn(j)<=mu(1)
                % Hardening branch
                R_dyn(j,i) = alpha_1(i)*mu_dyn(j)^beta_1(i);
            elseif mu_dyn(j)<=mu(2)
                % Softening branch
                R_dyn(j,i) = alpha_2(i)*mu_dyn(j)^2+beta_2(i)*mu_dyn(j)+gamma_2(i);
            elseif mu_dyn(j)<=mu(3)
                % Residual branch
                R_dyn(j,i) = alpha_3(i)*mu_dyn(j)+beta_3(i);
            elseif mu_dyn(j)<=mu(4)
                % Strength degradation branch
                R_dyn(j,i) = alpha_4(i)*mu_dyn(j)+beta_4(i);
            end
        end
    end

    %% Check that fits are monotonically increasing
    for i=1:3
        for j=1:length(mu_dyn)-1
            if R_dyn(j+1,i)<R_dyn(j,i)
                R_dyn(j+1,i)=R_dyn(j,i);
            end
        end
    end

    %% Adjust the discontinuities between brances
    for i=1:3
        % Differences
        % Hardening Initiation
        H(i)=(1-alpha_1(i)*1^beta_1(i));
        
        % Connection Hardening-Softening
        HS(i)=(alpha_1(i)*mu(1)^beta_1(i)) - (alpha_2(i)*mu(1)^2+beta_2(i)*mu(1)+gamma_2(i));
        
        % Connection Softening-Plateau
        SP(i)=(alpha_2(i)*mu(2)^2+beta_2(i)*mu(2)+gamma_2(i)) - (alpha_3(i)*mu(2)+beta_3(i));
        
        % Connection Plateau-Degradation
        PD(i)=(alpha_3(i)*mu(3)+beta_3(i)) - (alpha_4(i)*mu(3)+beta_4(i));
        
        % Adjust
        R_dyn(mu_dyn< mu(1),i)=R_dyn(mu_dyn< mu(1),i)+H(i);
        R_dyn(mu_dyn>=mu(1),i)=R_dyn(mu_dyn>=mu(1),i)+HS(i);
        R_dyn(mu_dyn>=mu(2),i)=R_dyn(mu_dyn>=mu(2),i)+SP(i);
        R_dyn(mu_dyn>=mu(3),i)=R_dyn(mu_dyn>=mu(3),i)+PD(i);
        
    end
    
    %% Add in a flatline point
    R_dyn(end+1,:)=R_dyn(end,:);
    mu_dyn(end+1)=mu_dyn(end)+5;
    
    %% Plot the diagram
    figure; hold on; grid on; box on;
    plot(mu_dyn,R_dyn(:,1),'-.g');
    plot(mu_dyn,R_dyn(:,2),'-r');
    plot(mu_dyn,R_dyn(:,3),'-.m');
    plot([0, 1, mu], [0, 1, R],'-k');
    legend('16%','50%','84%','SPO','location','southeast');
    xlabel('Ductility \mu');
    ylabel('Strength Ratio R');

    

