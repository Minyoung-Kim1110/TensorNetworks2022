clear

% iTEBD parameters
tau_ini = 1; % initial imaginary time step size
tau_fin = 1e-6; % final imaginary time step size
Nstep = 2e3;
% Local operators
[S,I] = getLocalSpace('Spin',1);
% Heisenberg interaction as two-site gate S*S'
H = contract(S,3,3,permute(conj(S),[2 1 3]),3,3);

% Initialize with random Lambda and Gamma
Lambda_init = cell(1,2);
Gamma_init = cell(1,2);
for Nkeep = (30:10:50)

    for itn = (1:numel(Lambda_init))
        Lambda_init{itn} = rand(Nkeep,1);
        Gamma_init{itn} = rand(Nkeep,Nkeep,size(I,2));
    end
    % Vidal ground state search
    taus = tau_ini*((tau_fin/tau_ini).^linspace(0,1,Nstep));
    [Lambda,Gamma,Eiters] = iTEBD_GS_Vidal(Lambda_init,Gamma_init,H,Nkeep,taus);
    
    % orthonormalize 
    new_Gamma = contract(Gamma{1}, 3, 2, diag(Lambda{1}), 2, 1, [1 3 2]);
    new_Gamma =contract(new_Gamma, 3, 2, Gamma{2}, 3, 1, [1 3 2 4]);
    new_Gamma = reshape(new_Gamma, size(new_Gamma, 1), size(new_Gamma, 2), size(new_Gamma, 3)*size(new_Gamma, 4));
    new_Lambda = Lambda{2};
    
    [Lambda, Gamma] = ortho_Orus_Ex(new_Lambda, new_Gamma);
    
    % define operator
    Sz_o = kron(S(:, :, 2), I);
    Sz_e = kron(I, S(:, :, 2));
    % define state 
    state = contract(diag(Lambda), 2, 1, Gamma, 3, 1);
    
    max_n = 30 ;
    correlation= zeros(2, max_n);
    
    % calculate middle ( state - identity - state)
    Ts = cell(1, max_n);
    for itN =(1:max_n)
        if itN==1
            Ts{itN} = contract(state, 3, 3, conj(state), 3, 3, [1 3 2 4]);
      
        else 
            Ts{itN} = contract(Ts{itN-1}, 4, [3, 4], Ts{1}, 4, [1, 2]);
        end 
    end
    
    % calculate spin-spin correlation 
    for itN = (10:max_n)
        % odd 
        T_o = updateLeft([], [], state, Sz_o, 2, state);
        T_e = updateLeft([], [], state, Sz_e, 2, state);
        norm = updateLeft([],[], state, [],[], state);
        
        if itN>1
            T_o = contract(T_o, 2, [1, 2], Ts{itN-1}, 4, [2, 1], [2, 1]);
            T_e = contract(T_e, 2, [1, 2], Ts{itN-1}, 4, [2, 1], [2, 1]);
            norm = contract(norm, 2, [1,2], Ts{itN-1}, 4, [2, 1], [2, 1]);
        end
        
        identity_3= reshape(eye(size(T_o, 1)), size(T_o, 1), size(T_o, 1), 1); 
        
        T_o = updateLeft(T_o, 2, state, Sz_o, 2, state);
        T_o = contract(T_o, 2, [1,2], identity_3,3, [1,2]); % zip
        T_e = updateLeft(T_e, 2, state, Sz_e, 2, state);
        T_e = contract(T_e, 2, [1,2], identity_3,3, [1,2]); % zip

        norm = updateLeft(norm, 2, state, [], [], state);
        norm = contract(norm, 2, [1,2], identity_3, 3, [1,2]); % zip
        correlation(1, itN) = T_o/norm;
        correlation(2, itN) = T_e/norm;
        
    end
    x = 2* (10:max_n)';
    yo = log(abs(correlation(1, 10:end)))';
    ye = log(abs(correlation(2, 10:end)))';
    
    bo = x\yo;
    be = x\ye;

    xi_o = - 1/bo;
    xi_e = - 1/be;
    
    disp([newline 'check xi for odd bond convergence: Nkeep = ',sprintf('%i',Nkeep), ', xi= ',sprintf('%.4g',xi_o) ]);
    disp(['check xi for even bond convergence: Nkeep = ',sprintf('%i',Nkeep), ', xi= ',sprintf('%.4g',xi_e) newline]);
end

figure;
hold on;
plot(2 * (2:max_n), abs(correlation(1, 2:end)),'LineWidth',1);
plot(2 * (2:max_n), abs(correlation(2, 2:end)),'LineWidth',1);

hold off;
set(gca, 'YScale', 'log', 'LineWidth', 1, 'FontSize', 13);
xlabel('n');
ylabel('correlation');
legend(['odd bond'],['even bond']);
grid on;


