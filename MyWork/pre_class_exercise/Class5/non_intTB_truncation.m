function Gs = non_intTB_truncation(t, Nkeep)
    N = numel(t)+1;
    [F,Z,I] = getLocalSpace('Fermion');
    H = zeros(2); % initialize Hamiltonian
    Aprev = getIdentity(1, 2, I, 2, [1,3,2]); 
    Fprev = updateLeft([],[],Aprev,F,3,Aprev);
    U = cell(1, N); % left uniary 
    Ds = cell(1, N);
    U{1}=1;
    
    for itN = (2:N) 
        % rank-3 identity tensor for the current iteration
        % It spans Hilbert space 
        Anow = getIdentity(Aprev, 2, I,2,[1 3 2]);
        % contract the Hamiltonian up to the last iteration with ket and bra
        % tensors -> Hamiltonian with new basis 
        H = updateLeft(H,2,Anow,[],[],Anow);
        ZF = contract(Z,2,2,F,3,1);

        % hopping from the last site to the current site
        Hhop = (-t(itN-1))*updateLeft(Fprev,3,Anow,permute(conj(ZF),[2 1 3]),3,Anow);
        H = H + Hhop + Hhop';
        
        % get energies 
        [V, D] = eig((H+H')/2);
        
        if Nkeep<numel(diag(D))
            % sort by energy(eigen values of H) 
            [Es, idx] = sort(diag(D), 'ascend');
            D = D(idx, idx); V = V(:, idx);
            energy_criteria = Es(Nkeep); % criteria = Nkeep th energy 
            % check degeneracy
            degeneracy = 1;
            while abs(energy_criteria - Es(Nkeep+degeneracy))<eps
                degeneracy = degeneracy+1;
            end
            Nkeep_revised = Nkeep+degeneracy-1 ; % update Nkeep 
            % Truncate 
            V(:, Nkeep_revised+1:end) = []; 

            Anow =contract(Anow, 3, 2, V, 2, 1, [1,3,2]);
            H = D(1:Nkeep_revised, 1:Nkeep_revised);
            Ds{itN}=D(1:Nkeep_revised, 1:Nkeep_revised);
        else 
            H = D;
            Ds{itN}=D;
            U{itN}=V;
            Anow = contract(Anow, 3, 2, V, 2, 1, [1,3,2]);
        end
        
        
        % update operator in new basis for the next iteration
        Fprev = updateLeft([],[],Anow,F,3,Anow);
        Aprev = Anow;
        
    end
    
    Es = sort(eig((H+H')/2),'ascend');
    Gs = Es(1);
end

