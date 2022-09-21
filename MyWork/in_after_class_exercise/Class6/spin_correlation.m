AKLT = zeros(2,2,3);
% local spin S_z = +1
AKLT(1,2,1) = sqrt(2/3);
% local spin S_z = 0
AKLT(1,1,2) = -1/sqrt(3);
AKLT(2,2,2) = +1/sqrt(3);
% local spin S_z = -1
AKLT(2,1,3) = -sqrt(2/3);

[S,I] = getLocalSpace('Spin',1);
Sz = S(:,:,2); % spin-z

L = 50; % number of sites

result = zeros(2, 2, L);
analytic_result = zeros(2,2,L);
for alpha = (1:2)
    for beta = (1:2)

        for n =(1:L)
            M = cell(1,L); % MPS
            M(:) = {AKLT};

            M{1} = M{1}(alpha,:,:);
            M{end} = M{end}(:,beta,:);
            
            [M, Sv] = canonForm(M, n ,[],0);
            
            T = updateLeft([], [], M{n}, Sz, 2, M{n});
            T =  diag(Sv)/norm(Sv)* T * diag(Sv)/norm(Sv);
            T = updateLeft(T, 2, M{n+1}, Sz, 2, M{n+1});
            result(alpha, beta, n) = sum(diag(T));
            analytic_result(alpha, beta, n)  =analytic_spin_spin_correlation(alpha, beta, L, n)
        end 
    end 
end 
