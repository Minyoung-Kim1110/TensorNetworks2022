clear
N =3; % number of spins
%t = exp(1i*(1:N-1)); % hopping amplitudes
t = ones(1, N-1);
[F,Z,I] = getLocalSpace('Fermion');
H = 0; % initialize Hamiltonian
Aprev = 1; % identity for the vacuum

for itN = (1:N)
 % rank-3 identity tensor for the current iteration
 Anow = getIdentity(Aprev,2,I,2,[1 3 2]);
 % contract the Hamiltonian up to the last iteration with
 % ket and bra tensors
 H = updateLeft(H,2,Anow,[],[],Anow);
  if itN > 1
 ZF = contract(Z,2,2,F,3,1);
 Hhop = (-t(itN-1))*updateLeft(Fprev,3,Anow, ...
 permute(conj(ZF),[2 1 3]),3,Anow);
 % hopping from the last site to the current site
 H = H + Hhop + Hhop';
 end
 % update operator for the next iteration
 Fprev = updateLeft([],[],Anow,F,3,Anow);
 Aprev = Anow; % to be used for the next iteration
end

Es = sort(eig((H+H')/2),'ascend');
disp(Es(1:7).');
