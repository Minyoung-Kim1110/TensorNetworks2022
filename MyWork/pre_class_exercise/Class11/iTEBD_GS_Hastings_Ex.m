function [Lambda,Bs,Eiter] = iTEBD_GS_Hastings_Ex (Lambda,Bs,H,Nkeep,taus)
% < Description >
%
% [Lambda,Bs,Eiter] = iTEBD_GS_Hastings_Ex (Lambda,Bs,H,Nkeep,taus)
%
% The iTEBD (infinite time-evolving block decimation) method to find the
% ground state of an infinite one-dimensional system, by applying imaginary
% time evolutions. This function implements M. Hastings' version proposed
% in Hastings2009 [M. Hastings, J. Math. Phys. 50, 095207 (2009)]. Here we
% consider a unit cell of two sites.
%
% < Input >
% Lambda : [1 x 2 cell] Lambda{1} and Lambda{2} contain the singular values
%       at odd and even bond, respectively, as column vectors. An odd
%       (even) bond sits just on the right of an odd (even) site.
% Bs : [1 x 2 cell] Bs{..} are "supposedly" right-normalized tensors that
%       are given by the contraction of Gamma*Lambda type. Bs{1} and Bs{2}
%       act on an odd and even site, respectively. In terms of the
%       conventions used in iTEBD_GS_Vidal (which implements Vidal's
%       original iTEBD), Bs{k} corresponds to Gamma{k}*diag(Lambda{k}).
%       Note that, if we initialize the tensors with random elements, they
%       are of course not right-normalized. But as the imaginary time
%       evolution goes, they will converge to the right-normalized forms.
% H : [tensor] Two-site interaction Hamiltonian. Its leg convention is as
%       below:
%
%    2         4        [ 2 (4) is to be contracted with the third leg
%    ^         ^          of Bs{1} (Bs{2}) ]
%    |   ...   |
%   [     H     ]
%    |   ...   |
%    ^         ^
%    1         3
%
% Nkeep : [integer] Maximum bond dimension.
% taus : [numeric] Vector of imaginary time step sizes. Each "outer"
%       iteration consists of two imaginary-time evolutions, the first for
%       odd bonds and the second for even bonds. Both time evolutions
%       within the m-th outer iteration take the times size taus(m).
%
% < Output >
% Lambda, Bs : [1 x 2 cells each] Cell arrays of Lambda and approximately
%       right-normalized tensors, repectively, after the imaginary time
%       evolution.
% Eiter : [(numel(taus) x 2 x 2 matrix] Eiter(m,n,k) is the measured energy
%       for an odd (k = 1) or even (k = 2) bond after odd (n = 1) or
%       even (n = 2) bonds are updated, at the m-th "outer" iteration.
%
% Written by S.Lee (Jun.18,2017)
% Updated by S.Lee (Jun.19,2017)
% Updated by S.Lee (Jun.04,2019): Revised for Sose 2019.
% Updated by S.Lee (Sep.30,2022): Revised for the course at SNU.


tobj = tic2;

Lambda = Lambda(:);
Bs = Bs(:);
Nstep = numel(taus);
ldim = size(H,1); % local space dimension
Skeep = 1e-8;

% % % check the integrity of input
if any([numel(Lambda) numel(Bs)] ~= 2)
    error('ERR: # of sites per unit cell should be 2.');
end

if ndims(H) > 4
    error('ERR: H should be rank-4.');
elseif any(ldim ~= [size(H,2) size(H,3) size(H,4)])
    error('ERR: All the legs of H should have the same size.');
end

for it = (1:2)
    if ~isvector(Lambda{it})
        error(['ERR: Lambda{',sprintf('%i',it),'} should be vector.']);
    elseif numel(Lambda{it}) ~= size(Bs{it},2)
        error(['ERR: Dimensions for Lambda{',sprintf('%i',it),'} and Bs{', ...
            sprintf('%i',it),'} do not match.']);
    elseif numel(Lambda{mod(it,2)+1}) ~= size(Bs{it},1)
        error(['ERR: Dimensions for Lambda{',sprintf('%i',mod(it,2)+1), ...
            '} and Bs{',sprintf('%i',it),'} do not match.']);
    elseif size(Bs{it},3) ~= ldim
        error(['ERR: The third leg of Bs{',sprintf('%i',mod(it)), ...
            '} should be of size equal to the leg of H.']);
    end
end
% % % 

% show information
disptime(['iTEBD ground state search (Hastings'' version): Nkeep = ',sprintf('%i',Nkeep), ...
    ', # of imag. time steps = ',sprintf('%.4g',Nstep)]);

% energy expectation value at each step
Eiter = zeros(Nstep,2,2);

% diagonalize the Hamiltonian to exponentiate
Hmat = reshape(permute(H,[1 3 2 4]),(ldim^2)*[1 1]); % matrix representation
[VH,DH] = eig((Hmat+Hmat')/2);
DH = diag(DH);

for it1 = (1:Nstep)
    % exponentiate the matrix representation of Hamiltonian
    expH = VH*diag(exp(-taus(it1)*DH))*VH';

    % reshape matrix -> rank-4 tensor
    expH = reshape(expH,ldim*ones(1,4)); % [ ldim, ldim, ldim, ldim], [1, 3, 2, 4]
    
    for it2 = (1:2) % 1 (2): update odd (even) bonds
        % % % % TODO (start) % % % %

        % It is also important to normalize the rank-4 tensor (without the
        % Lambda tensor on its left) with the same normalization factor
        % used to normalize the singular value vector; otherwise the result
        % can diverge!
        % T = contract(diag(Lambda{3-it2}), 2, 2, Bs{it2}, 3, 1);
        T = contract(Bs{it2}, 3, 2, Bs{3-it2}, 3, 1, [1 2 4 3]);
        Phi = contract(expH, 4, [3, 4], T, 4, [2, 3], [3 1 2 4]); % why 3 1 2 4?

        eHT = contract(diag(Lambda{3-it2}), 2, 2, Phi, 4, 1);

        [UT,Lambda{it2},Bs{3-it2}] = svdTr(eHT,4,[1 2],Nkeep,Skeep);
        UT = permute(UT,[1 3 2]);
        Bs{3-it2} = permute(Bs{3-it2},[1 3 2]);

        normalization_factor = norm(Lambda{it2});
        Lambda{it2} = Lambda{it2}/normalization_factor; 
        Phi = Phi/normalization_factor ; 
        Bs{it2} = contract(Phi, 4, [3, 4], conj(Bs{3-it2}), 3, [3, 2], [1 3 2]);
        


        % Measure energy per bond; consider the following ket:
        % Lambda{2}*Bs{1}*Bs{2}*Bs{1}
        %
        % The physical legs of the first two B tensors, i.e., Bs{1} and
        % Bs{2}, contract to the legs of the two-site gate that updates
        % an odd bond; those of the latter B tensors, i.e., Bs{2} and
        % Bs{1}, updates an even bond.
        
        ketA = contract(diag(Lambda{2}), 2, 2, Bs{1}, 3, 1); % Lambda{2}*Bs{1}
        ketB = Bs{2}; % Bs{2}
        ketC = Bs{1}; % Bs{1}
        

        % contract the nearest-neighbor interaction term for an odd bond,
        % i.e., its 1st and 2nd legs act on an odd site
        T = contract(ketA,3,2,ketB,3,1,[1 2 4 3]);
        Ho = contract(H,4,[2 4],T,4,[2 3],[3 1 2 4]);
        Ho = contract(conj(T),4,(1:3),Ho,4,(1:3));
        Ho = contract(Ho,2,2,ketC,3,1);
        Ho = contract(conj(ketC),3,(1:3),Ho,3,(1:3));

        % contract the nearest-neighbor interaction term for an even bond,
        % i.e., its 1st and 2nd legs act on an even site
        T = contract(ketB,3,2,ketC,3,1,[1 2 4 3]);
        He = contract(H,4,[2 4],T,4,[2 3],[3 1 2 4]);
        He = contract(conj(T),4,(2:4),He,4,(2:4));
        He = contract(ketA,3,2,He,2,2);
        He = contract(conj(ketA),3,[1 3 2],He,3,(1:3));

        % compute the squared norm of the rank-5 ket tensor as a denominator for normalization
        ovl = updateLeft([],[],ketA,[],[],ketA);
        ovl = updateLeft(ovl,2,ketB,[],[],ketB);
        ovl = contract(ovl,2,2,ketC,3,1);
        ovl = contract(conj(ketC),3,(1:3),ovl,3,(1:3));

        % assign normalized energy values
        Eiter(it1,it2,:) = [Ho He]/ovl;

        % % % % TODO (end) % % % %
    end

    if (mod(it1,500) == 0) && (it1 < Nstep)
        disptime(['#',sprintf('%i/%i',[it1 Nstep]),', E = ',sprintf('%.8g',mean(Eiter(it1,end,:),3))]);
    end
end

disptime(['#',sprintf('%i/%i',[it1 Nstep]),', E = ',sprintf('%.8g',mean(Eiter(end,end,:),3))]);

% check performance
toc2(tobj,'-v');
chkmem;


end