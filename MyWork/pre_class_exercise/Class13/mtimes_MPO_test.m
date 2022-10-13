clear
% simple test 
Nkeep = 10;

B = cell(1, 3);
B{1} = rand(4, 4, 1, 3);
B{2} = rand(4, 4, 3, 3);
B{3} = rand(4, 4, 3, 1);

A = cell(1, 3); 
A{1} = rand(4, 4, 1, 3);
A{2} = rand(4, 4, 3, 3);
A{3} = rand(4, 4, 3, 1);

% mtimes_MPO is accurate up to normalized tensor 
%for it = (1:3)
%    mat = A{it};
%    A{it} = A{it}/norm(mat(:));
%    mat = B{it};
%    B{it} = B{it}/norm(mat(:));
%end 


true = cell(1, 3);
for it = (1:3)
    true{it} = contract(A{it}, 4, 1, B{it}, 4, 2 ,[ 4 1 2 5 3 6]);
    true{it} = reshape(true{it}, size(true{it}, 1),size(true{it}, 2), ...
        size(true{it}, 3)*size(true{it}, 4), size(true{it}, 5)*size(true{it}, 6));
    
    if it ~=3 && numel(size(true{it})) ~= 4 
        error("error in reshaping contracted MPO");
    end 
end 


for it = (2:-1:1)
    mat = contract(true{it}, 4, 4, true{it+1}, 4, 3); 
    [u, s, vd] = svdTr(mat, 6, [1 2 3], Nkeep, []);
    true{it} = contract(u, 4, 4, diag(s), 2, 1);
    true{it+1} = permute(vd, [2 3 1 4]);
    %size(true{it})
    %size(true{it+1})
end
mat = true{1};
mat = contract(mat, 4, 4, true{2}, 4, 3, [1 4 2 5 3 6]); 
mat = contract(mat, 6, 6, true{3}, 4, 3, [ 1 2 6 3 4 7 5 8]);

%tobj = tic2;
Nsweep = 4; 
errors = zeros(1, 4);
for Nsweep = (1:4)
    cal = mtimes_MPO_Ex(B,A,Nkeep,Nsweep);
    mat2 = cal{1};
    mat2= contract(mat2, 4, 4, cal{2}, 4, 3, [1 4 2 5 3 6]); 
    mat2= contract(mat2, 6, 6, cal{3}, 4, 3, [ 1 2 6 3 4 7 5 8]);
   errors(Nsweep) = max(abs(mat - mat2),[], 'all')
end 
%toc2(tobj,'-v');
%chkmem;


