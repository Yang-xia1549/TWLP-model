function [X, G, Core, Out,B] = TWLP(F, Omega, opts)
if isfield(opts, 'tol');         tol   = opts.tol;              end
if isfield(opts, 'maxit');       maxit = opts.maxit;            end
if isfield(opts, 'rho');         rho   = opts.rho;              end
if isfield(opts, 'R');           max_R = opts.R;                end
if isfield(opts,'alpha');        alpha = opts.alpha;  end
if isfield(opts,'miu');        miu = opts.miu;  end
if isfield(opts,'omiga');        omiga = opts.omiga;  end
if isfield(opts,'lamda');        lamda = opts.lamda;  end
if isfield(opts,'p');        p = opts.p;  end
if isfield(opts,'eta');        eta = opts.eta;  end
if isfield(opts,'Xtrue');     Xtrue=opts.Xtrue;end
%
num_padarray = 0;% if num_padarray=0, no_rank_inc_adaptive
R = max(max_R-num_padarray, 2);
Ndim = ndims(F); 
Nway = size(F);
X = F;

Factors_dims = factor_dims(Nway, R);
Max_Factors_dims = factor_dims(Nway, max_R);
rng('default')
G = cell(Ndim,1);
A = cell(Ndim,1);
Asum=cell(Ndim,1);
for i=1:Ndim
    G{i}=rand(Factors_dims(i,:));
    A{i}=rand(Factors_dims(i,:));
    Asum{i}=rand(Factors_dims(i,:));
end
Core = rand(R(2,:));
B = rand(R(2,:));
C_old = Core;
Out.RSE = [];

A_sum=cell(Ndim,4);
for i=1:Ndim
    for j=1:4
        A_sum{i,j}=zeros(size(A{i}));
    end
end
r_change=0.001;
for k = 1:maxit
    X_old = X;
    % Update G_k, k=1,2,...,N.
    for num = 1:Ndim
        GCrest = tenremat(circ_tnprod_rest(G,Core,num), Ndim);  % Q is the right part of the relation equation
        TempA  = tenmat_sb(X,num)*GCrest'+rho*my_Unfold(G{num},size(G{num}),2)+2*lamda*my_Unfold(A{num},size(A{num}),2);
        TempB  = (GCrest*GCrest')+(2*lamda+rho)*eye(size(GCrest,1),size(GCrest,1));
        G{num} = my_Fold(TempA/TempB,size(G{num}),2);
    end
    % Update the core tensor C.
    Girest = tnreshapemat(order_tnprod_rest(G), Ndim);
     if k<=200 || numel(Core)>numel(C_old) || mod(k,20)==0
   
        TempC = reshape(X,[1,prod(Nway)])*Girest'+rho*reshape(Core,[1,numel(Core)])+omiga*reshape(B,[1,numel(B)]);
        TempD = (Girest*Girest')+(omiga+rho)*eye(size(Girest,1),size(Girest,1));
        TempE = TempC/TempD;
        Core = reshape(TempE,size(Core));
           %% update B
    TempF=(omiga*Core+rho*B)/(omiga+rho);
    B   =  solve_Lp_w( TempF, miu/(omiga+rho), p );
     end

  %% % Update X
    X = (rho*X_old+reshape(TempE*Girest,Nway))/(1+rho);
    X(Omega)= F(Omega); 
  

   %% update A
 
    for num = 1:Ndim
            tn=0;
            Asum{num}=Asum{num}*0;
    for i=1:4
    if alpha(num,i)>0

           A_sum{num,i} = fold(Pro2TraceNorm(unfold((lamda*G{num}+rho*A_sum{num,i})/(lamda+rho), size(G{num}),i,2), eta*alpha(num,i)/(lamda+rho)), size(G{num}), i, 2);
    tn=tn+1;
    Asum{num}=Asum{num}+A_sum{num,i};
    end    
    end
     if tn==0
        tn=1;
     end
 
     A{num}=Asum{num}/tn;
   
    end
    %% check the convergence
    rse = norm(X(:)-X_old(:)) / norm(X_old(:));
    Out.RSE = [Out.RSE, rse];

    if k == 1 || mod(k, 20) == 0
        fprintf('inc_TW-TC: iter = %d   RSE=%.10f   \n', k, rse);
    end

    if rse < tol 
        break;
    end
    
    C_old = Core;
    rank_inc = double(Factors_dims<Max_Factors_dims);
    if rse<r_change && sum(rank_inc(:))~=0
        [G, Core,A,Asum,B,A_sum] = rank_inc_adaptive(G, Core,A, Asum,B,rank_inc, Ndim,A_sum);
        Factors_dims  = Factors_dims+rank_inc;
        r_change = r_change*0.1;
    end
end
function T = fold(Tk,dims,k,d)
        dims = circshift(dims,d-k); 
        T = reshape(Tk,dims);        
        T = shiftdim(T,4+d-k);      
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cover the orignal tensor T into matrix T_<k(d)>
function  Tk = unfold(T,dims,k,d)
     T = shiftdim(T,4+k-d);     
     dims = circshift(dims,d-k); 
     Tk = reshape(T,prod(dims(1:d)),[]);    
end


end



function [G, Core,A,Asum,B,A_sum] = rank_inc_adaptive(G, Core,A,Asum,B, rank_inc, N,A_sum)
    % increase the estimated rank until max_R
    for j = 1:N
    G{j} = padarray(G{j}, rank_inc(j,:), rand(1), 'post');
    A{j} = padarray(A{j}, rank_inc(j,:), rand(1), 'post');
    Asum{j} = padarray(Asum{j}, rank_inc(j,:), rand(1), 'post');
    for a=1:4
    A_sum{j,a}=padarray(A_sum{j,a}, rank_inc(j,:), rand(1), 'post');
    end
    end
    Core = padarray(Core, rank_inc(:,3), rand(1), 'post');
    B= padarray(B, rank_inc(:,3), rand(1), 'post');
end
