
function [S] = SSSI_L2p_ADMM(B,L,VertConn,sparseweight,varargin)
%% Description: Extended sources reconstruction based on structured sparse regularization and L2p-norm
% Input:
%         B(d_b x T):               M/EEG Measurement
%         L(d_b x d_s):             Leadfiled Matrix
%         VertConn:                 Cortex Connectivity Condition
%         sparseweight:             sparseweight for the sparse variation
% Output:
%         S:                        Estimated Sources

[Nsensor,Nsource] = size(L);
T = size(B,2);
V = VariationEdge(VertConn);
tol = 1e-3;
QUIET = 1;
rou_update = 0;
cost = 0;
rou1 = 1e-10;%6e-11'
rou2 = rou1/100;
p=0.7;

if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'transform'
                transform = varargin{i+1};
            case 'tol'
                tol = varargin{i+1};
            case 'lam'
                lam = varargin{i+1};
            case 'alp'
                alp = varargin{i+1};
            case 'p'
                p = varargin{i+1};
        end
    end
end       
  Edge = VariationEdge(VertConn);      
if strcmp(transform, 'Variation')
    V = VariationEdge(VertConn);
elseif strcmp(transform, 'Laplacian')
    NVertConn = sum(VertConn,2);
    V = bsxfun(@minus,spdiags(ones(Nsource,1),0,Nsource,Nsource),bsxfun(@times,bsxfun(@rdivide,VertConn,NVertConn),0.95*ones(Nsource,1)));
elseif strcmp(transform, 'Laplacian+Variation')
    NVertConn = sum(VertConn,2);
    V = bsxfun(@minus,spdiags(ones(Nsource,1),0,Nsource,Nsource),bsxfun(@times,bsxfun(@rdivide,VertConn,NVertConn),0.95*ones(Nsource,1)));
    V = [opts.laplacian*V;VariationEdge(VertConn)];
elseif strcmp(transform,'Sparse+Laplacian')
    NVertConn = sum(VertConn,2);
    V = bsxfun(@minus,spdiags(ones(Nsource,1),0,Nsource,Nsource),bsxfun(@times,bsxfun(@rdivide,VertConn,NVertConn),0.95*ones(Nsource,1)));
    V = [sparseweight*sparse(1:Nsource,1:Nsource,1);V];
elseif strcmp(transform,'Sparse+Variation')
    V = [sparseweight*sparse(1:Nsource,1:Nsource,1);VariationEdge(VertConn)];
elseif strcmp(transform, 'Sparse')
    V = sparse(1:Nsource,1:Nsource,1);
end

ADMM_iter = 600;

% Initial values
TMNE = MNE(B,[],L,[],'MNE');
S_MNE = TMNE*B;
S = S_MNE;

U = V*S;  U_old = U;
Z = S; Z_old = Z;
W = zeros(size(V,1),T);
C = zeros(size(V,2),T);

rou1_old = rou1;
rou2_old = rou2;
S_old = S;
tmp = sqrt(sum((V*L'*B).^2,2));
lam = (1/(max(tmp)*0.01))^(-p)
alp = 0.01
kesi = max(U, [], 2)*1e-3;
kesi2 = max(S, [], 2)*1e-3;
% kesi = std(U,0,2).^2;
% kesi2 = std(S,0,2).^2;

% approximation of inverse (Woodbury inversion lemma)
Lambda_MAX = eigs(V'*V+rou2/rou1*speye(Nsource),1,'lm');
LLt = L*L'; LtB = L'*B;
mu = 0.9/(rou1*Lambda_MAX);
Precision = mu*speye(Nsource) - mu^2*L'/(eye(Nsensor) + mu*LLt)*L;

tic 
%% ADMM
alpha = 0.6;
for iter_ADMM = 1:ADMM_iter
%----------------------------S update----------------------------%
    S = Precision*(LtB + (1/mu)*S - rou1*V'*(V*S - U + W) - rou2*(S - Z + C));
%----------------------------U update----------------------------%
    VS = V * S;
    VS_hat = alpha*VS + (1-alpha)*U_old;
    
    w = 1;
    U_inter_old = zeros(size(V,1),T);
    for iter_U = 1:10%U_iter
        U = repmat(1./(rou1/2 +2*w),1,size(U,2)).*(VS_hat + W)*rou1/2;
        rr = sqrt(kesi+sum(U.^2,2));
        %             w = lam^(-p)*rr.^(p-2);
        w = p*lam*rr.^(p-2);
        w(w<max(w)*1e-6) = max(w)*1e-6;
        if  norm((U - U_inter_old),'fro')/norm(U,'fro') < 1e-6 || sum(w) == 0
            break;
        end
        U_inter_old = U;
    end
%----------------------------Z update----------------------------%
    S_hat = alpha*S + (1-alpha)*Z_old;

    w = 1;
    Z_inter_old = zeros(size(S,1),T);
    for iter_Z = 1:10
        Z = repmat(1./(rou2/2+2*w),1,size(Z,2)).*(S_hat + C)*rou2/2;
        rr = sqrt(kesi2+sum(Z.^2,2));
        %             w = lam^(-p)*rr.^(p-2);
        w = p*lam*alp*rr.^(p-2);
        w(w<max(w)*1e-6) = max(w)*1e-6;
        if  norm((Z - Z_inter_old),'fro')/norm(Z,'fro') < 1e-6 || sum(w) == 0
            break;
        end
        Z_inter_old = Z;
    end
%----------------------------W update----------------------------%
    W = W + (VS_hat - U);
%----------------------------C update----------------------------%
    C = C + (S_hat - Z);
%-------------------------stop criterion-------------------------%
    primerror1 = norm(VS - U,'fro');
    dualerror1 = norm(rou1*V'*(U - U_old),'fro');
    U_old = U;

    primerror2 = norm(S - Z,'fro');
    dualerror2 = norm(rou2*(Z - Z_old),'fro');
    Z_old = Z;

    tol_prim1 = 1e-6*max(norm(U,'fro'),norm(VS,'fro'));
    tol_dual1 = 1e-6*rou1*norm(V'*W,'fro');
    tol_prim2 = 1e-6*max(norm(Z,'fro'),norm(S,'fro'));
    tol_dual2 = 1e-6*rou2*norm(C,'fro');

    Rprimerror1 = primerror1/max(norm(U,'fro'),norm(VS,'fro'));
    Rdualerror1 = dualerror1/norm(rou1*V'*W,'fro');
    Rprimerror2 = primerror2/max(norm(Z,'fro'),norm(S,'fro'));
    Rdualerror2 = dualerror2/norm(rou2*C,'fro');

    if ~QUIET && mod(iter_ADMM,10) == 0
        fprintf('ADMM : iteration: %g, Data-Fit : %g, PrimError: %g, DualError: %g\n', ier_ADMM, norm(B - L*S,'fro')/norm(B,'fro'), primerror1, dualerror1);
    end

    if (primerror1 < tol_prim1 && dualerror1 < tol_dual1) || (primerror2 < tol_prim2 && dualerror2 < tol_dual2)
        break;
    end
%---------------------------rou update---------------------------%
    if rou_update && mod(iter_ADMM,10) == 0
        ratio1 = -1;
        ratio2 = -1;
        if Rdualerror1~=0
            ratio1 = sqrt(Rprimerror1/Rdualerror1);
        end
        if Rdualerror2~=0
            ratio2 = sqrt(Rprimerror2/Rdualerror2);
        end

        tau_max = 2;
        if ratio1>=1 && ratio1<tau_max, tau1 = ratio1;
        elseif ratio1>1/tau_max && ratio1<1, tau1 = 1/ratio1;
        else tau1 = tau_max;
        end
        if ratio2>=1 && ratio2<tau_max, tau2 = ratio2;
        elseif ratio2>1/tau_max && ratio2<1, tau2 = 1/ratio2;
        else tau2 = tau_max;
        end

        if Rprimerror1 > 10 * Rdualerror1
            rou1 = tau1*rou1; W = W./tau1;
        elseif Rdualerror1 > 10 * Rprimerror1
            rou1 = rou1/tau1; W = W.*tau1;
        end
        if Rprimerror2 > 10 * Rdualerror2
            rou2 = tau2*rou2; C = C./tau2;
        elseif Rdualerror2 > 10 * Rprimerror2
            rou2 = rou2/tau2; C = C.*tau2;
        end
        if ~QUIET
            fprintf('rou = %g, Rprimerror = %g, Rdualerror = %g\n',rou, Rprimerror1, Rdualerror1);
        end
        if rou1 ~= rou1_old || rou2 ~= rou2_old
            Lambda_MAX = eigs(V'*V+rou2/rou1*speye(Nsource),1,'lm');
            mu = 0.9/(rou1*Lambda_MAX);
            Precision = mu*speye(Nsource) - mu^2*L'/(eye(Nsensor) + mu*LLt)*L;
        end
        rou1_old = rou1;
        rou2_old = rou2;
    end
    if ~mod(iter_ADMM,100)
        SS{iter_ADMM/100} = S;
    end
    if ~mod(iter_ADMM,150)
%         kesi = std(VS,0,2).^2;
%         kesi2 = std(S,0,2).^2;
        kesi = max(VS, [], 2)*1e-3;
        kesi2 = max(S, [], 2)*1e-3;
    end
end
toc
S = SS;

function V = VariationEdge(VertConn)
Nsource = size(VertConn,1);
Nedge = numel(find(VertConn(:)~=0))/2;
V = sparse(Nedge,Nsource);
edge = 0;
for i = 1:Nsource
    idx = find(VertConn(i,:)~=0);
    idx = idx(idx<i);
    for j = 1:numel(idx)
        V(j+edge,i) = 1;
        V(j+edge,idx(j)) = -1;
    end
    edge = edge + numel(idx);
end

    

function Z = prox(Y,p,lam,rou)
k = 3;
temp = lam/rou;
thr = (2*temp*(1-p))^(1/(2-p)) + temp*p*(2*temp*(1-p))^((p-1)/(2-p));
Z = abs(Y);
Z(Z<=thr) = 0;
ind = find(Z>0);
for i = 1:k
    Z(ind) = abs(Y(ind)) - temp*p*Z(ind).^(p-1);
end
Z = sign(Y).*Z;
