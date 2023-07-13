function [final_result] = My_main_test(X, Y, nC, IterMax, p, beta, lambda1, lambda2, lambda3, rho1, rho2, anchor_rate)
% X为输入，Y为标签，nC为类别数，lambda2为Sv范数前参数，lambda3为trace()前方参数自适应更新
% beta张量的权重
nV = length(X);     %视图数
N = size(X{1},1);   %样本数
M = fix(N*anchor_rate);  %fix朝零方向舍入，锚点数
betav = repmat(1/nV, [1,nV]);  %初始化权重向量
%% ===initialization=== %%
for k = 1:nV
    % 需要重新初始化
    MM{k} = zeros(N, M);
    W{k} = zeros(N, M);
    T{k} = zeros(N, M);
    U{k} = zeros(N, M);
    V{k} = zeros(N, M);
    TUV{k} = zeros(N, M);
end
iter = 0;
Isconverg = 0;
F = zeros(N+M,nC);
Fn = zeros(N,nC);
Fm = zeros(M,nC);
sX = [N, M, nV];  %锚点图
max_rho1 = 10e10;
max_rho2 = 10e10;
pho_rho = 1.1;
final_result = zeros(1,7);

%%
tic;
% 初始化Sv
opt1. style = 1;
opt1. IterMax = IterMax;
opt1. toy = 0;
[S] = Construct_Graph(X,nC,anchor_rate, opt1,10);

while(Isconverg == 0)
    %% update F
    ZZ = 0;
    for v = 1: nV
        Dm{v} = diag(sum(S{v},1)+eps); %eps为一个很小的值，按列求和得到行向量，作为对角线元素
        Dn{v} = diag(sum(S{v},2)+eps); %按行求和得到列向量，作为对角线元素
        ZZ = ZZ + (1/betav(v))*(Dn{v}^-0.5)*S{v}*(Dm{v}^-0.5);  %E
    end
    [uu, ~, vv] = svd(ZZ);
    Fn = uu(:,1:nC)*(2^-0.5);
    Fm = vv(:,1:nC)*(2^-0.5);
    F=[Fn; Fm];
    
    %%  update J{v}
    for v =1:nV
        SW{v} = S{v} + W{v}./rho1;
    end
    SW_tensor = cat(3,SW{:,:});
    x = SW_tensor(:);  %多维张量 SW_tensor 转换为一个向量 x
    [myj, ~] = wshrinkObj_weight_lp(x, 100*lambda1*beta./rho1,sX, 0,3,p);
    J_tensor = reshape(myj, sX);
    for k=1:nV
        J{k} = J_tensor(:,:,k);
    end
    
    %%  update P{v}
    for v =1:nV
        [AA,BB,CC] = svd(S{v} + MM{v}./rho2,'econ');  %返回紧凑型奇异值分解结果
        a = diag(BB)-lambda2/rho2;
        a(a<0)=0;
        P{v} = AA*(diag(a))*CC';
    end
    
    %% update S
    for v =1:nV
        T{v} = (1/betav(v))*(Dm{v}^-0.5)*Fm*Fn'*(Dn{v}^-0.5);
        U{v} = J{v} - W{v}./rho1;
        V{v} = P{v} - MM{v}./rho2;
        TUV{v} = (rho1*U{v} + rho2*V{v} + 2*lambda3*T{v}')./(rho1 + rho2);
    end
    
    for v = 1:nV
        for i = 1:N
            TUV_v = TUV{v}(i,:);
            S{v}(i, :) = EProjSimplex_new(TUV_v, 1);
        end
    end
    
    %% update lambda2
    clear SS;
    SS = (1/betav(1))*S{1};
    for v = 2: nV
        SS = SS + (1/betav(v))*S{v};
    end
    Sum_betav = sum(1./betav);
    SS = SS./Sum_betav;
    
	%%聚类
    try
        [Flabel] = coclustering_bipartite_fast1(SS, nC, IterMax);   
    catch
        Isconverg  = 1;
    end
    result = ClusteringMeasure1(Y,Flabel);   %评价
    for n_result = 1:length(result)
        fprintf('%f ' ,result(n_result))
    end
    if (sum(result) - sum(final_result))>0
        final_result = result;
    end
    fprintf('\n')
    
    Final_S = SS;
    [~, ev1, ~] = svd(Final_S);
    ev = diag(ev1);   %提取对角线元素
    fn1 = sum(ev(1:nC));  %前 nC 个奇异值的和
    fn2 = sum(ev(1:nC+1));
    if fn1 < nC-0.0000001
        lambda3 = 2*lambda3;
    elseif fn2 > nC+1-0.0000001
        lambda3 = lambda3/2;
    else
       break
    end    
    
    %% update betav{v}
    sum_h = 0;
    for v = 1:nV
        Dm{v} = diag(sum(S{v},1)+eps);
        Dn{v} = diag(sum(S{v},2)+eps);
        Dv = blkdiag(Dn{v},Dm{v}); %分块对角矩阵，总度矩阵
        tmp1 = zeros(N+M);   %邻接矩阵
        tmp1(1:N,N+1:end) = S{v};
        tmp1(N+1:end,1:N) = S{v}';
        Ls = eye(M+N) - (Dv^-0.5)*tmp1*(Dv^-0.5);
        h(v) = trace(F'*Ls*F);
        sum_h = sum_h + h(v)^0.5;
    end
    for v = 1:nV
        betav(v) = h(v)^0.5/sum_h;
    end
    
    %% update W M
    for v = 1:nV
        W{v} = W{v} + rho1*(S{v} - J{v});
        MM{v} = MM{v} + rho2*(S{v} - P{v});
    end
    
    %% update rho
    rho1 = min(rho1*pho_rho, max_rho1);
    rho2 = min(rho2*pho_rho, max_rho2);
    
    %%
    if (iter > IterMax)
        Isconverg  = 1;
    end
    iter = iter + 1;
end
toc;

fprintf('Final_iter:%d\n',iter)
