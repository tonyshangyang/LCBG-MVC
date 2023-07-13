clear all;
dataname='Handwritten4';
addpath([pwd, '/funs']);
addpath([pwd, '/datasets']);

load(strcat(dataname,'.mat'));

nC = length(unique(Y));
p = [0.9];
lambda1 = 1/(sqrt(size(X{1},1)*length(X)));
lambda2 = 1/(size(X{1},1));
lambda3 = [0.02];
rho1 = [0.1];
rho2 = [0.1];
anchor_rate = 0.5;
beta = ones(length(X),1);
IterMax = 50;


%% 数据初始化A
% nV = length(X);
% disp('------Data preprocessing------');
% tic
% for v=1:nV
%     a = max(X{v}(:));
%     X{v} = double(X{v}./a);
% end
% toc

%% 数据初始化B
nV = length(X);
disp('------Data preprocessing------');
tic
for v=1:nV
    XX = X{v};
    for n=1:size(XX,1)
        XX(n,:) = XX(n,:)./norm(XX(n,:),'fro');
    end
    X{v} = double(XX);
end
toc


filename=['result-' dataname '.txt'];
fid = fopen(filename,'a');
time_start = clock;
for num1 = 1:length(p)
    for num2 = 1:length(lambda1)
        for num3 = 1:length(lambda2)
            for num4 = 1:length(lambda3)
                for num5 = 1:length(rho1)
                    for num6 = 1:length(rho2)
                        for num7 = 1:length(anchor_rate)
                            [final_result] = My_main_test(X, Y, nC, IterMax, p(num1), beta, lambda1(num2), lambda2(num3), lambda3(num4), rho1(num5), rho2(num6),anchor_rate(num7));
                            for n_result = 1:length(final_result)
                                fprintf(fid, '%f ' ,final_result(n_result));
                            end
                            fprintf('p=%f lambda1=%f lambda2=%f lambda3=%f rho1=%f rho2=%f anchor_rate=%f\n', p(num1), lambda1(num2), lambda2(num3), lambda3(num4), rho1(num5), rho2(num6),anchor_rate(num7));
                            fprintf(fid, 'p=%f lambda1=%f lambda2=%f lambda3=%f rho1=%f rho2=%f anchor_rate=%f\n', p(num1), lambda1(num2), lambda2(num3), lambda3(num4), rho1(num5), rho2(num6),anchor_rate(num7));
                        end
                    end
                end
            end
        end
    end
end
fclose(fid);


time_end = clock;
fprintf('总耗时:%f秒。\n',etime(time_end,time_start))
