function [ anchor, ind2, score ] = gen_anchor_std_en(X, m)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[n,d] = size(X);
X_std = std(X,[],2);
% vm = min(X,[],2);
% Xm = ones(n,1)*vm;
% X = X-vm;
% for i=1:d
%     maxd=max(X(:,i));
%     X(:,i)=X(:,i)./maxd;
% end
% score = sum(X, 2);
score = X_std.^2;
score(:,1) = score/max(score);
[~,ind(1)] = max(score);
for i=2:m
    for j=1:n
        A_1 = score(ind(i-1),i-1);
        A_2 = score(j,i-1);
        Co(j,:)=(1 + norm(A_1-A_2,2)^(0.5))^(-1);
    end
    pho = Co/max(Co);
    score(:,i) = score(:,i-1).*(ones(n,1)-pho);
    score(:,i) = score(:,i)/max(score(:,i));
    [~,ind(i)] = max(score(:,i));
end
ind2 = sort(ind,'ascend');
anchor = X(ind2,:);

% for i=1:4
% % idd=find(score(:,i)>0.98);
% % figure;
% % plot(X(:,1),X(:,2),'.b', 'MarkerSize', 10); hold on;
% % plot(X(idd,1),X(idd,2),'.r', 'MarkerSize', 10); hold on;
% figure; plot(score(:,i),'-o');
% % axis equal;
% set(gcf,'Position',[400,100,700,600],'color','w');
% set(gca,'fontsize',16);
% set(gca,'linewidth',0.8);
% % saveas(gcf,strcat('C:\Users\opt\Desktop\mydisk\ongoing\Fast multiview CLR\Latex_revised2\response\s',num2str(i),'.pdf'));
% end


% id=zeros(1,m);
% for i=1:m
%     if ind(i)<=200
%         id(i)=1;
%     elseif ind(i)<=400
%         id(i)=2;
%     else
%         id(i)=3;
%     end
% end
end