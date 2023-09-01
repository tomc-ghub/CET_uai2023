%% load data
% load ./data/Test_UAI_CPAG_A_pC100.mat    % matfile with results d=3.0, 7 sizes, 100 graphs per size
% ACC_A = ACC;
% load ./data/Test_UAI_CPAG_B_pC100.mat     % idem d=5.0
% ACC_B = ACC;

%% get avg/max per N per density Alg per 

for j = 1:7
  i = (j-1)*100 + 1;
  % density = 3.0 (version A)
  % org
  avg_A1(j,:) = sum(ACC_A.R1(i:(i+99),:),1)/100; 
  max_A1(j,:) = max(ACC_A.R1(i:(i+99),:),[],1); 
  % new
  avg_A2(j,:) = sum(ACC_A.R2(i:(i+99),:),1)/100; 
  max_A2(j,:) = max(ACC_A.R2(i:(i+99),:),[],1); 
  
  % density = 5.0 (version B)
  avg_B1(j,:) = sum(ACC_B.R1(i:(i+99),:),1)/100; 
  max_B1(j,:) = max(ACC_B.R1(i:(i+99),:),[],1); 
  %
  avg_B2(j,:) = sum(ACC_B.R2(i:(i+99),:),1)/100; 
  max_B2(j,:) = max(ACC_B.R2(i:(i+99),:),[],1); 
end;

%% add average time per stage
% org: columns (6,9,12,15,18,21 / 24)*100 
Time_A1 = (avg_A1(:,[6:3:21])./( sum(avg_A1(:,[6:3:21]),2)))*99.999;
Time_B1 = (avg_B1(:,[6:3:21])./( sum(avg_B1(:,[6:3:21]),2)))*99.999;

% new: columns [6,9,12,18,21 / 33
Time_A2 = (avg_A2(:,[6,9,12,15,18,21])./( sum(avg_A2(:,[6,9,12,15,18,21]),2)))*99.999;
Time_B2 = (avg_B2(:,[6,9,12,15,18,21])./( sum(avg_B2(:,[6,9,12,15,18,21]),2)))*99.999;


%% graphs
figure(5); clf; loglog(avg_A1(:,1),avg_A1(:,24),'-or','LineWidth',2, ...
  'DisplayName','avg. CPAG-from-Graph (d=3.0)');
hold on; loglog(max_A1(:,1),max_A1(:,24),'--r','LineWidth',1, ...
  'DisplayName','idem max.');
hold on; loglog(avg_B1(:,1),avg_B1(:,24),'-om','LineWidth',2, ...
  'DisplayName','avg. CPAG-from-Graph (d=5.0)');
hold on; loglog(max_B1(:,1),max_B1(:,24),'--m','LineWidth',1, ...
  'DisplayName','idem max.');

hold on; loglog(avg_A2(:,1),avg_A2(:,33),'-ob','LineWidth',2, ...
  'DisplayName','avg. new G-to-CPAG (d=3.0)');
hold on; loglog(max_A2(:,1),max_A2(:,33),'--b','LineWidth',1, ...
  'DisplayName','idem max.');
hold on; loglog(avg_B2(:,1),avg_B2(:,33),'-og','LineWidth',2, ...
  'DisplayName','avg. new G-to-CPAG (d=5.0)');
hold on; loglog(max_B2(:,1),max_B2(:,33),'--g','LineWidth',1, ...
  'DisplayName','idem max.');
%
title('Running time CPAG reconstruction','FontSize',16);
xlabel('Number of nodes N','FontSize',14);
ylabel('Time T (sec)','FontSize',14);
ax = gca;
ax.XMinorTick = 'on';
ax.XAxis.TickValues = [10,20,40,60,80,100,200];
ax.XAxis.Exponent = 1;
legend('Location','northwest','FontSize',12);



%% grouped version
X = categorical({'10','20','40','60','80','100','200'});
X = reordercats(X,{'10','20','40','60','80','100','200'});
figure(1); clf;
subplot(2,2,1) 
bar(X,Time_A1,'stacked');
title('CPAG-from-Graph (d=3.0)') 
xlabel 'Number of nodes N'; 
ylabel 'Time per stage (%)'; 
subplot(2,2,2) 
bar(X,Time_B1,'stacked');
title('CPAG-from-Graph (d=5.0)') 
xlabel 'Number of nodes N'; 
ylabel 'Time per stage (%)'; 
subplot(2,2,3) 
bar(X,Time_A2,'stacked');
title('New Graph-to-CPAG (d=3.0)') 
xlabel 'Number of nodes N'; 
ylabel 'Time per stage (%)'; 
subplot(2,2,4) 
bar(X,Time_B2,'stacked');
title('New Graph-to-CPAG (d=5.0)') 
xlabel 'Number of nodes N'; 
ylabel 'Time per stage (%)'; 


