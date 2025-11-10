fig = figure(1); %, clf
set(gcf, 'Units', 'pixels', 'Position', [100 100 600 600]);

s1 = subplot(2,1,1);
set(s1, 'Units', 'pixels', 'Position', [60 342 510 200])  % Top subplot
plot(X,p,'k','LineWidth',lU);
title(['Case ',num2str(caseNum),': \rho bar = ',num2str(pbar),', a = ',num2str(aexp), ...
 ', Steps = ', num2str(iters), ', t = ', num2str(round(t, 2))],'FontSize',14)
	xlabel('\itx/t', 'FontSize',16), ylabel('\itp','FontSize',16)


    
s2 = subplot(2,1,2);
set(s2, 'Units', 'pixels', 'Position', [60 70 510 200])
plot(X,u,'k','LineWidth',lU), ylabel('\itu','FontSize',16)
xlabel(['\itx/t; \rm Data \itL \rm= (',num2str(pL),',', num2str(uL), ...
	'), \itR \rm= (',num2str(pR),',',num2str(uR),')' ],'FontSize',8)
    

lU = lU + 0.05; % increase line width after each iteration
