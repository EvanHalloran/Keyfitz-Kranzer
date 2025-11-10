% graphs p and u
%   lU: variable to control line width, increases after each iteration
%   fig: figure object with all p,u plots
%   fig2: figure object with latest p,u plot

fig = figure(1); %, clf
hold on
subplot(2,2,1), hold on
plot(X,p,'k','LineWidth',lU);

if a_s ==0
title(['                                     Case ',num2str(caseNum),': \rho bar = ',num2str(pbar),', a = ',num2str(aexp), ...
 ', Steps = ', num2str(iters), ', t = ', num2str(round(t, 2))],'FontSize',16)
else
title(['                                                    Case ',num2str(caseNum),': a(t) = ', num2str(a_s) ', \rho bar = ',num2str(pbar),', a = ',num2str(aexp), ...
 ', Steps = ', num2str(iters), ', t = ', num2str(round(t, 2))],'FontSize',12)
end
	xlabel('\itx/t', 'FontSize',16), ylabel('\itp','FontSize',16)

    
subplot(2,2,3), hold on
plot(X,u,'k','LineWidth',lU), ylabel('\itu','FontSize',16)
xlabel(['\itx/t; \rm Data \itL \rm= (',num2str(pL),',', num2str(uL), ...
	'), \itR \rm= (',num2str(pR),',',num2str(uR),')' ],'FontSize',8)

hold off
subplot(2,2,2)
plot(X,p,'k','LineWidth',lU)
	xlabel('\itx/t', 'FontSize',16), ylabel('\itp','FontSize',16)


subplot(2,2,4), hold off
plot(X,u,'k','LineWidth',lU), ylabel('\itu','FontSize',16)
xlabel(['\itx/t; \rm Data \itL \rm= (',num2str(pL),',', num2str(uL), ...
	'), \itR \rm= (',num2str(pR),',',num2str(uR),')' ],'FontSize',8)


lU = lU + 0.05; % increase line width after each iterationgraphs.m
