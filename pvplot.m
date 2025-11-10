% pvplot Script
% this is just a script to print out the locus graphs. automatically calls 
% initvars to get initial variables :)
% however, this script only allows for ONE POINT at a time, not a vector

initvars;

% create figure and get curves for shocks. t=1 
fig4 = figure(4);
plotX = linspace(0,10,100);

u1 = uL .* ones(size(plotX));
u2 = (((pL^aexp) - (plotX.^aexp))./((plotX.^aexp) - (pbar^aexp))) .* (uL + (a_s.*1)) + uL;

% plot up curves and all points (in red)
hold on;

% u2 becomes a vertical line when pbar = pL
plot(plotX,u1,'-r');
if pbar == pL
    xline(pbar, '-b');
else
    plot(plotX,u2,'-b');
end

plot(pR,uR,'.r');

% asymptote (?) A 
plot(plotX,-(a_s.*1).*ones(length(plotX)), ':k')

% asymptote B
B = -((pL/pbar)^aexp)*uL + (uL - a_s.*1);
plot(plotX,B.*ones(length(plotX)), ':m')

% pbar asymptote 
xline(pbar, ':k')

%y-axis and x-axis are solid black lines
 xline(0);
 yline(0);

ylim([-10,10])
