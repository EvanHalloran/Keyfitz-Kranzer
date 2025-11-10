% Autogen Script
% NOTE: You only need to run this script. This calls everything else.
% Iterates through laxfried desired amount of times
%   numInner: controls precisely how many times laxfried is run for each
%             point
%   numOuter: controls how many points are run, only need this variable if 
%             using the random equations, if not change the outer for loop
%             to be 1:length(randomVR) 
% 

clear all;
clc; 
close all;

% vars to control loop
numInner = 20; % numInner controls how many times laxfried is run
numOuter = 3; % numOuter controls how many points/ figures are generated

% NOTE: the following takes random points but you can comment this out and preset
% randomVR and randomVR as vectors with your desired coordinates instead of
% using the random point generation
%
% for N random numbers between [a, b], use r = a + (b-a).*rand(1,N) from
% MATLAB documentation 

randomPR = [1];

%4.5 4.5 4.5 7.5 3.01 2.2 1.7 3.4 4
%5.*rand(1,numOuter); randomVR = floor(randomVR .* 10) ./ 10; % generates numOuter numbers between 0 and 5
   
    % randomUR = [2];

%2.5 1.7 1.9 1.9 1.9 0.2 -1.2 1.05 1.39
%-5 + 10.*rand(1,numOuter); randomWR = floor(randomWR .* 10) ./ 10; % generates numOuter numbers between -5 and 5

% generate numOuter figures
for j = 1:length(randomPR)
    
    initvars; % initialize variables

    % create filename base
    nameFile = ['Case' num2str(caseNum) 'a(t)' num2str(a_s) 'L(' num2str(pL) ',' num2str(uL) ')R(' num2str(pR) ',' num2str(uR) ')'];


%    pR = randomPR(j); uR = randomUR(j); % assign right states each iteration

    % setup taken from initvars since it reruns every time 
    qL = pL*uL;
    qR = pR*uR;

    UL=[pL qL];
    UR=[pR qR];
    
    p = pL * ones(1,lx0); 
    q = qL * ones(1,lx0);
    u = uL * ones(1,lx0);
%meow
    %Will need to add right endpoints and finish constructing inital variables
    for i = 3:lx0
        p(i) = p(i) + (pR - pL);
        q(i) = q(i) + (qR - qL);
        u(i) = u(i) + (uR - uL);
    end
    %
    % end of setup

    % loop through laxfried numInner times
    for i = 1:numInner
        laxfried; % call laxfried
    end

    % save files as png
    saveas(fig,[nameFile '.png']);

end

% create figure and get curves for shocks. t=1 
figure;
plotX = linspace(0,10,100);

u1 = uL .* ones(size(plotX));
u2 = (((pL^aexp) - (plotX.^aexp))./((plotX.^aexp) - (pbar^aexp))) .* (uL + (a_s.*t_graph)) + uL;

% plot up curves and points
hold on;

% u2 becomes a vertical line when pbar = pL
plot(plotX,u1,'-r');
if pbar == pL
    xline(pbar, '-b');
else
    plot(plotX,u2,'-b');
end


% plots the right state point that we r looking at
plot(pR,uR,'.r')

% asymptote (?) A 
plot(plotX,-(a_s.*t_graph).*ones(length(plotX)), ':k')

% asymptote B (it's purple lol) 
B = -((pL/pbar)^aexp)*uL + (uL - a_s.*t_graph);
plot(plotX,B.*ones(length(plotX)), ':m')

% pbar asymptote 
xline(pbar, ':k')

%y-axis and x-axis are solid black lines
 xline(0);
 yline(0);

 % this plots the path it takes
 plot(p,u,"magenta")

ylim([-10,10])
