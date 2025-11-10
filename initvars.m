%This script provides the initial conditions used in the LLF method. The
%purpose of each initalized variable is commented beside it. To adjust for
%a different problem, simply edit parameters and change the left and right
%inital states.

% Figures are initalized here
figure(1), clf
% figure(2), clf

% Grid Variables
dx = 1; 
lx0 = 4; %Inital vector size
CFL = 0.5; %Numerical stability condition, controls dt/dx
renorm = 100; 
start = 0; 
totalTime = 1000; % total steps 

% initial line width for graphing
lU = 0.25;

% Values of Constants
pbar = 5;
aexp = -1.5;
a_s = 0;


% this is for the graph:
t_graph = 1;


% Initial Conditions given here
%uL corresponds to the left state of a characterstic, while uR corresponds
%to the right states
pL = 3;
uL = -3;
pR = 4;
uR = -4;  

% Determines case for graph titles
% Cases are determined by which of pL or pR is bigger, as well as the value
% of the constant aexp. There are 24 in total.
caseNum = 0;

if aexp < -1
    if uL < 0
        if pL < pbar
            caseNum = 1; % case 1
        elseif pL == pbar
            caseNum = 2; % case 2
        else
            caseNum = 3; % case 3
        end % pL > pbar
    else % uL > 0
        if pL < pbar
            caseNum = 4; % case 4
        elseif pL == pbar
            caseNum = 5; % case 5
        else % pL > pbar
            caseNum = 6; % case 6
        end
    end
elseif aexp == -1
    if uL < 0
        if pL < pbar
            caseNum = 7; % case 7
        elseif pL == pbar
            caseNum = 8; % case 8
        else % pL > pbar
            caseNum = 9; % case 9
        end
    else % uL > 0
        if pL < pbar
            caseNum = 10; % case 10
        elseif pL == pbar
            caseNum = 11; % case 11
        else % pL > pbar
            caseNum = 12; % case 12
        end
    end
elseif aexp > -1 && aexp < 0
    if uL < 0
        if pL < pbar
            caseNum = 13; % case 13
        elseif pL == pbar
            caseNum = 14; % case 14
        else % pL > pbar
            caseNum = 15; % case 15
        end
    else % uL > 0
        if pL < pbar
            caseNum = 16; % case 16
        elseif pL == pbar
            caseNum = 17; % case 17
        else % pL > pbar
            caseNum = 18; % case 18
        end
    end
elseif aexp > 0
    if uL < 0
        if pL < pbar
            caseNum = 19; % case 19
        elseif pL == pbar
            caseNum = 20; % case 20
        else % pL > pbar
            caseNum = 21; % case 21
        end
    else % uL > 0
        if pL < pbar
            caseNum = 22; % case 22
        elseif pL == pbar
            caseNum = 23; % case 23
        else % pL > pbar
            caseNum = 24; % case 24
        end
    end
else
    error('error lol');
end
