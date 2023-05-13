%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Project #1
%       Ethan Swan-Begnaud
%       Jason Rodriguez
%       Nicolas Kudr
%       Ian Coble
%       Marlon Coates
%   ME 2543--Simulations Methods
%   Spring 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QUESTION 1

clear all; close all; clc

Qtot = 1.671;                       % Cubic Feet per Second
L = [10 10 10];                     % Feet
D = [1/6 1/6 1/6];                  % Feet
e = [0.00085 0.00085 0.00085];      % Feet
p = 1.94;                           % Slugs per Cubic Foot
m = 2.05*(10^-5);                   % Pound-Force-Seconds per Square Foot
A = pi*(D./2).^2;                    % ft^2

Guess = [0.5 0.5 0.5 70];
% Our initial guess for Q1, Q2, Q3, and DeltaP, respectively.

oldGuess = zeros(1,4);
% Initialize f

while abs(Guess(1)-oldGuess(1)) > 0.001

    Re = (p.*Guess(1:end-1).*D)./(A.*m);

    oldGuess = Guess;

    if Re(1) <= 2300
        f =  64./Re;
    else
        f = 0.25.*(log10(((e./D)/3.7)+5.74./(Re.^(0.9)))).^-2;
    end 
    func = @(Q)equations(Q,Qtot,p,f,L,D,A);
    Guess = fsolve(func,Guess);

end

clc
sprintf(['The three flow rates are %.2f, %.2f, and %.2f GPM, respectively. ' ...
    'The total change in pressure through the network is %.2f PSI.'], ...
    Guess(1,1)*(448.831),Guess(1,2)*(448.831),Guess(1,3)*(448.831),Guess(1,4)/144)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QUESTION 2

clear all; close all; clc

Qtot = 1.671;                       % Cubic Feet per Second
L = [20 10 30];                     % Feet
D = [2/12 2.5/12 1.5/12];                  % Feet
e = [0.00085 0.00085 0.00085];      % Feet
p = 1.94;                           % Slugs per Cubic Foot
m = 2.05*(10^-5);                   % Pound-Force-Seconds per Square Foot
A = pi*(D./2).^2;                    % ft^2

Guess = [0.5 0.5 0.5 70];
% Our initial guess for Q1, Q2, Q3, and DeltaP, respectively.

oldGuess = zeros(1,4);
% Initialize f

while abs(Guess(1)-oldGuess(1)) > 0.001

    Re = (p.*Guess(1:end-1).*D)./(A.*m);

    oldGuess = Guess;
    
    if Re(1) <= 2300
        f =  64./Re;
    else
        f = 0.25.*(log10(((e./D)/3.7)+5.74./(Re.^(0.9)))).^-2;
    end 
    func = @(Q)equations(Q,Qtot,p,f,L,D,A);
    Guess = fsolve(func,Guess);

end

clc
sprintf(['The three flow rates are %.2f, %.2f, and %.2f GPM, respectively. ' ...
    'The total change in pressure through the network is %.2f PSI.'], ...
    Guess(1,1)*(448.831),Guess(1,2)*(448.831),Guess(1,3)*(448.831),Guess(1,4)/144)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 3

clear all; close all; clc

L = [20 10 30];                     % Feet
D = [2/12 2.5/12 1.5/12];           % Feet
e = [0.00085 0.00085 0.00085];      % Feet
p = 1.94;                           % Slugs per Cubic Foot
m = 2.05*(10^-5);                   % Pound-Force-Seconds per Square Foot
A = pi*(D./2).^2;                   % ft^2

cleanFoP = zeros(2,26);

for i = 1:26
    Qtot = 100/448.831 + (i-1)*((1400/25)/448.831);   % Cubic Feet per Second

Guess = [0.5 0.5 0.5 70];
% Our initial guess for Q1, Q2, Q3, and DeltaP, respectively.

oldGuess = zeros(1,4);
% Initialize f

while abs(Guess(1)-oldGuess(1)) > 0.001

    Re = (p.*Guess(1:end-1).*D)./(A.*m);

    oldGuess = Guess;
    
    if Re(1) <= 2300
        f =  64./Re;
    else
        f = 0.25.*(log10(((e./D)/3.7)+5.74./(Re.^(0.9)))).^-2;
    end 
    func = @(Q)equations(Q,Qtot,p,f,L,D,A);
    Guess = fsolve(func,Guess);

end   
    cleanFoP(1,i) = Guess(4)/144;
    cleanFoP(2,i) = Qtot*(448.831);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 4.a

fouledFoP_A = zeros(2,26);

for i = 1:26
    Qtot = 100/448.831 + (i-1)*((1400/25)/448.831);   % Cubic Feet per Second

Guess = [0.5 0.5 0.5 70];
% Our initial guess for Q1, Q2, Q3, and DeltaP, respectively.

oldGuess = zeros(1,4);
% Initialize f

while abs(Guess(1)-oldGuess(1)) > 0.001

    Re = (p.*Guess(1:end-1).*D)./(A.*m);

    oldGuess = Guess;
    
    if Re(1) <= 2300
        f =  64./Re;
    else
        f = 0.25.*(log10((((e.*1.25)/D)/3.7)+5.74./(Re.^(0.9)))).^-2;
    end 
    func = @(Q)equations(Q,Qtot,p,f,L,D,A);
    Guess = fsolve(func,Guess);

end   
    fouledFoP_A(1,i) = Guess(4)/144;
    fouledFoP_A(2,i) = Qtot*(448.831);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 4.b

fouledFoP_B = zeros(2,26);

for i = 1:26
    Qtot = 100/448.831 + (i-1)*((1400/25)/448.831);   % Cubic Feet per Second

Guess = [0.5 0.5 0.5 70];
% Our initial guess for Q1, Q2, Q3, and DeltaP, respectively.

oldGuess = zeros(1,4);
% Initialize f

while abs(Guess(1)-oldGuess(1)) > 0.001

    Re = (p.*Guess(1:end-1).*D)./(A.*m);

    oldGuess = Guess;
    
    if Re(1) <= 2300
        f =  64./Re;
    else
        f = 0.25.*(log10((((e.*1.35)/D)/3.7)+5.74./(Re.^(0.9)))).^-2;
    end 
    func = @(Q)equations(Q,Qtot,p,f,L,D,A);
    Guess = fsolve(func,Guess);

end   
    fouledFoP_B(1,i) = Guess(4)/144;
    fouledFoP_B(2,i) = Qtot*(448.831);
end

clc

hold on
grid on
CleanGraph = plot(cleanFoP(1,1:end),cleanFoP(2,1:end),"Black", ...
    fouledFoP_A(1,1:end),fouledFoP_A(2,1:end),"R", ...
    fouledFoP_B(1,1:end),fouledFoP_B(2,1:end),"B");
xlabel('Pressure(psi)')
ylabel('Total Flow Rate (gpm)')
title('Total Flow Rate vs Pressure')
legend("Clean Pipes","Fouled A","Fouled B")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions

function F = equations(Q,Qtot,p,f,L,D,A)
    F(1) = Q(1) + Q(2) + Q(3) - Qtot;
    F(2) = p*f(1)*(L(1)/D(1))*((Q(1)^2)/(2*(A(1)^2))) - Q(4);
    F(3) = p*f(2)*(L(2)/D(2))*((Q(2)^2)/(2*(A(2)^2))) - Q(4);
    F(4) = p*f(3)*(L(3)/D(3))*((Q(3)^2)/(2*(A(3)^2))) - Q(4);
end














