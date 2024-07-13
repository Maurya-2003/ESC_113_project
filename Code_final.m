%ESC113 Project Group-2
%Calculation of different products in a hypothetical process of NO2 
%Let's define the initial conditions first...
%We initially consider an environment with 1M NO and 1M O2

%Though the step of the reaction (Formation of (NO)2) is a very fast
%process we assume it to be of a significantly observable pace in this
%hypothetical system to ease our analysis and to achieve a more efficient
%solution.

k1=0.1;

k2=0.001;

k3=0.001;

%Method-1 (Explicit Euler Method)
%Lets define The Initial Concentrations...

Conc_new=[1;0;1;0];%[NO,N2O2,O2,NO2]

%We assume the reaction duration to be 4 seconds and the step size for our
%method to be 0.1s.
h=0.1;

T=40;

N=T/h+1;

Conc_new=zeros(4,N);

Conc_new(:,1)=[1;0;1;0];

Conc_new
for i=2:N

    %As per the rate law we had obtained 4 differential equations for 4
    %unknowns which are mentioned below;

    dca=-k1*Conc_new(1,i-1)^2+k2*Conc_new(2,i-1);
    dcb=k1*Conc_new(1,i-1)^2-k2*Conc_new(2,i-1)-k3*Conc_new(3,i-1)*Conc_new(2,i-1);
    dcc=-k3*Conc_new(2,i-1);
    dcd=k3*Conc_new(2,i-1)*Conc_new(3,i-1);

    %Finally After Calculating the derivatives iterate over to the next
    %concentration of products in time.

    Conc_new(1,i)=Conc_new(1,i-1)+h*dca;
    Conc_new(2,i)=Conc_new(2,i-1)+h*dcb;
    Conc_new(3,i)=Conc_new(3,i-1)+h*dcc;
    Conc_new(4,i)=Conc_new(4,i-1)+h*dcd;
    %Conc_old=Conc_new;
end

Concexp=Conc_new(:,10);
Conc_new

%Plotting the data for Explicit Euler's Method...

plot(0:h:T,Conc_new(1,:),"--",0:h:T,Conc_new(2,:),"--",0:h:T,Conc_new(3,:),"--",0:h:T,Conc_new(4,:),"--");

title(["Data Analysis by Euler's Explicit Method"])
xlabel('time');
ylabel('concentration');

legend("[NO]","[N2O2]","[O2]","[NO2]")

figure

% Method-2 (Predictor-Corrector Method)

%The method is just the primitive version of RK2 (Runge Kutta 2) hence we
%shall see the computation of only one of them...

k1=0.1;k2=0.001;k3=0.001;

%Lets define The Initial Concentrations...

Conc_old=[1;0;1;0];%[NO,N2O2,O2,NO2]

%We assume the reaction duration to be 40 seconds and the step size for our
%method to be 0.1s.

h=0.1;  %Step Size

T=40;   %Time duration for the reaction...

N=T/h+1;  %Number of steps in the reaction...

Conc_new=zeros(4,N);  %Array for storing the concentration of the 4 components of the reaction...

Conc_new(:,1)=[1;0;1;0];  %Initializing the initial concentration for the reaction...

Conc_exp=zeros(4,N);  %Defining the concentrations obtained using the explicit methods...

Conc_new  
for i=2:N

    %As per the rate law we had obtained 4 differential equations for 4
    %unknowns which are mentioned below;

    dca=-k1*Conc_new(1,i-1)^2+k2*Conc_new(2,i-1);
    dcb=k1*Conc_new(1,i-1)^2-k2*Conc_new(2,i-1)-k3*Conc_new(3,i-1)*Conc_new(2,i-1);
    dcc=-k3*Conc_new(2,i-1);
    dcd=k3*Conc_new(2,i-1)*Conc_new(3,i-1);

    %Finally After Calculating the derivatives iterate over to the next
    %concentration of products in time (For getting the explicit version).

    Conc_new(1,i)=Conc_new(1,i-1)+h*dca;
    Conc_new(2,i)=Conc_new(2,i-1)+h*dcb;
    Conc_new(3,i)=Conc_new(3,i-1)+h*dcc;
    Conc_new(4,i)=Conc_new(4,i-1)+h*dcd;

    %After getting the explicit version we shall use the derivative
    %equations again to get the experession for Predictor-corrector(RK2)

    dca=dca-k1*Conc_new(1,i)^2+k2*Conc_new(2,i);
    dcb=dcb+k1*Conc_new(1,i)^2-k2*Conc_new(2,i)-k3*Conc_new(3,i)*Conc_new(2,i);
    dcc=dcc-k3*Conc_new(2,i);
    dcd=dcd+k3*Conc_new(2,i)*Conc_new(3,i);

    %Finally the actual value ofter the computation is assigned to the
    %array of Concentration of different products over time

    Conc_new(1,i)=Conc_new(1,i-1)+h*dca/2;
    Conc_new(2,i)=Conc_new(2,i-1)+h*dcb/2;
    Conc_new(3,i)=Conc_new(3,i-1)+h*dcc/2;
    Conc_new(4,i)=Conc_new(4,i-1)+h*dcd/2;

    %Conc_old=Conc_new;
end

Concpc=Conc_new(:,10); %Concentration of components at 10th second by
                       %Predictor Corrector method

Conc_new

%Plotting the data for Predictor-Corrector method...

plot(0:h:T,Conc_new(1,:),"--",0:h:T,Conc_new(2,:),"--",0:h:T,Conc_new(3,:),"--",0:h:T,Conc_new(4,:),"--");

title(["Data Analysis by Predictor-Corrector Method"])

legend("[NO]","[N2O2]","[O2]","[NO2]")
xlabel('time');
ylabel('concentration');

figure


%Method-3 (RK 4 Method)

%This method is a more accurate version of the Runge Kutta solving
%methods and it produces a global truncation error of O(h^4)

Conc_new=[1;0;1;0];%[NO,N2O2,O2,NO2]

%We assume the reaction duration to be 40 seconds and the step size for our
%method to be 0.1s.

h=0.1;   %Step Size

T=40; %Time duration for the reaction

N=T/h+1; % Number of steps;

Conc_new=zeros(4,N); %Array for storing the concentration of the 4 different components of the reaction

Conc_new(:,1)=[1;0;1;0]; %Defining the initial Concentration;

y = Conc_new(:,1);  %Initializing the initial concentration...

for i=2:N

      K1 = getd(y);         %Using the recuresive process we find the values of K1

      K2 = getd(y+(h/2)*K1); %K2

      K3 = getd(y+(h/2)*K2); %K3

      K4 = getd(y+h*K3);     %K4 and using these values getting the final value of the next iterative step...

      y = y + h/6*(K1+2*K2+2*K3+K4); %Storing the final value in the variable...

    Conc_new(:,i) = y;
end

Concrk=Conc_new(:,10);


Conc_new

Errexp=computenorm(Concrk,Concexp) %Calculating error in Explicit Method treating RK4 as threshold
Errpc=computenorm(Concrk,Concpc)  %Calculating error in Predictor-Corrector treating RK4 as threshold
disp("Predictor-Corrector is a more accurate method since Errexp>Errpc")
% plot(0:h:4,Conc_new);

%Plotting the data for Runge Kutta -4 method

plot(0:h:T,Conc_new(1,:),"--",0:h:T,Conc_new(2,:),"--",0:h:T,Conc_new(3,:),"--",0:h:T,Conc_new(4,:),"--");
title(["Data Analysis by Fourth Order Runge Kutta"])

legend("[NO]","[N2O2]","[O2]","[NO2]")
xlabel('time');
ylabel('concentration');

function D = getd(y)
 k1=0.1;k2=0.001;k3=0.001;
 D = zeros(4,1);
 D(1)=-k1*y(1)^2+k2*y(2);
 D(2)=k1*y(1)^2-k2*y(2)-k3*y(3)*y(2);
 D(3)=-k3*y(2);
 D(4)=k3*y(2)*y(3);
end

% function R=getR(x,y)
%     h=0.1;
%     k1=0.1;k2=0.001;k3=0.001;
%     R(1,1)=x(1)-y(1)-h*(-k1*x(1)^2+k2*x(2));
%     R(2,1)=x(2)-y(2)-h*(k1*x(1)^2-k2*x(2)-k3*x(3)*x(2));
%     R(3,1)=x(3)-y(3)-h*(-k3*x(3)*x(2));
%     R(4,1)=x(4)-y(4)-h*(k3*x(3)*x(2));
% end
% 
% function J=getJ(x)
%     h=0.1;
%     k1=0.1;k2=0.001;k3=0.001;
%     J=zeros(4,4);
%     J(1,1)=1+2*h*k1*x(1);
%     J(1,2)=-h*k2*x(2);
%     J(2,2)=1+h*k2+h*k3*x(3);
%     J(2,3)=h*k3*x(2);
%     J(3,2)=h*k3*x(3);
%     J(3,3)=h*k3*x(2);
%     J(4,2)=-k3*h*x(3);
%     J(4,3)=-h*k3*x(2);
%     J(4,4)=1;
% end
function val=computenorm(xold,xnew)
    val=0;
    for i=1:size(xold)
        val=val+(xnew(i)-xold(i))^2;
    end
    val=sqrt(val);
end