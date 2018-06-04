% pozar_02_exercise_29.m
%
% Tranmission line with attenuation constant alpha=.5/lambda dB/m lambda:wave-length
% matched generator, load not matched.
% Results: Input Power, Power on Load, Power loss.
% Smith Chart plot showing how reflection coefficient no longer moves along a circle.
%
% MATLAB script author: John Bofarull. 
% Send comments to: jgb2014@live.co.uk
%
% v1.0 Jun3rd18:
% 

% here gamma= alpha+1j*beta, complex (with non null real part) propagation constant [POZAR] pg20.
% in other exercises gamma is used as reflection coefficient

clear all;close all;clc

str65='pozar_02_exercise_29_vars.mat';
rel1=load_results_if_available(str65)

switch rel1
    
    case 0
        return                                                          % results already stored in .mat file.
        
    case 1

Z0=50           % ohm
Zgen=50
Vgen=10       % Volt, phasor
ZL=100

syms lambda 
alpha_dB=.5/lambda  % [dB/lambda]
Np2dB=10*log10(exp(1)^2)
alpha_Np=alpha_dB/Np2dB

beta=2*pi/lambda  % [m^-1]
L=2.3*lambda        % length of transmission line [m]

gamma=alpha_Np+1j*beta
gamma_L=gamma*L
gamma_L_rad=double(gamma_L)  % =   0.132398642847158 +14.451326206513048i

% sometimes for gamma Np/m and degree are mixed, for this exercise
% the solutions manual shows gamma=.1325 +1j*108º
% 108 degree is the remainder of imag(gamma_degree)/360

14.4513*180/pi-floor(14.4513*180/pi/360)*360

% reflection coefficients on both sides of the transmission line
refl_Load=(ZL-Z0)/(ZL+Z0)
refl_gen=refl_Load*exp(-2*gamma_L_rad)

Zin=Z0*(ZL+Z0*tanh(gamma_L_rad))/(Z0+ZL*tanh(gamma_L_rad))

% Pin: TL input power
Pin=(abs(Vgen*Zin/(Zgen+Zin)))^2/(2*Z0)*(1-(abs(refl_gen))^2)*exp(2*real(gamma_L_rad))

% PLoad: power reaching load
PLoad=(abs(Vgen*Zin/(Zgen+Zin)))^2/(2*Z0)*(1-(abs(refl_Load))^2)

% Ploss: lost power
Ploss=Pin-PLoad

% Comment: 
% maximum available power from generator
Pmax_gen=.5*(abs(Vgen))^2/(4*real(Zgen))  % .25 W

% Alternatively, perturbation method pg82
% alpha_dB and alpha_Np are per metre.

alpha=double(alpha_Np*L)
% alpha =   0.132398642847158

Pin*exp(-2*alpha)


% Comment 2: discrepancy with solutions manual: V0plus~=Vgen/2*exp(-alpha*L)

V0plus=Vgen/2*exp(-real(gamma_L_rad))
absV0plus=abs(V0plus)

% from  2.89a in [POZAR] pg81

Vin=Vgen*Zin/(Zin+Zgen)*1/(exp(gamma_L_rad)+refl_gen*exp(-gamma_L_rad))
abs(Vin)


% This Pin is the same as the above used

Pin=(abs(Vgen*Zin/(Zgen+Zin)))^2/(2*Z0*(1-(abs(refl_gen))^2)*exp(2*real(gamma_L_rad)))


% And it's not the other apparently correct

Vgen*Zin/(Zin+Zgen)  % =  3.965319190462363 + 0.751739611040107i

abs(Vgen*Zin/(Zin+Zgen))  % =   4.035947066681602


Loss_dB=double(alpha_dB*L)

10*log10(Pin)-10*log10(PLoad)  % = 1.367657186248119 dB

10*log10(exp(2*alpha))

% 0.21 dB missing, the perturbation method is not exact

% Vgen*Zin/(Zgen+Zin)*(1+refl_gen) should be used, 2.89a in pg81

Zin=Z0*(ZL+Z0*tanh(gamma_L_rad))/(Z0+ZL*tanh(gamma_L_rad))
refl_gen=(Zin-Z0)/(Zin+Z0)
Vin=Vgen*Zin/(Zin+Zgen)*(1+refl_gen)

abs(Vin)   %  =  3.257773745011164

% Comment: Zin with lossless TL would be
Zin_lossless=Z0*(ZL+1j*Z0*tan(imag(gamma_L_rad)))/(Z0+1j*ZL*tan(imag(gamma_L_rad)))

% same as
Zin_lossless=Z0*(ZL+1j*Z0*tan(2*pi*2.3))/(Z0+1j*ZL*tan(2*pi*2.3))


%% Now with Smith Chart:

Z0=50; ZL=100;

sm1=smithchart; ax=gca; hold all; 
refl_ZL=(ZL-Z0)/(ZL+Z0); 

syms lambda
alpha_dB=.5/lambda  % [dB/lambda]
Np2dB=10*log10(exp(1)^2)
alpha_Np=alpha_dB/Np2dB
beta=2*pi/lambda
L=2.3*lambda

g=alpha_Np+1j*beta    % g=gamma=alpha+1j*beta
g_L=g*L
g_L=double(simplify(g_L))

% angle to run along TL: beta*L=pi rad means 360º around Smith Chart  
% so to match the imag(g_L) angle with Smith chart angle, double it.
N1=100
da=2*pi/N1  % angle differential for 2*pi around Smith Chart = lambda/2 

amount_full_turns=floor(14.4513*180/pi/360)

% rads left when no full turns left
14.4513-floor(14.4513*180/pi/360)

amount_da_rem=floor((14.4513-floor(14.4513*180/pi/360))/da)

% total amount angle steps needed to achieve N1 resolution

N2=N1*amount_full_turns+amount_da_rem

a0=angle(refl_ZL)
a=double(linspace(a0,2*imag(g_L),N2));

refl_mod=linspace(refl_ZL,refl_ZL*exp(-real(g_L)),N2)
refl_in=refl_mod(end).*cos(a(end))+1j*refl_mod(end).*sin(a(end))

plot(ax,real(refl_ZL),imag(refl_ZL),...
                                              'o','Color',[1 0 0],...
                                              'LineWidth',2,...
                                              'MarkerEdgeColor','b',...
                                              'MarkerFaceColor',[.8 .2 .2],...
                                              'MarkerSize',7)                                     % ZL 
                                                                    
 plot(ax,refl_mod(end).*cos(a(end)),refl_mod(end).*sin(a(end)),...
                                               'o','Color',[1 0 0],...
                                               'LineWidth',2,...
                                               'MarkerEdgeColor','b',...
                                               'MarkerFaceColor',[.2 .8 .2],...
                                               'MarkerSize',7)                                   % Zin at generator 
    
[x_swr,y_swr]=Smith_plotGammaCircle(ax,ZL,Z0,[1 .4 .4])

Zin_gen=Z0*(1+refl_in)/(1-refl_in)

[x_swr,y_swr]=Smith_plotGammaCircle(ax,Zin_gen,Z0,[.2 .1 .2])

SWR_Load=(1+abs(refl_ZL))/(1-abs(refl_ZL))
SWR_gen=(1+abs(refl_in))/(1-abs(refl_in))

str_swr_load=['SWR Load = ' num2str(SWR_Load)]
str_swr_gen=['SWR Generator = ' num2str(SWR_gen)]

legend(ax,'ZL','Zgen',str_swr_load,str_swr_gen)   

%%

% legend(ax,'off')

% repeat ZL Zgen points on new Smith chart to avoid the SWR circles to
% overlap the spiral that is now the imedance path across the Smith chart.

figure;sm1=smithchart; ax=gca; hold all; 

plot(ax,real(refl_ZL),imag(refl_ZL),...
                                              'o','Color',[1 0 0],...
                                              'LineWidth',2,...
                                              'MarkerEdgeColor','b',...
                                              'MarkerFaceColor',[.8 .2 .2],...
                                              'MarkerSize',7)                                     % ZL 
                                                                    
plot(ax,refl_mod(end).*cos(a(end)),refl_mod(end).*sin(a(end)),...
                                               'o','Color',[1 0 0],...
                                               'LineWidth',2,...
                                               'MarkerEdgeColor','b',...
                                               'MarkerFaceColor',[.2 .8 .2],...
                                               'MarkerSize',7)                                   % Zin at generator 
    

Smith_plotRefLine2PhaseCircle(ax,ZL,Z0,[.6 1 .6])
Smith_plotRefLine2PhaseCircle(ax,Zin_gen,Z0,[.6 1 .6])

plot(ax,refl_mod.*cos(a),refl_mod.*sin(a),'-','LineWidth',.5,'Color',[1 0 0]) 

axis([-0.4008   0.3832  -0.4103  0.3737])    

                                                      
                                                                    
% Comment 4: error in solutions manual, the value of tanh(gamma*L) does not correspond to the tanh of the calculated gamma*L
% It's actually 'amplified' by 1/exp(alpha*L) but it's not mentioned anywhere.
tanh(g_L)
tanh(real(g_L)+1j*(imag(g_L)-4*pi))


Vgen*Zin/(Zin+Zgen)*1/(exp(1j*imag(gamma_L_rad))+refl_Load*exp(-1j*imag(gamma_L_rad)))
abs(Vgen*Zin/(Zin+Zgen)*1/(exp(1j*imag(gamma_L_rad))+refl_Load*exp(-1j*imag(gamma_L_rad))))

Vgen*Z0/(Z0+Zgen)*exp(-1j*imag(gamma_L_rad))/(1-refl_Load*refl_gen*exp(-1j*2*imag(gamma_L_rad)))
abs(Vgen*Z0/(Z0+Zgen)*exp(-1j*imag(gamma_L_rad))/(1-refl_Load*refl_gen*exp(-1j*2*imag(gamma_L_rad))))





save(str65)   
        
end
