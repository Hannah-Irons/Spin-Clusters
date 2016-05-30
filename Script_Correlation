% The is an example script file to calculate teh correlation functions of our spin systems
%
% This script calculates;
% 1. The real space correlation functions in the XX spin direction for a spin chain in the ground state
% 2. The real space correlation function in the XX spin direction for a spin cluster at some finite temperture.
%
clear all;

source "Hamiltonians.m"  % source required function files
source "Operators.m"

N=8; Gamma=0.2; Delta=0; % Delta is zero so XY-model
h=[0:0.05:2];    
sigmax=[0,1;1,0];        % measuring for correlation in the spin X direction only, can do any combination.

  for hh=1:length(h)
   [states_array,energies_array]=DiagHmat(N,Gamma,Delta,h(hh));     % can change this is periodic boundary condition function

   delta=1;
%   GSCorr=zeros(N,N);
%     for delta=1:N
     jj=1;
      for deltap=1:N 
        
          GSCorr(deltap,hh)=gsCorr(delta,deltap,sigmax,sigmax,N,states_array(:,1));      
      jj++;   
 %     endfor
     endfor
     endfor
   

R=[0:N-1];

data=[R',abs(GSCorr)]    % I take the absolute value because the system is antiferromagetic and I'm looking for overall flatness in teh correlation functions

save "obN8CorrelatorXXG02fieldrange.dat" data        % save the dta to a data file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% or for a finite tempertare correlator
clear all;

source "Hamiltonians.m"  % source required function files
source "Operators.m"

N=4; Gamma=0.2; Delta=0;     
sigmax=[0,1;1,0];
h=0.98;
T=0.01;                 % temperture is in terms of the interaction energy, can do a range of temps in this script

   [states_array,energies_vec]=DiagpbHmat(N,Gamma,Delta,h);    % this is for PBC can change function to do chain.
   partZ=Z(T,energies_vec);
    FinCorr=zeros(1,N);
    R=zeros(1,N);
    jj=1;
   delta=1;
%delta=1;
   for tt=1:length(T)          % if you want a range of temps
      for deltap=1:N 
          R(jj)=deltap-delta;
          FinCorr(delta,deltap)=finiteCorr(delta,deltap,sigmax,sigmax,N,states_array,energies_vec,T(tt),partZ);
          jj++;              
      end
   end

data=[R',abs(FinCorr(1,:)')]   % save for the first temp, but in this case there is only one temp.

save "N4CorrelationXXT001G02H098.dat" data  % create data file.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
