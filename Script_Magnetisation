% This is an example script to calculate Magnetisation of a spin chain or cluster. Originally written for Octave.
%
% This script will;
% 1. calculate the magnetisation of a spin chain in the ground state.
% 2. calculate the magnetisation of a spin cluster (with periodic boundary conditions)in the ground state.
% 3. calculate the magnetisation of a spin chain for a range of tempertures.
% 4. calculate the magnetisation of a spin cluster (with periodic boundary conditions)for a range of tempertures.
% 5. saves data to .dat files. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

source "Hamiltonians.m"
source "Operators"

N=8; Delta=0;     % set desired parameters. Delta=0 means XY-model
h=[0:0.01:1.5];   % range of field values with small step size
Gamma=[0:0.01:1];  % range of in plane anisotropy with smal step size
 
sigmax=[0,1;1,0];             % Pauli matrices, can measure for a given direction of magnetisation
sigmay=[0,-i;i,0];
sigmaz=[1,0;0,-1];

   for gg=1:length(Gamma)
      for hh=1:length(h);
         [statesPB,energiesPB]=DiagpbHmat(N,Gamma(gg),Delta,h(hh));
         [states,energies]=DiagHmat(N,Gamma(gg),Delta,h(hh));   
         GSPB=statesPB(:,1);
         GS=states(:,1);;
         delta=1;
            OpMagPB=magSite(sigmaz,N,delta,GSPB);
            OpMag=magSite(sigmaz,N,delta,GS);
            MagPB(gg,hh)=OpMagPB;                     % can stop and save data here for surface plots for anistropy and field
            Mag(gg,hh)=OpMag;        
         
        DiffMag=Mag-MagPB;                     % we take a difference because I was looking for sensitivity to edge conditions
      endfor
   endfor
data=DiffMag;

save "N8GSmagnetisationDiff.dat" data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% or for a range of temperatures
clear all

source "Hamiltonians.m"   % source required function files
source "Operators.m"

N=4; Delta=0; h=[0:0.01:2]; T=[0.01,0.05,0.1,0.2];  % set parameters

Gamma=0.5; 

sigmax=[0,1;1,0];                        % Pauli Matrices
sigmay=[0,-i;i,0];
sigmaz=[1,0;0,-1];
delta=1;

for hh=1:length(h)

      [statesPB,energiesPB]=DiagpbHmat(N,Gamma,Delta,h(hh));
      [states,energies]=DiagHmat(N,Gamma,Delta,h(hh));
      GSPB=statesPB(:,1);
      GS=states(:,1);;
      OpMagPB=magSite(sigmaz,N,delta,GSPB);
      OpMag=magSite(sigmaz,N,delta,GS);
      MagPBgs(hh)=OpMagPB;  
      Maggs(hh)=OpMag;        

endfor


for tt=1:length(T)
   for hh=1:length(h)
         [statesPB,energiesPB]=DiagpbHmat(N,Gamma,Delta,h(hh));
         [states,energies]=DiagHmat(N,Gamma,Delta,h(hh)); 
         partZPB=Z(T(tt),energiesPB); 
         partZ=Z(T(tt),energies); 
         OpMagPB=magT(sigmaz,N,delta,statesPB,energiesPB,T(tt),partZPB);
         OpMag=magT(sigmaz,N,delta,states,energies,T(tt),partZ);
         MagPB(hh,tt)=OpMagPB;  
         Mag(hh,tt)=OpMag;
   endfor
endfor

         
data=[h',MagPBgs',Maggs',MagPB,Mag];              % a variety of plot choices from this data

save "N4G05magnetsationTemp.dat" data;
