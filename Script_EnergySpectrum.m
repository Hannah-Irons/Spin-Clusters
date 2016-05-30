%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This script file is an example in using the hamiltonain function file to build a spin cluster 
% hamiltonian matrix and solve it for a range of applied field values to explore its energy spectrum.
% Orginally scripted in octave.
%
% This Script will;
% 1. Set the paramters of the Hamiltonian.
% 2. Set the range and step size of the applied field.
% 3. Calls in the appropraite functions to build and solve the Hamiltonian matrices dependent on field.
% 4. Extracts the energy values for all energy levels for each value of the field.
% 5. Outputs the data in a saved data file. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

source "Hamiltonains.m" % source in function file that builds and solves the Hamiltonains

% set the parameters

N=4; % system size
Delta=0.5; % out of plane anisotropy
Gamma=0.5; % in plane anisotropy
h=[0:0.01:2]; % field range and step size, in vector format.
 
% This is for a four spin system XYZ heisenberg model.

% Use a loop to scan across the field range, same approach used to scan across any varaible, i.e. aniostropy.
   
   for hh=1:length(h);
      for Level=1:2^N
         [states_array,energies_array]=DiagpbHmat(N,Gamma,Delta,h(hh)); % we are using Periodic Boundary Conditions
         Eng(hh,Level)=energies_array(Level);
      endfor
   endfor
data=[h',Eng];

save "N4G02D05EnergySpectrum.dat" data; % I tend to save data files with the paramter information in the title.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
