% This is a function file list functiosn to create operators that can be used to calulate measurable quantities for the
% spin system by using the eiegn values and vectors found using the hamiltonian functions. 
%
% These functions will
% 1. create the operator for any given pair of spins in a chain or cluster.
% 2. preform that operation for any given pair projected 100% on to any state for a range of states.
% 3. evaluates the partition function.
% 4. gives the correlation function for the spin system at any temperture using boltzmann distubaution across all states.
% 5. caluclates the ground state correlation function for the spin system.
%%%%%%%%%
% 6. This function file also calculates the ground state magnetsation
% 7. And caluclate the magnetisation for any temperature using boltzmann distubation.
%
% note: tempertaure and field values are in terms of the interaction energy between spins.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%
% FUNCTION: corrPair 
%           Constructs a realspace correlator operator for a given pair in a finite system.    
%           To be called in by pureCorrPoint
% INPUT:... delta   Positive Integer. Position i
%           deltap  Positive Integer. Position j
%           sigma_1 and sigma_2 pauli matrices. what type of correlator do you want? xx, xy etc
%           N       Positive Integer. size of Finite system.
% OUTPUT... correlatorPair   matrix (2^N,2^N) for a set correlation between i and j sites. 
function correlatorPair=corrPair(delta,deltap,sigma_1,sigma_2,N)

corr_1=zeros((2^N),(2^N));   % initiate 
corr_2=zeros((2^N),(2^N));
 

      corr_1(:,:)=kron(eye(2^(delta-1),2^(delta-1)),(kron(0.5*sigma_1,eye(2^(N-delta),2^(N-delta)))));        % position i
    
      corr_2(:,:)=kron(eye(2^(deltap-1),2^(deltap-1)),(kron(0.5*sigma_2,eye(2^(N-deltap),2^(N-deltap)))));       % position j

         correlatorPair(:,:)=corr_1(:,:)*corr_2(:,:);   

end
%%%%%%%%%%%
%%%%%%%%%%%
% FUNCTION: pureCorr
%           Evaluates the realspace correlation value for a given pair of sites i and j. 
%           assumes the system is in a pure state. to then be fed into finiteCorr for finite 
%           temperatture results. Groundstate given here for T=0 when using states_array(:1) %           only. Calls in operator function (above).
% INPUT:... delta   Positive Integer. Position i
%           deltap  Positive Integer. Position j
%           sigma_1 and sigma_2 pauli matrices. what type of correlator do you want? xx, xy etc
%           N       Positive Integer. size of Finite system.
%           states_array. matrix (2^N,2^N) are the eigenstates calculated from diagonalising 
%           the Hamiltonian. Using function DiagHmat. 
%           NOTE: must take states_array as a matrix, can't take one vector at a time. this 
%           feature is integrated into finiteCorr function as to not require this one. 
% OUTPUT:...pureCorrPoint   vector (2^N) correlator value per state for a given pair of sites. 
function pureCorrPoint=pureCorr(delta,deltap,sigma_1,sigma_2,N,states_array)

   correlatorPair=corrPair(delta,deltap,sigma_1,sigma_2,N);   % compute operator pair
                                           % compute expectation value for each state

      for ii=1:2**N
               pureCorrPoint(ii,:)=(conj(states_array(:,ii)'))*(correlatorPair*states_array(:,ii));
      endfor

end  
%%%%%%%%%%%
%%%%%%%%%%%
% FUNCTION: Z
%           Evaluates the Partition Function for a given set of eigen values.
% INPUT:... T               Positive number. Temperature as a percentage of J.
%           energies_array  vector 2^N.      Eigenvalues of the system.
% OUTPUT:...partitionFunction Positive number. 
function partitionFunction=Z(T,energies_array)

   [x,y] = size (energies_array);

        partitionFunction = 0;
        for ii=1:x
                partitionFunction = partitionFunction + exp(-(energies_array(ii)/T));
        endfor


end
%%%%%%%%%%%% 
%%%%%%%%%%%% 
% FUNCTION: finiteCorr
%           Evaluates the correlators in pair for a given temperature T. calls in pureCorr.
% INPUT:... delta   Positive Integer. Position i
%           deltap  Positive Integer. Position j
%           sigma_1 and sigma_2 pauli matrices. what type of correlator do you want? xx, xy etc
%           N       Positive Integer. size of Finite system.
%           states_array. matrix (2^N,2^N) are the eigenstates calculated from diagonalising 
%           the Hamiltonian. Using function DiagHmat.
%           T               Positive number. Temperature as a percentage of J.
%           energies_array  vector 2^N.      Eigenvalues of the system. 
%           partZ   Positive number. Partition Function
% OUTPUT:...finiteCorrelator. A number. the realspace correlation value for a given pair difined 
%           by delta and deltap position at a given temperature for teh original set of 
%           parameters. Gamma, Delta, field etc.
function finiteCorrelator=finiteCorr(delta,deltap,sigma_1,sigma_2,N,states_array,energies_vec,T,partZ)


 correlatorPair=corrPair(delta,deltap,sigma_1,sigma_2,N);
%[x,y]=size(energies_vec);

finCorr=0;
	for ii=1:2**N 
	   pureCorrPoint(ii,:)=(conj(states_array(:,ii)'))*(correlatorPair*states_array(:,ii));
       value=(pureCorrPoint(ii)*exp((-energies_vec(ii))/T))/partZ;  
       finCorr=finCorr+value;

	endfor
   finiteCorrelator=finCorr;
end
%%%%%%%%%%%%
%%%%%%%%%%%%
% FUNCTION: gsCorr
%           Evaluates the correlators in pair for a given temperature T. calls in pureCorr.
% INPUT:... delta   Positive Integer. Position i
%           deltap  Positive Integer. Position j
%           sigma_1 and sigma_2 pauli matrices. what type of correlator do you want? xx, xy etc
%           N       Positive Integer. size of Finite system.
%           state   Vector (2**N) the groundstate
% OUTPUT:...gsCorrelator. Groundstate correlator value for pair of delta deltap.
function gsCorrelator=gsCorr(delta,deltap,sigma_1,sigma_2,N,state)

   correlatorPair=corrPair(delta,deltap,sigma_1,sigma_2,N);
	
   gsCorrelator=(conj(state'))*(correlatorPair*state);

end
%%%%%%%%%%%%
%%%%%%%%%%%%
%% FUNCTION: magOp
%            This function constructs the magnestisation operator for a given system size, site and 
%            direction
% INPUT:...  sigma                 pauli matrix, identifies direction to measure magnetism along
%            N                     System size
%            delta                 a given site
% OUTPUT:... Magnetisation Operator matrix
%%%%%%%%%%%%
function operatorMag=magOp(sigma,N,delta)

   operatorMag=zeros(2^N,2^N);

   operatorMag=kron(eye(2^(delta-1)),(kron((sigma),eye(2^(N-delta),2^(N-delta)))));
	     

endfunction
%%%%%%%%%%%%
%%%%%%%%%%%%
%% FUNCTION: magExGS
%            This function evaluates the expectation value of the magnetisation at a given site 
%            for the groundstate 
% INPUT:...  sigma                 pauli matrix, identifies direction to measure magnetism along
%            N                     System size
%            delta                 a given site
%            state                groundstate
% OUTPUT:... Magnetisation at a given site. sum over all sites for net magnetisation.
%%%%%%%%%%%%
function siteMag=magSite(sigma,N,delta,state)

   operatorMag=magOp(sigma,N,delta);
   siteMag=(state')*(operatorMag*(state));

endfunction
%%%%%%%%%%%%
%%%%%%%%%%%%
%% FUNCTION: magT
%            This function evaluates the expectation value of the magnetisation at a given site 
%            for the groundstate 

% INPUT:...  sigma                 pauli matrix, identifies direction to measure magnetism along
%            N                     System size
%            delta                 a given site
%            state                groundstate
% OUTPUT:... Magnetisation at a given site. sum over all sites for net magnetisation.
%%%%%%%%%%%%
function siteMagT=magT(sigma,N,delta,states_array,energies_vec,T,partZ)

   operatorMag=magOp(sigma,N,delta);
   magtemp=0;
   for ii=1:2**N
   siteMag(ii)=(states_array(:,ii)')*(operatorMag*(states_array(:,ii)));
   value=(siteMag(ii)*exp((-energies_vec(ii))/T))/partZ;
   magtemp=magtemp+value;
   end
   siteMagT=magtemp;

endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
