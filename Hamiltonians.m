% The following file is a function file used to create the Heisenberg Hamiltonian for small finite-sized 
% spin 1/2 systems with open or periodic boundary conditions.These functions were used and written in Octave. 
%
% This Function file contains
% 1. Functions that build up matrices to describe the spin interactions S(alpha)_1S(alpha)_2 dependent on the SIZE of the system.
% 2. Builds the Hamiltonian matrix uisng the spin constrauction and can apply a field in the z direction and add aniostropy in any direction
% 3. Will solve that Hamiltonain using exact diagonilastion for a range or field values or anisotropy.
% 4. Final output is an array of eigen values and corresponding eigen vectors for the chosen system. The first is always the ground state.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: spinConstruction. 
%           This function creates spin matrices that describe the interaction and field
%           componants in the Hamiltonian for a spin chain.
% INPUT:... sigma     Pauli Matrix (2x2). for x y z
%           N         Positive Interger. Size of chain. has limiting factors based on active memory to hold teh size of the matrix.
%
% OUTPUT:...Spin Matrix (2^N,2^N) for S(alpa)i*S(alpha)i+1 (interaction between neighbouring spins)
%           Field Matrix (2^N,2^N) transverse field componant.         
%%%%%%%%
function [Spin_matrix,F_matrix]=spinConstruct(sigma,N)

% spin component;

    Spin_matrix=zeros((2^N),(2^N));    % initiate

    for ii=1:N-1                       % for each pair of interactions for an OPEN CHAIN

            interactionMat=kron(eye(2^(ii-1),2^(ii-1)),(kron(sigma,eye(2^(N-ii),2^(N-ii)))))*kron(eye(2^(ii),2^(ii)),(kron(sigma,eye(2^(N-ii-1),2^(N-ii-1)))));		
	        Spin_matrix=Spin_matrix+interactionMat;

    endfor

% field componant. NOTE: this assumes an applied field in the Z direction

   F_matrix=zeros(2^N,2^N);

   for jj=1:N

	    fieldElement=kron(eye(2^(jj-1)),(kron(([1,0;0,-1]),eye(2^(N-jj),2^(N-jj)))));
	    F_matrix=F_matrix+fieldElement;

   endfor
end
%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%
% FUNCTION: PBspinConstruction : same as above for PERIODIC BOUNDARY CONDITIONS
%           This function creates spin matrices that describe the interaction and field
%           componants in the Hamiltonian 
% INPUT:... sigma     Pauli Matrix (2x2). for x y z
%           N         Positive Interger. Size of chain. has limiting factors. N=13 is max.
%
% OUTPUT:...Spin Matrix (2^N,2^N) for S(alpa)i*S(alpha)i+1
%           Field Matrix (2^n,2^n) transverse field componant.         
%%%%%%%%
%%%%%%%%%
function [pbSpin_matrix,F_matrix]=PBspinConstruct(sigma,N)

Spin_matrix=zeros((2^N),(2^N));

for ii=1:N-1	

	interactionMat=kron(eye(2^(ii-1),2^(ii-1)),(kron(sigma,eye(2^(N-ii),2^(N-ii)))))*kron(eye(2^(ii),2^(ii)),(kron(sigma,eye(2^(N-ii-1),2^(N-ii-1)))));	
                  
	Spin_matrix=Spin_matrix+interactionMat;

endfor
  pb=kron(eye(2^(N-1),2^(N-1)),sigma)*kron(sigma,eye(2^(N-1),2^(N-1)));
  pbSpin_matrix=Spin_matrix+pb;

% field componant. NOTE: this assumes a transverse field in the Z direction

   F_matrix=zeros(2^N,2^N);

   for jj=1:N

	    fieldElement=kron(eye(2^(jj-1)),(kron(([1,0;0,-1]),eye(2^(N-jj),2^(N-jj)))));
	    F_matrix=F_matrix+fieldElement;

   endfor

end	
%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%
%%%%%%%%%
% FUNCTION: Hmat
%          This function calls in function spinConstruct three times for SxSx, SySy and SzSz.
%          This means that it can be applied to the xyz model also. It takes the spin 
%          interaction matrices and field componant and evaluates them in the Heisenberg
%          Hamiltonian.
% INPUT:...N      Positive Integer. Size of finite Chain.
%          Gamma  Positive number 0<Gamma<1. Anisotropy for x and y
%          Delta  Positive number. Anisotropy for z. taken as zero for xy model.
%          h      Positive number. Effective field. i.e h=1 is the critical field.
% OUTPUT...The Hamiltonian matrix (2^N,2^N)
function Hmatrix=Hmat(N,Gamma,Delta,h)
% combime for the Hamiltonian. INPUT: Gamma, Delta and field call in above function, require input N

   [SxSx,F_matrix]=spinConstruct([0,1;1,0],N);   % call in for x,y,z
   [SySy,F_matrix]=spinConstruct([0,i;-i,0],N);
   [SzSz,F_matrix]=spinConstruct([1,0;0,-1],N);

   Vx=SxSx*(1+Gamma);     % add anisotropy for x,y and z 
   Vy=SySy*(1-Gamma);     % gamma is the in-plane anisotropy
   Vz=SzSz*Delta;         % Delta is the out of plane anisotropy

   Hmatrix=0.25*(Vx+Vy+Vz)-0.5*h*F_matrix;       % construct heisenberg hamiltonian matrix with field. 

end
%%%%%%%%%%
%%%%%%%%%%
% FUNCTION: pbHmat: same as above for PERIODIC BOUNDARY CONDITIONS
%          This function calls in function spinConstruct threee times for SxSx, SySy and SzSz.
%          This means that it can be applied to the xyx model also. It takes the spin 
%          interaction matrices and field componant and evaluates them in the Heisenberg
%          Hamiltonian.
% INPUT:...N      Positive Integer. Size of finite Chain.
%          Gamma  Positive number 0<Gamma<1. Anisotropy for x and y
%          Delta  Positive number. Anisotropy for z. taken as zero for xy model.
%          h      Positive number. Effective field. i.e h=1 is the critical field.
% OUTPUT...The Hamiltonian matrix (2^N,2^N)
function Hmatrix=pbHmat(N,Gamma,Delta,h)
% combime for the Hamiltonian. INPUT: Gamma, Delta and field call in above function, require input N

   [SxSx,F_matrix]=PBspinConstruct([0,1;1,0],N);   % call in for x,y,z
   [SySy,F_matrix]=PBspinConstruct([0,i;-i,0],N);
   [SzSz,F_matrix]=PBspinConstruct([1,0;0,-1],N);

   Vx=SxSx*(1+Gamma);     % add anisotropy for x,y and z
   Vy=SySy*(1-Gamma);
   Vz=SzSz*Delta;

   Hmatrix=(1/4)*(Vx+Vy+Vz)-(1/2)*h*F_matrix;       % construct hamiltonian matrix with field. 

end
%%%%%%%%%%
%%%%%%%%%%
% FUNCTION: DiagHmat
%           This calls in Hmat so would be the only required function to call by a script to 
%           generate eigenvalues and eigen vecors for a system for a given range of h. h can be 
%           a singular number or a vector. It does this through exact diagonilisation
% INPUT:... N      Positive Integer. Size of Finite chain
%           Gamma  Positive number. 0<Gamma<1. Anisotropy for x and y
%           Delta  Positive number. Anisotropy for z. taken as zero for xy model.
%           h_vec  Positive number or vector. Effective field. i.e h=1 is the critical field.
% OUTPUT... eVec   Matrix (2^N,2^N) of eigen vectors.
%           eVal   Vector 2^N of eigen values.
function [eVec,eVal]=DiagHmat(N,Gamma,Delta,h_vec)
for ii=1:length(h_vec) % h_vec can be a single field value of a vector range of them.
	
	h=h_vec(ii); 
	[e_vector,e_value]=eig(Hmat(N,Gamma,Delta,h));
	eVec(:,:,ii)=e_vector;
	eVal(:,ii)=diag(e_value);
		
endfor
end
%%%%%%%%%%
%%%%%%%%%%
% FUNCTION: DiagpbHmat: same as above for PERIODIC BOUNDARY CONDITIONS
%           This calls in Hmat so would be the only required function to call by a script to 
%           generate eigenvalues and eigen vecors for a system for a given range of h. h can be 
%           a singular number or a vector. It does this through exact diagonilisation
% INPUT:... N      Positive Integer. Size of Finite chain
%           Gamma  Positive number. 0<Gamma<1. Anisotropy for x and y
%           Delta  Positive number. Anisotropy for z. taken as zero for xy model.
%           h_vec  Positive number or vector. Effective field. i.e h=1 is the critical field.
% OUTPUT... eVec   Matrix (2^N,2^N) of eigen vectors.
%           eVal   Vector 2^N of eigen values.
function [eVec,eVal]=DiagpbHmat(N,Gamma,Delta,h_vec)
for ii=1:length(h_vec) % h_vec can be a single field value of a vector range of them.
	
	h=h_vec(ii); 
	[e_vector,e_value]=eig(pbHmat(N,Gamma,Delta,h));
	eVec(:,:,ii)=e_vector;
	eVal(:,ii)=diag(e_value);
		
endfor
end
%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
