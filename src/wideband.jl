load("Sims")
using Sims

##
## Wideband line model simulation
##
##
## Based on Matlab code from:
##
## Octavio Ramos-Leaños, Jose Luis Naredo and Jose Alberto
## Gutierrez-Robles, "An Advanced Transmission Line and Cable Model in
## Matlab for the Simulation of Power-System Transients," MATLAB - A
## Fundamental Tool for Scientific Computing and Engineering
## Applications - Volume 1, Edited by Vasilios N. Katsikis, ISBN
## 978-953-51-0750-7, InTech,2012.
## http://dx.doi.org/10.5772/48530
##
## Copyright 2012 Octavio Ramos-Leaños, Jose Luis Naredo and Jose Alberto Gutierrez-Robles 
## Creative Commons Attribution License (http://creativecommons.org/licenses/by/3.0)
##
##
## Translated into Julia by Tom Short, EPRI, 2012
##


const j = 1.0im
mysort(a) = sort((x,y) -> abs(x) == abs(y) ? imag(x) < imag(y) : abs(x) < abs(y),a)

function LineParameters(Mu,Eo,Rsu,Geom,Ncon,Ns,w)
  # Function to compute the distances between conductor
  Dij,dij,hij =Height(Geom)
  Zg = zeros(Complex128, Ncon,Ncon,Ns)
  Zt = zeros(Complex128, Ncon,Ncon,Ns)
  Zc = zeros(Complex128, Ncon,Ncon,Ns)
  Yg = zeros(Complex128, Ncon,Ncon,Ns)
  Zcd = zeros(Complex128, Ncon,Ns)
  Zaf = zeros(Complex128, Ncon,Ns)
  P = (1./sqrt(j*w*Mu/Rsu)) # Complex depth
  Pmatrix = log(Dij./dij) # Potential Coeff. Matrix
  Pinv = inv(Pmatrix) # Inverse of Pmatrix
  # Loop to compute matrices at all frequencies
  for kl = 1:Ns
    # Geometric impedance
    Zg[:,:,kl] = (j*w[kl]*Mu/(2*pi))*Pmatrix
    # Earth impedance
    for km = 1:Ncon
      for kn = 1:Ncon
        if km == kn
          Zt[km,km,kl] = (j*w[kl]*Mu/(2*pi))*
          log(1+P[kl]./(0.5*hij[km,km]))
        else
          num = hij[km,kn]^2 + 4*P[kl]*hij[km,kn] +
          4*P[kl]^2 + dij[km,kn]^2
          den = hij[km,kn]^2 + dij[km,kn]^2
          Zt[km,kn,kl] = (j*w[kl]*Mu/(4*pi))*
          log(num/den)
        end
      end
    end
    # Geometric admittance
    Yg[:,:,kl] = (j*w[kl]*2*pi*Eo)*Pinv
  end
  # Conductor impedance
  for kd = 1:Ncon
    Rcon = Geom[kd,4] # conductor radii in m.
    Nhaz = Geom[kd,5] # # of conductor in bundle
    Rpha = Geom[kd,7] # Resistivity in Ohm-m.
    Zcd[kd,:] = (1/Nhaz)*Rpha./(pi.*Rcon.^2)
    Zaf[kd,:] = (1/Nhaz)*(1+j).*(1./(2.*pi.*Rcon)) .*
    sqrt(0.5.*w.*Mu.*Rpha)
    Zc[kd,kd,:] = sqrt(Zcd[kd,:].^2 + Zaf[kd,:].^2)[:]
  end
  # Outputs
  ZT = Zg + Zt + Zc  # Total impedance
  YT = Yg  # Total admittance
  (Zg,Zt,Zc,Yg,ZT,YT)
end

function YcFit(Yc,f,Ns,Ncon)
    # Trace of characteristic admittance matrix
    Ytrace = zeros(Complex128, Ns)
    for k = 1:Ns
      Ytrace[k] = trace(Yc[:,:,k])
    end
    Npol = 6 # Number of poles
    Ps = InitialPoles(f,Npol) # Set initial poles
    s = j*2*pi*f # Vector of values of variable "s"
    Ka=2 # 1.-Strictly proper, 2.-Proper, 3.-Improper
    local YcPoles
    for khg=1:20
      # Fit the trace of Yc (Poles)
      YcPoles=Poles(Ytrace,s,Ps,Ns,Ka)
      Ps=YcPoles
    end
    # Residues and constant term for Yc from poles of trace of Yc
    YcResidues = zeros(Complex128, Ncon, Ncon, Npol)
    YcConstant = zeros(Complex128, Ncon, Ncon)
    YcProportional = zeros(Complex128, Ncon, Ncon)
    for k = 1:Ncon
      for l = 1:Ncon
        Hs = Yc[k,l,:][:] # k-l term of admittance
        C,D,E = Residue(Hs,s,YcPoles,Ns,Ka)
        YcResidues[k,l,:] = C # k-l residues term
        YcConstant[k,l] = D # k-l constant term
        YcProportional[k,l]=E #k-l proportional term
      end
    end
    (YcPoles,YcResidues,YcConstant, YcProportional)
end

function ABYZLM(Z,Y,Lo,w)
  Lm, M = eig(Y*Z) # Eigenvalues of YZ
  Lm = diagm(Lm)
  Minv = inv(M) # Inverse of eigenvectors matrix
  Yc = inv(Z)*(M*sqrt(Lm)*Minv) # Characteristic Admittance
  Gamma = sqrt(diag(Lm)) # Propagation constants.
  Vm = w./imag(Gamma) # Modal Velocities
  Hmo = diag(expm(-sqrt(Lm)*Lo)) # Modal propag. Matrix H
  (Yc,Vm,Hmo)
end

function findrowmin(x::Matrix)
    nrow, ncol = size(x)
    res = zeros(Int, nrow)
    for i in 1:nrow
        res[i] = findmin(x[i,:])[2]
    end
    res
end

#
# Function to calculate the modal delay
#
# Input
#
#    H     ==> Modal H
#    f     ==> Vector of frequencies
#    long  ==> Line lenght
#    V     ==> Modal velocities
#    Ns    ==> Number of samples
#    Nc    ==> Number of conductors
#
# Outputs
#
#    taumin  ==> minimum delay of each mode
#
# Call the function
#
#    [taumin]=ModeDelay(H,f,long,V,Ns,Nc)
#
#

function ModeDelay(H,f,long,V,Ns,Nc)
  
  w = 2*pi*f  # Angular frequency rad/seg
  
  # Minimal value Hm and position Hp of each propagation mode
  Hm = min(abs(H[:,2:Ns]),(), 2)
  Hp = findrowmin(abs(H[:,2:Ns]))

  Td = zeros(Nc)
  phase1 = zeros(Nc)
  for k=1:Nc
      # Time delay
      Td[k]=long./V[k,Hp[k]]
      # Minimal phase angle
      k1 = log(abs(H[k,Hp[k]+1])/abs(H[k,Hp[k]-1]))
      k2 = log((w[Hp[k]+1])/(w[Hp[k]-1]))
      phase1[k]= (pi/2)*k1/k2
  end
  
  # Minimal phase delay in seconds
  taumin   = Td + phase1./w[Hp]
  taumin
end

#
# Idempotent matrices of Y*Z
#
# Input
#
#    Z      ==> Total impedance
#    Y      ==> Total admittance
#    L      ==> Line lenght
#    f      ==> Vector of frequencies
#    Delays ==> Modal delays
#    Hm     ==> Modal Hm
#
# Outputs
#
#    Idem ==> Idempotent asociated to lambda 1,2 and 3
#    Idem = [Id1 Id2 Id3]
#
# Call the function
#
#    [Idem] = HmIdem(Z,Y,L,f,Delays,Hm)
#
#
function HmIdem(Z,Y,L,f,Delays,Hm)

  s  = j*2*pi*f       # Vector of the variable "s"
  t1 = Delays[1]   # Time delay of mode 1
  t2 = Delays[2]   # Time delay of mode 2
  t3 = Delays[3]   # Time delay of mode 3
  
  Lambda,M = eig(Y*Z)   # Eigenvectors and eigenvalues of Y*Z
  Mi       = inv(M)      # Inverse of the eigenvectors
  
  C1 = M[:,1];   R1 = Mi[1,:]  # First column and row of the matrices M and Mi respectively
  C2 = M[:,2];   R2 = Mi[2,:]  # Second column and row of the matrices M and Mi respectively
  C3 = M[:,3];   R3 = Mi[3,:]  # Third column and row of the matrices M and Mi respectively
  
  Id1 = (C1*R1)*Hm[1]*exp(s*t1)  # First Idempotent
  Id2 = (C2*R2)*Hm[2]*exp(s*t2)  # Second Idempotent
  Id3 = (C3*R3)*Hm[3]*exp(s*t3)  # Third Idempotent
  
  Idem = cat(3, Id1, Id2, Id3)
  Idem
end


function HkFitTrace(Hm,Vm,ZL,YL,f,Ns, lenght,Ncon)
  # Minimum phase of each mode
  md = ModeDelay(Hm.',f,lenght,Vm.',Ns,Ncon)
  # Computing Idempotents
  HkIdem = zeros(Complex128, Ncon, Ncon, Ncon, Ns) 
  for k=1:Ns
    # Function to calculate Idempotents of Y*Z
    Hk = HmIdem(ZL[:,:,k],YL[:,:,k],lenght,f[k], md,Hm[k,:])
    HkIdem[:,:,:,k] = Hk # Idempotents
  end
  TraceHk = zeros(Complex128, 3, Ns)
  for m = 1:3
    for k=1:Ns
      TraceHk[m,k] = trace(HkIdem[:,:,m,k])
    end
  end
  s = j*2*pi*f # Vector of the variable "s"
  Ka =1 #1.-Strictly proper, 2.-Proper, 3.-Improper
  Npol = 5 # Number of poles
  ## local Ps
  Ps = InitialPoles(f,Npol) # Set the initial poles
  HkPoles = zeros(Complex128, 3, Npol)
  for m = 1:3
    Hk = TraceHk[m,:][:]
    for khg=1:10
      HkPol=Poles(Hk,s,Ps,Ns,Ka)
      Ps=HkPol
    end
    HkPoles[m,:]=Ps
  end
  # Residues for Idempotent matrices of
  # Hm from the poles of each trace.
  HkResidues = zeros(Complex128, Ncon, Ncon, 3, Npol)
  HkConstant = zeros(Complex128, Ncon, Ncon, 3)
  HkProportional = zeros(Complex128, Ncon, Ncon, 3)
  for m = 1:3
    for k = 1:Ncon
      for l = 1:Ncon
        Hs = HkIdem[k,l,m,:][:] # k-l term
        C,D,E = Residue(Hs,s,HkPoles[m,:],Ns,Ka)
        HkResidues[k,l,m,:] = C # k-l-m term
        HkConstant[k,l,m] = D # k-l-m constant
        HkProportional[k,l,m] = E # k-l-m prop
      end
    end
  end
  (HkPoles,HkResidues,HkConstant, HkProportional,md)
end


function Height(Geom)
  Ls = Geom[max(Geom[:,1]),1]
  Req = zeros(int(Ls))
  # Equivalent bundle radii
  k4 = sqrt(2*(Geom[:,6]/2).^2)
  for nc = 1: Ls
    if Geom[nc,5]==1
      Req[nc] = Geom[nc,4]
    else
      Req[nc] = (Geom[nc,4].*Geom[nc,5].*k4[nc].^(Geom[nc,5]-1)).^(1./Geom[nc,5])
    end
  end
  dij = zeros(int(Ls), int(Ls))
  hij = copy(dij)
  Dij = copy(dij)
  # Direct and image distances among conductors
  for xl = 1:Ls
    for yl = 1:Ls
      if xl==yl
        dij[xl,yl]=Req[xl]
        y1=Geom[yl,3]
        hij[xl,yl]=2*y1
        Dij[xl,yl]=hij[xl,yl]
      else
        x=abs(Geom[yl,2]-Geom[xl,2])
        y=abs(Geom[yl,3]-Geom[xl,3])
        dij[xl,yl]=sqrt(x^2 + y^2)
        y1=Geom[xl,3]
        y2=Geom[yl,3]
        hij[xl,yl]=y1+y2
        x=abs(Geom[yl,2]-Geom[xl,2])
        y=hij[xl,yl]
        Dij[xl,yl]=sqrt(x^2 + y^2)
      end
    end
  end
  (Dij,dij,hij)
end

function InitialPoles(f,Npol)
  even = int(floor(Npol/2)) # # of complex initial poles
  p_odd = Npol/2 - even # Auxiliary variable to check if number
  # of initial poles is odd
  disc = p_odd != 0 # 0 for even Nr of initial poles & 1 ¡V for
  # odd Nr.
  # Set a real pole in case of disc == 1
  if disc == 0 # Even Nr of initial poles
    pols = []
  else # Odd Nr of initial poles
    pols = [(max(f)-min(f))/2]
  end
  # Set the complex initial poles
  bet = linspace(min(f),max(f),even)
  for n=1:length(bet)
    alf=-bet[n]*1e-2
    pols=[pols, (alf-j*bet[n]), (alf+j*bet[n]) ]
  end
  pols
end

function Poles(Fs,s,Pi,Ns,Ka)
  Np = length(Pi) # Length of vector containing starting poles
  CPX = imag(Pi).!=0 # 0 for real pole and 1 for complex pole
  rp = 0 # Initialize the index for real poles
  cp = 0 # Initialize the index for complex poles
  RePole = Complex128[] # Initialize the vector of real poles
  CxPole=Complex128[]#Initialize the vector of complex poles
  # Loop to separate real poles and complex poles
  for k = 1:Np
    if CPX[k] # Real Pole
      rp = rp + 1
      append!(RePole, [Pi[k]])
    else # Complex pole
      cp = cp + 1
      append!(CxPole, [Pi[k]])
    end
  end
  Lambda = Pi
  RePole = mysort(RePole) # Sort real poles
  CxPole = mysort(CxPole) # Sort complex poles
  Lambda = [RePole, CxPole] # Concatenate poles
  I = diag(ones(1,Np)) # Unit matrix
  A = [] # Poles
  B = ones(Ns) # the weight factor
  C = [] # Residues
  D = zeros(1) # Constant term
  E = zeros(1) # Proportional term
  KQA = ones(Ns,1)
  cpx = imag(Lambda).!=0 # 0 if pole is real and 1 if pole is complex.
  dix = zeros(Int,Np) # Initializes vector of pole types
  if cpx[1] # If the first pole is complex
    dix[1]=1 # real part
    dix[2]=2 # imag part
    k=3 # continue dix for third position
  else
    k=2 # If the first pole is real continue dix for the second position
  end
  # complete the classification of the poles
  for m=k:Np
    if cpx[m]!=0 # If the pole is complex
      if dix[m-1]==1
        dix[m]=2 # If the previous position has the real part put 2
      # to identifies the imag part
      else
        dix[m]=1 # 1 for the real part of a complex pole
      end
    end
  end
  # Creates matriz A divided in four parts A = [A1 A2 A3 A4]
  # A1 = Dk
  # A2 = B.*ones(Ns,1)
  # A3 = B.*s
  # A4 = -Dk*Fs
  Dk=zeros(Complex128,Ns,Np) # Initialize matrix with zeros
  for m=1:Np # Iterative cycle for all poles
    if dix[m]== 0 # For a real pole
      Dk[:,m] = B./(s-Lambda[m])
    elseif dix[m]== 1 # For the real part
      Dk[:,m]=B./(s-Lambda[m]) +
      B./(s-Lambda[m]')
    elseif dix[m]== 2 # For the imag part
      Dk[:,m] = j.*B./(s-Lambda[m-1]) -
      j.*B./(s-Lambda[m-1]')
    end
  end
  # Creates work space for matrix A
  A1 = Dk
  A2 = B.*ones(Ns,1)
  A3 = B.*s
  A4 = zeros(Complex128,Ns,Np) 
  for col = 1:Np
    A4[:,col] = -(Dk[:,col].*Fs)
  end
  # Asigns values to A
  if Ka == 1
    A = [A1 A4] # Strictly proper rational fitting
  elseif Ka == 2
    A = [A1 A2 A4] # Proper rational fitting
  elseif Ka == 3
    A = [A1 A2 A3 A4] # Improper rational fitting
  else
    error("Ka need to be 1, 2 or 3")
  end
  # Creates matrix b = B*Fs
  b = B.*Fs
  # Separating real and imaginary part
  Are = real(A) # Real part of matrix A
  Aim = imag(A) # Imaginary part of matrix A
  bre = real(b) # Real part of matrix b
  bim = imag(b) # Imaginary part of matrix b
  An = [Are, Aim] # Real and imaginary part of A
  bn = [bre, bim] # Real and imaginary part of b
  # Routine to applies the Euclidian norm to An
  Xmax, Ymax = size(An)
  Euclidian = zeros(Ymax)
  for col=1:Ymax
    Euclidian[col]=norm(An[:,col],2)
    An[:,col]=An[:,col]./Euclidian[col]
  end
  # Solving system
  Xn = An\bn
  Xn = Xn./Euclidian
  # Put the residues into matrix C
  if Ka == 1
    C = Xn[Np+1:Ymax] # Strictly proper fitting
  elseif Ka == 2
    C = Xn[Np+2:Ymax] # Proper rational fitting
  elseif Ka == 3
    C = Xn[Np+3:Ymax]# Improper rational fitting
  else
    disp("Ka need to be 1, 2 or 3")
  end
  # C complex when the residues are complex
  C = complex128(C)
  for m=1:Np
    if dix[m]==1
      alfa = C[m] # real part of a complex pole
      betta = C[m+1] # imag part of a complex pole
      C[m] = alfa + j*betta # the complex pole
      C[m+1] = alfa - j*betta # the conjugate
    end
  end
  # Now calculate the zeros for sigma
  BDA = zeros(Complex128, Np, Np)
  KQA = ones(Complex128, Np)
  # Loop to calculate the zeros of sigma which are the new poles
  for km = 1:Np
    if dix[km]== 0 # For a real pole
      BDA[km,km] = Lambda[km]
    elseif dix[km]== 1 # For a cp with - imag part
      BDA[km,km] = real(Lambda[km])
      BDA[km,km+1] = imag(Lambda[km])
      KQA[km] = 2.
      Aux = C[km]
      C[km] = real(Aux)
    elseif dix[km]== 2 # For a cp with + imag part
      BDA[km,km] = real(Lambda[km])
      BDA[km,km-1] = imag(Lambda[km])
      KQA[km] = 0.
      C[km] = imag(Aux)
    end
  end
  ZEROS = BDA - KQA*C'
  POLS,dum = eig(ZEROS)
  #Forcing (flipping) unstable poles to make them stable
  uns = real(POLS).>0
  POLS[uns] = POLS[uns]-2*real(POLS[uns])
  # Sort poles in ascending order. First real poles and then complex poles
  CPX = imag(POLS).!=0 # Set to 0 for a real pole and to1 for a
  #complex pole
  rp = 0 # Initialize index for real poles
  cp = 0 # Initialize index for complex poles
  RePole = Complex128[] # Initialize the vector of real poles
  CxPole = Complex128[] # Initialize the vector of cp
  # Loop to separate real and complex poles
  for k = 1:Np
    if CPX[k] == 0 # Real Pole
      rp = rp + 1
      append!(RePole, [POLS[k]])
    elseif CPX[k] == 1 # Complex pole
      cp = cp + 1
      append!(CxPole, [POLS[k]])
    end
  end
  RePole = mysort(RePole) # Sort real poles
  CxPole = mysort(CxPole) # Sort complex poles
  # For conjugate pairs store first the one with positive imag part
  CxPole = (CxPole.')'
  NewPol = [RePole, CxPole]
  NewPol
end

function Residue(Fs,s,Pi,Ns,Ka)
  Np = length(Pi)
  CPX = imag(Pi).!=0 
  rp = 0 # Initialize the index for real poles
  cp = 0 # Initialize the index for complex poles
  RePole = Complex128[] # Initialize the vector of real poles
  CxPole = Complex128[] # Initialize the vector of cp
  # Loop to separate real and complex poles
  for k = 1:Np
    if !CPX[k] # Real Pole
      rp = rp + 1
      append!(RePole, [Pi[k]])
    else # Complex pole
      cp = cp + 1
      append!(CxPole, [Pi[k]])
    end
  end
  RePole = mysort(RePole) # Sort real poles
  CxPole = mysort(CxPole) # Sort complex poles
  CxPole = (CxPole.')'
  Lambda = [RePole, CxPole]
  I = diagm(ones(Np)) # Unit diagonal matrix
  A = [] # Poles
  B = ones(Ns) # weight factor
  C = [] # Residues
  D = zeros(1) # Constant term
  E = zeros(1) # Proportional term
  cpx = imag(Lambda).!=0 # 0 for rp and 1 for cp
  dix = zeros(Int,Np) # Vto identifies poles
  if cpx[1] # If the first pole is complex
    dix[1]=1 # put 1 in dix[1] for the real part
    dix[2]=2 # put 2 in dix[2] for the imag part
    k=3 # continue dix for the third position
  else
    k=2 # If the first pole is real continue dix for the second
  # position
  end
  # complete classification of the poles
  for m=k:Np
    if cpx[m] # If the pole is complex
      if dix[m-1]==1
        dix[m]=2 # If the previous position has the real part, set to # 2 to identify the imag part
      else
        dix[m]=1 # put 1 for the real part of a cp
      end
    end
  end
  # Output matrices:
  Dk=zeros(Complex128,Ns,Np)
  for m=1:Np
    if dix[m]==0 # Real pole
      Dk[:,m] = B./(s-Lambda[m])
    elseif dix[m]==1 # Complex pole, 1st part
      Dk[:,m] = B./(s-Lambda[m]) + B./(s-Lambda[m]')
    elseif dix[m]==2 # Complex pole, 2nd part
      Dk[:,m] = j.*B./(s-Lambda[m-1]) - j.*B./(s-Lambda[m-1]')
    end
  end
  # Creates work space for matrices A and b
  AA1=Dk
  AA2=B.*ones(Ns)
  AA3=B.*s
  if Ka == 1
    AA = [AA1] # Strictly proper rational fit
  elseif Ka == 2
    AA = [AA1 AA2] # Proper rational fit
  elseif Ka == 3
    AA = [AA1 AA2 AA3] # Improper fit
  else
    error("Ka must be 1, 2 or 3")
  end
  bb = B.*Fs
  AAre = real(AA) # Real part of matrix A
  AAim = imag(AA) # Imaginary part of matrix A
  bbre = real(bb) # Real part of matrix b
  bbim = imag(bb) # Imaginary part of matrix b
  AAn = [AAre, AAim] # Real and imag part of A
  bbn = [bbre, bbim] # Real and imag part of b
  (Xmax, Ymax) = size(AAn)
  Euclidian = zeros(Ymax)
  for col=1:Ymax
    Euclidian[col]=norm(AAn[:,col],2)
    AAn[:,col]=AAn[:,col]./Euclidian[col]
  end
  # Solving system X
  Xxn=AAn\bbn
  X=Xxn./Euclidian
  # Putting residues into matrix C
  C=complex128(X[1:Np])
  # C is complex when the residues are complex
  for m=1:Np
    if dix[m]==1
      alfa = C[m] # real part of a complex pole
      betta = C[m+1] # imag part of a complex pole
      C[m] = alfa + j*betta # the complex pole
      C[m+1] = alfa - j*betta # the conjugate
    end
  end
  # Outputs
  if Ka == 1
    A = Lambda.' # Poles
    C = C # Residues
    D = 0 # Constant term
    E = 0 # Proportional term
  elseif Ka == 2
    A = Lambda.' # Poles
    C = C # Residues
    D = X[Np+1] # Constant term
    E = 0 # Proportional term
  elseif Ka == 3
    A = Lambda.' # Poles
    C = C # Residues
    D = X[Np+1] # Constant term
    E = X[Np+2] # Proportional term
  end
  (C,D,E)
end


  
function Wideband(Vx::ElectricalNode, Vy::ElectricalNode,
                  geom::Matrix, lenght::Real)

    # setup
    Ncon = int(geom[max(geom[:,1]),1]) # # of cond
    Mu = 4*pi*1E-7 # Henry's/meters
    Eo = (1/(36*pi))*1E-9 # Farads/meters
    Rsu = 100 # Earth resistivity Ohm-m
    Ns = 500 # Number of samples
    f = logspace(-2.0, 6.0, Ns) # Vector of log spaced Frequencies
    w = 2*pi*f # Vector of freqs in radian/sec.
    
    # Per unit length parameters
    Zg,Zt,Zc,Yg,ZL,YL = LineParameters(Mu,Eo,Rsu,geom,Ncon,Ns,w)
    
    # Modal Parameters
    Yc = zeros(Complex128, Ncon,Ncon,Ns)
    Vm = zeros(Ns,Ncon)
    Hm = zeros(Complex128, Ns,Ncon)
    for k=1:Ns
      Yc[:,:,k],Vm[k,:],Hm[k,:] = ABYZLM(ZL[:,:,k],YL[:,:,k],lenght,w[k])  
    end
    
    # Characteristic Admittance Fitting
    YcPoles,YcResidues,YcConstant,YcProportional = YcFit(Yc,f,Ns,Ncon)
    NpYc = size(YcPoles, 1)
    
    # Hk fit
    HkPoles,HkResidues,HkConstant,HkProportional,md = HkFitTrace(Hm,Vm,ZL,YL,f,Ns,lenght,Ncon)
    NpH = size(HkPoles, 2)
    Ng  = size(HkResidues, 3)
    
    # Unknown currents
    Ix = Current(zeros(Ncon))
    Iy = Current(zeros(Ncon))
    Ishx = Current(zeros(Ncon))
    Ishy = Current(zeros(Ncon))
    Iauxx = Current(zeros(Ncon))
    Iauxy = Current(zeros(Ncon))
    Irx = Current(zeros(Ncon))
    Iry = Current(zeros(Ncon))
    # Unknown model states
    Wx = Unknown(zeros(Ncon, NpYc))
    Wy = Unknown(zeros(Ncon, NpYc))
    Xx = Unknown(zeros(Ncon, NpH, Ng))
    Xy = Unknown(zeros(Ncon, NpH, Ng))
     
    # Model equations
    {
     RefBranch(Vx, Ix)
     RefBranch(Vy, Iy)
     # sum the currents
     Ix - Ishx + Iauxx
     Iy - Ishy + Iauxy
     # State-space model per Octavio et al.'s book chapter
     # eq 39:
     YcConstant * Vx + sum(Wx, 2) - Ishx
     YcConstant * Vy + sum(Wy, 2) - Ishy
     # eq 40:
     {YcPoles[i] * Wx[:,i] + YcResidues[:,:,i] * Vx  -  der(Wx[:,i])   for i in 1:NpYc}
     {YcPoles[i] * Wy[:,i] + YcResidues[:,:,i] * Vy  -  der(Wy[:,i])   for i in 1:NpYc}
     # eq 43:   
     sum(sum(Xx,3),2) - Iauxx   
     sum(sum(Xy,3),2) - Iauxy   
     # eq 44:  
     {HkPoles[k,i] * Xx[:,k,i] + HkResidues[:,:,k,i] * delay(Iry, md[k])  -  der(Xx[:,k,i])  for k in 1:Ng, i in 1:NpH}[:]
     {HkPoles[k,i] * Xy[:,k,i] + HkResidues[:,:,k,i] * delay(Irx, md[k])  -  der(Xy[:,k,i])  for k in 1:Ng, i in 1:NpH}[:]
     # find the reflected current
     Irx + Iauxx - 2 * Ishx
     Iry + Iauxy - 2 * Ishy
    }
     
end

function ex_LineEnergization()
    nv = Voltage(zeros(3), "Vs")
    ns = Voltage(zeros(3), "Vs")
    nr = Voltage(zeros(3), "Vr")
    g = 0.0
    
    # Line Geometry
    # column 1-conductor number
    # column 2-- x position of each cond in m
    # column 3-- y position of each cod in m
    # column 4-- radii of each conductor
    # column 5-- number of conductor in bundle
    # column 6-- distance between conductors in bundle
    # column 7-conductor resistivity
    # column 8-conductor relative permitivity
    Geom=[1  0 20 0.0153 3 0.4 2.826e-8 1e3
          2 10 20 0.0153 3 0.4 2.826e-8 1e3
          3 20 20 0.0153 3 0.4 2.826e-8 1e3]
    lenght = 150e3 # Line lenght
    {
     SineVoltage(nv, g, 600.0, 60.0, [0, -2/3*pi, 2/3*pi])
     Resistor(nv, ns, fill(600., 3))
     Wideband(ns, nr, Geom, lenght)
     Resistor(nr, g, fill(600., 3))
    }
end

  ## load("Winston/src/Plot")
  ## using Plot
  ## plot(vt,Vi[1,a1:a2],"-r",vt,Vi[2,a1:a2],"-g",vt,Vi[3,a1:a2],"-b",
  ##      vt,Vf[1,a1:a2],";r",vt,Vf[2,a1:a2],";g",vt,Vf[3,a1:a2],";b",) | print
       

m = ex_LineEnergization()
f = elaborate(m)
s = create_sim(f)
y = sim(s, 0.2)
