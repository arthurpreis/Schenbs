
using PyPlot, Statistics, LinearAlgebra


kB = 1.38064852e-23        # m^2kgs^(-2)K^(-1). Boltzman constant
hbar = 6.62607004e-34/2π;  # m^2kg/s. Reduced Plank's constant


R = 0.325                                              # m. Sphere radius
MS = 1150                                              # kg. Sphere mass
V = 4π/3*R^3                                           # m^3. Sphere volume
ρ = 8065.68                                            # kg/m^3. Sphere density
f0n =[3172.5, 3183.0, 3213.6, 3222.9, 3240.0]          # measured natural frequencies
ω0n = 2π*f0n
f0 = mean(f0n)                                         # mean of natural frequencies, transducers frequencies
ω0 = 2π*f0                                             # angular frequencies
#ω0n = ω0*ones(5) # teste as frequências degeneradas retirar depois
Q = 1.23e6                                             # mechanical Q
T = 4.2                                                # Sphere temperature 
α = 2.86239                                            # radial component α(r) in r=R 
χ = 0.60138            
Meff = 5/6*(χ/2)*V*ρ;                                  # kg. effective mass of the sphere
β = ω0/2Q 


MR2 = 0.0000123                                       # kg. second ressonator mass
MR1 = sqrt(Meff*MR2)                                  # kg. first ressonator mass                          
μ = sqrt(MR1/Meff)                                    # mass ratio
ν = sqrt(MS/Meff)
ωR1 = ω0
ωR2 = ω0
QR1 = 10^6
QR2 = 10^5;


nR = 6                          # number of transducers
nm = 5                          # number of sphere modes
kS = MS*ω0n.^2                   # kg/s^2. sphere spring constant
kR1 = MR1*ωR1^2                 # kg/s^2. first ressonator spring constant
kR2 = MR2*ωR2^2                 # kg/s^2. second ressonator spring constant
HS=MS*ω0/Q                      # coeficiente de amortecimento da esfera em kg m/s
HR1=MR1*ωR1/QR1                 # coeficiente de amortecimento do primeiro ressonador em kg m/s
HR2=MR2*ωR2/QR2;                # coeficiente de amortecimento do segundo ressonador em kg m/s


phi = (1+√5)/2   # golden ratio
x = [-(sqrt(3)*phi*sqrt(5*phi^2+4*phi+1)-2*phi-1)/(2*sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1)),
     -(phi^2+phi)/(sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1)),
     (sqrt(3)*phi*sqrt(5*phi^2+4*phi+1)+2*phi+1)/(2*sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1)),
     (sqrt(3)*sqrt(5*phi^2+4*phi+1)-phi^2)/(2*sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1)),
     -(sqrt(3)*sqrt(5*phi^2+4*phi+1)+phi^2)/(2*sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1)),
     phi^2/(sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1))]
y = [(phi*sqrt(5*phi^2+4*phi+1)+2*sqrt(3)*phi+sqrt(3))/(2*sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1)),
     -(sqrt(3)*phi^2+sqrt(3)*phi)/(sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1)),
     -(phi*sqrt(5*phi^2+4*phi+1)-2*sqrt(3)*phi-sqrt(3))/(2*sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1)),
     -(sqrt(5*phi^2+4*phi+1)+sqrt(3)*phi^2)/(2*sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1)),
     (sqrt(5*phi^2+4*phi+1)-sqrt(3)*phi^2)/(2*sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1)),
     (sqrt(3)*phi^2)/(sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1))]
z = [phi/(sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1)),
     -(phi^2-2*phi-1)/(sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1)),
       phi/(sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1)),
      (2*phi^2+phi)/(sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1)),
      (2*phi^2+phi)/(sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1)),
      (phi^2+2*phi+1)/(sqrt(phi^2+1)*sqrt(5*phi^2+4*phi+1))]


ϕ = atan.(y,x)*180/π


θ = acos.(z)*180/π


B = zeros(5,6)
B[1,:] = sqrt(15/16π)*(x.^2-y.^2)
B[2,:] = sqrt(15/16π)*(2*x.*y)
B[3,:] = sqrt(15/16π)*(2*y.*z)
B[4,:] = sqrt(15/16π)*(2*x.*z)
B[5,:] = sqrt(15/16π)*(3*z.^2 .- 1)/(sqrt(3));


BT = B'
B*BT/(3/2π)


F3T = 3206.3                                  # frequência do transdutor 
w3T = 2π*F3T                                  # frequência angular do transdutor
w3Tsq = w3T^2                                 # quadrado da frequência angular do transdutor

FO1 = 3172.5
FO2 = 3240.0                                  # f- e f+ do transdutor
wO1 = 2π*FO1        
wO2 = 2π*FO2

Fpump = 10e9                                  # frequência da bomba
wpump = 2π*Fpump                              # frequência angular da bomba
wpumpsq = wpump^2                             # quadrado da frequência angular da bomba
Betae = 0.3                                   # fator de acoplamento elétrico
dfdx = 7.26e14                                # shift de frequência com a distância
Pinc = 5e-9                                   # potência incidente do oscilador
Tamp = 8                                      # temperatura de ruído do amplificador
Qe = 3e5                                      # electrical quality factor
Spm = 10^(-13.0)                              # densidade espectral do ruído de fase
Sam = 10^(-14.0)                              # densidade espectral do ruído de amplitude
Lamp = 0.5
Sseries = Lamp*(Tamp*kB/Pinc)*(Fpump/2/Qe/dfdx)^2  # densidade espectral do ruído de "series"

Qfac = (Qe^2/(1+4*Qe^2*ω0^2/wpump^2))^0.5
beta0 = (4*Betae*Pinc)/((1+Betae)*MR2*ω0^2*wpump)*((2*Qfac*dfdx)/(Fpump))^2
w3 = w3T*(1-(3*3^0.5*beta0)/(16))^0.5

Dw = 50000*2π
DUO1 = (Dw+wO1)/wpump  
DLO1 = (Dw-wO1)/wpump
D0 = Dw/wpump

ab1 = w3^2/wO1
ac1 = wO1^2/wpump^2 
ad1 = (1+4*Qe^2*ac1)
bb1 = (beta0*ab1*ad1)
bb3 = 2*(1+4*Qe^2*D0^2)
bbO1 = bb1/bb3  
bcU1 = (2*Qe*DUO1)^2+1
bcL1 = (2*Qe*DLO1)^2+1 
bcUO1 = 1/bcU1
bcLO1 = 1/bcL1  

ZmR = -1*bbO1*(bcUO1-bcLO1)
ZmI = -im*(bbO1*2*Qe*DUO1*(bcUO1-bcLO1))

Sbaction = Pinc^2/2/wpump^2*(2*Qe*dfdx/Fpump)^2*Sam;       # double sided  


eye(n) = one(zeros(n,n)) 


N = [eye(nm)/ν        zeros(nm, nR)       zeros(nm, nR);
     zeros(nR,nm)     eye(nR)/μ           zeros(nR,nR);
     zeros(nR,nm)     zeros(nR,nR)        eye(nR)/μ^2]
iN = N^-1;


M = [ν^2*eye(nm) zeros(nm,nR) zeros(nm,nR);
     μ^2*α*BT    μ^2*eye(nR)  zeros(nR, nR);
     μ^4*α*BT    μ^4*eye(nR)  μ^4*eye(nR)]
iM = M^-1;


K = [ν^2*diagm(ω0n.^2 ./ ω0.^2)   -μ^2*α*B                 zeros(nm, nR);
     zeros(nR,nm)                  μ^2*eye(nR)            -μ^4*eye(nR);
     zeros(nR,nm)                  zeros(nR,nR)            μ^4*eye(nR)];


H = K;


P = [eye(nm)      -α*B           zeros(nm,nR);
     zeros(nR,nm)  eye(nR)      -eye(nR);
     zeros(nR,nm)  zeros(nR,nR)  eye(nR)]
iP = P^-1;


My = N*M*N
Ky = N*K*N
Py = N*P;


Kz = My^-1*Ky
Pz = My^-1*Py;


sum(abs.(Kz-Kz'))


for n in 1:5
    Kz[n,n] = ω0n[n]^2 / ω0^2
end


Kz = (Kz+Kz')/2      # eigen is very sensitive to very small asymmetries       
ω2, U = eigen(Kz)
f17 =ω0*sqrt.(ω2)/2π


D = U'*Kz*U
# Inverse of J
function iJ(w)
    iJ = zeros(17,17)*im
    for k in 1:17
        iJ[k,k] = -Meff*w^2+2im*β*Meff*ω0*D[k,k]+Meff*ω0^2*D[k,k]
    end
    iJ
end
function J(w)
    J = zeros(17,17)*im
    for k in 1:17
        J[k,k] = 1/(-Meff*w^2+2im*β*Meff*ω0*D[k,k]+Meff*ω0^2*D[k,k])
    end
    J
end


D


L(w) = N*U*J(w)*U'*iN*iM*P
iL(w) = iP*M*N*U*iJ(w)*U'*iN       # inverse of L
Y(w) = im*w*L(w)
Sxx(w) = 4kB*T/w^2*real(Y(w))
iSxx(w) = Sxx(w)^-1
L31(w) = L(w)[12:17,1:5]
R33(w) = iSxx(w)[12:17,12:17]
hS(w) = sqrt(1/abs(tr(L31(w)'*R33(w)*L31(w))))/(1/2*MS*χ*R*w^2)


af = 3100:0.1:3300
ahs = zeros(length(af))
for (i,f) in enumerate(af)
    ahs[i] = hS(2π*f)
end
figure(1)
    semilogy(af,ahs)
    title("Schenberg sensitivity")
    xlabel(L"$f$ [Hz]")
    ylabel(L"$h_S$ [Hz]$^{-1/2}$")
    grid(true, which="minor")
    grid(true, which="both")

