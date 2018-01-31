import numpy
import sys

# Cosmo params
OM = 0.25648
h = 0.71903

b1 = 1.666
b2 = 1.666
nbar = 1.046e-3
z = 1
dzbin = 0.1
A = 1000         # deg^2
kmax = 0.20*h    # Mpc^-1, no h
sigmaR1 = 5.0    # Mpc, radial smearing
sigmaR2 = 5.0    # Mpc, radial smearing

# Over-ride
z = float(sys.argv[1])
dzbin = float(sys.argv[2])
b1 = b2 = float(sys.argv[3])
nbar = float(sys.argv[4])
A = float(sys.argv[5])

# 2 tracer
db=0.0
b1 -= db; b2 += db

# finger of God lengths: dz=(1+z)*0.001
dDdz = 2997.92458/h/numpy.sqrt(1-OM+OM*(1+z)**3)
sigmaR1 = sigmaR2 = (1+z)*float(sys.argv[6])*dDdz

# Turn off FoG
#sigmaR1=sigmaR2=0.01

# Function to get Fisher per mode for alpha, beta, Pv, log s^2_1, log s^2_2
# nP for each pair at mu=0
def getFmat(b,alpha,beta,ks1,ks2,nP1,nP2):
  # Number of mu bins
  nmu = 25
  # Initialize
  # pars = (alpha, beta, ln Pv, ln s_1^2, ln s_2^2)
  F = numpy.zeros((5,5))

  # Normalize to Pv = 1 so
  # n1 = nP1*beta^2
  # n2 = nP2*(beta/alpha)^2

  for imu in range(nmu):
    mu = (imu+.5)/nmu
    T1 = numpy.exp(ks1**2*mu**2)
    T2 = numpy.exp(ks2**2*mu**2)
    # Cov = [ (1/beta+mu^2)^2 Pv + 1/n1          (1/beta+mu^2)(alpha/beta+mu^2) Pv ]
    #       [ (1/beta+mu^2)(alpha/beta+mu^2) Pv  (alpha/beta+mu^2)^2 Pv + 1/n2     ]
    C = numpy.zeros((2,2))
    C[0][0] = (1/beta+mu**2)**2 + 1./nP1/beta**2*T1
    C[1][1] = (alpha/beta+mu**2)**2 + alpha**2/nP2/beta**2*T2
    C[0][1] = C[1][0] = (1/beta+mu**2)*(alpha/beta+mu**2)
    Cinv = numpy.linalg.inv(C)

    # Now parameter derivatives
    Cp = numpy.zeros((5,2,2))
    # alpha-derivatives
    Cp[0][0][1] = Cp[0][1][0] = (1/beta+mu**2)/beta
    Cp[0][1][1] = 2*(alpha/beta+mu**2)/beta
    # beta-derivatives
    Cp[1][0][0] = 2*(1/beta+mu**2)*(-1./beta**2)
    Cp[1][0][1] = Cp[1][1][0] = (-1./beta**2)*(alpha/beta+mu**2) + (1./beta+mu**2)*(-alpha/beta**2)
    Cp[1][1][1] = 2*(alpha/beta+mu**2)*(-alpha/beta**2)
    # Pv-derivatives
    Cp[2][0][0] = (1/beta+mu**2)**2
    Cp[2][0][1] = Cp[2][1][0] = (1/beta+mu**2)*(alpha/beta+mu**2)
    Cp[2][1][1] = (alpha/beta+mu**2)**2
    # log ks derivatives
    Cp[3][0][0] = (1/beta+mu**2)**2*ks1**2*mu**2
    Cp[3][0][1] = Cp[3][1][0] = (1/beta+mu**2)*(alpha/beta+mu**2)*ks1**2*mu**2/2.
    Cp[4][0][1] = Cp[4][1][0] = (1/beta+mu**2)*(alpha/beta+mu**2)*ks2**2*mu**2/2.
    Cp[4][1][1] = (alpha/beta+mu**2)**2*ks2**2*mu**2

    # Build Fisher
    for j in range(5):
      for k in range(5):
        F[j,k] += numpy.matrix.trace(numpy.dot(numpy.dot(Cp[j,:,:],Cinv), numpy.dot(Cp[k,:,:],Cinv)))/2.
  F /= nmu

  return F


# --- main program here ---

# Get power spectrum
k_, Pk_ = numpy.loadtxt("wmap5Pk00_pk.dat", comments="#", unpack=True)
k_ *= h
Pk_ /= h**3    # Convert to just Mpc

# Get cosmological properties back then
OMz = OM/(OM + (1.-OM)/(1+z)**3)
f_r = OMz**.55
Gz = 2.5*OMz/(OMz**(4./7)-1+OMz+(1+0.5*OMz*(1+(1-OMz)/70)))
Gz /= 2.5*OM/(OM**(4./7)-1+OM+(1+0.5*OM*(1+(1-OM)/70)))
Gz /= 1+z

#print OMz, f_r, Gz
#print k_.shape
#print Pk_.shape
#print k_
#print Pk_

# Get survey volume
dDdz = 2997.92458/h/numpy.sqrt(1-OM+OM*(1+z)**3)
Dist = 2997.92458/h*z*( 0.1739275/numpy.sqrt(1-OM+OM*(1+z*0.069432)**3)
                       +0.3260725/numpy.sqrt(1-OM+OM*(1+z*0.3300095)**3)
                       +0.3260725/numpy.sqrt(1-OM+OM*(1+z*0.6699905)**3)
                       +0.1739275/numpy.sqrt(1-OM+OM*(1+z*0.930568)**3))
V = A/(180./numpy.pi)**2*dDdz*Dist**2*dzbin

# Galaxy properties
bave = (b1+b2)/2.
b = b1
alpha = b2/b1
beta = f_r/b

nk = 100
F = numpy.zeros((5,5))
for i in range(nk):
  k = kmax*(i+.5)/nk
  dk = kmax/nk
  Pk = numpy.interp(k, k_, Pk_)
  nP = nbar*Pk*bave**2*Gz**2
  nmodes = k**2*dk/(2.*numpy.pi**2)*V
  #print k, Pk, nP, nmodes
  thisF = getFmat(b,alpha,beta,k*sigmaR1,k*sigmaR2,nP/2.*(b1/bave)**2,nP/2.*(b2/bave)**2)
  #print thisF
  F += nmodes*thisF

# Fix FoG length pars?
#F[4,4] += 1e24
#F[3,3] += 1e24

#print F
CovPar = numpy.linalg.inv(F)
#print 'Covariances ='
#print CovPar
#for i in range(5):
#  print 'corr. with par',i,':', CovPar[2,i]/numpy.sqrt(CovPar[2,2]*CovPar[i,i])

#print 'sigma(fs8)/fs8 =',
print '{:7.5f}'.format(numpy.sqrt(CovPar[2,2])/2.)
