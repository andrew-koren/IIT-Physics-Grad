from pylab import *

# input parameters
p0=10000. # reference momentum [MeV/c]
m=938.272 # proton rest mass in [MeV]
c=299792458. # speed of light [m/s]
v0=p0*c/sqrt(m**2+p0**2) # reference velocity [m/s]
eps_xn=100. # normalized emittance [mm] (rms)
eps_yn=100.
beta_x=32100 # transverse amplitude function [mm] (take from MAD-X!)
beta_y=6150 # transverse amplitude function [mm] (take from MAD-X!)
alpha_x=0. # transverse alpha (I assume we start at the waist, e.g. in the middle of a quad)
alpha_y=0. # transverse alpha (I assume we start at the waist, e.g. in the middle of a quad)

# geometric emittance
eps_x=eps_xn*m/p0
eps_y=eps_yn*m/p0

# the third Twiss parameter
gamma_x=(1+alpha_x**2)/beta_x # transverse gamma
gamma_y=(1+alpha_y**2)/beta_y # transverse gamma

# calculate the ellipse semi-axes (a,b)
a_x=sqrt(2*eps_x/(gamma_x+beta_x-sqrt((gamma_x+beta_x)**2-4))) # major semi-axis
b_x=sqrt(2*eps_x/(gamma_x+beta_x+sqrt((gamma_x+beta_x)**2-4))) # minor semi-axis

a_y=sqrt(2*eps_y/(gamma_y+beta_y-sqrt((gamma_y+beta_y)**2-4))) # major semi-axis
b_y=sqrt(2*eps_y/(gamma_y+beta_y+sqrt((gamma_y+beta_y)**2-4))) # minor semi-axis

# calculate the orientation of the beam ellipse
# this is the trickiest part, and some degenerate cases need to be handled carefully
theta_x=0
if (gamma_x>beta_x):
    theta_x=pi/2
if (alpha_x!=0):
    theta_x=theta_x+1/2*atan(2*alpha_x/(gamma_x-beta_x))

theta_y=0
if (gamma_y>beta_y):
    theta_y=pi/2
if (alpha_y!=0):
    theta_y=theta_y+1/2*atan(2*alpha_y/(gamma_y-beta_y))

# calculate distribution sigmas from the beam ellipse parameters
# I call them "max" here to avoid confusion if the emiitance we start with
# is not the rms emittance
max_x=sqrt(eps_x*beta_x) # max deviation in x [mm]
max_x1=sqrt(eps_x*gamma_x) # max deviation in x'

max_y=sqrt(eps_y*beta_y) # max deviation in y [mm]
max_y1=sqrt(eps_y*gamma_y) # max deviation in y'

# This is a simplification for our "transverse only" case
# For longitudinal motion I had a similar calculation,
# removed it for simplicity
max_z=0 # max deviation in z [mm]
max_z1=0 # max deviation in z'
max_t=max_z/1e3/v0*1e9 # max t

max_px=max_x1*p0 # maximum horizontal momentum deviation [MeV/c]
max_py=max_y1*p0 # maximum vertical momentum deviation [MeV/c]
max_pz=max_z1*p0 # maximum longitudinal momentum deviation [MeV/c]

N=10000

# generate an upright ellipse with major semi-axes a_{x,y} and b_{x,y}
x=randn(N)*a_x
x1=randn(N)*b_x

y=randn(N)*a_y
y1=randn(N)*b_y

# rotate this ellipse by theta, thus obtaining the desired distribution
# (matching)
RX=[[cos(theta_x),-sin(theta_x)],[sin(theta_x),cos(theta_x)]]
RY=[[cos(theta_y),-sin(theta_y)],[sin(theta_y),cos(theta_y)]]

[x,x1]=dot(RX,[x,x1])
[y,y1]=dot(RY,[y,y1])
px=x1*p0 # switching from (x,x') to (x,px)
py=y1*p0

pz=sqrt(p0**2-px**2-py**2)


#%%
# This is a test to verify the generated distribution is consistent 
f=(gamma_x*x**2+2*alpha_x*x*x1+beta_x*x1**2<eps_x) # epsT*6 for 95% emittance
print(sum(f)/N) # check that the percentage of the beam is correct

# ellipse vs beam (transverse)
figure(0)
plot(x,x1,'.') # particle distribution
X = arange(-max_x,max_x+2*max_x/100,2*max_x/100) # horizontal coordinate
Y1 = -alpha_x/beta_x*X - 1/beta_x*sqrt((alpha_x**2-beta_x*gamma_x)*X**2+beta_x*eps_x) # lower branch of the ellipse
Y2 = -alpha_x/beta_x*X + 1/beta_x*sqrt((alpha_x**2-beta_x*gamma_x)*X**2+beta_x*eps_x) # Upper branch of the ellipse
plot(X,Y1,'k-',linewidth=3)
plot(X,Y2,'k-',linewidth=3)
grid('on')

#%%

# saving data for g4beamline
n=x.shape[0]
d=zeros((n,12))
d[:,0]=x
d[:,1]=y
d[:,2]=0
d[:,3]=px
d[:,4]=py
d[:,5]=pz
d[:,6]=0
d[:,7]=2212 # proton!
d[:,8]=arange(1,n+1)
d[:,9]=1 # track number
d[:,10]=0 # no parent
d[:,11]=1.0 # weight

dest=open('distrib.dat.test','w')
dest.write('# BLTrackFile output generated from for009\n')
dest.write('#x y z Px Py Pz t PDGid EventID TrackID ParentID Weight\n')
dest.write('#mm mm mm MeV/c MeV/c MeV/c ns - - - - -\n')

for line in d:
  dest.write('%g %g %g %g %g %g %g %i %i %i %i %g\n'%tuple(line))

dest.close()
