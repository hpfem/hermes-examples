from numpy import arctan, sqrt, pi, sin, cos
t = 0.1  # thickness
l = 0.7  # length 
a = sqrt(l**2 - (l - t)**2)
print "a =", a
alpha = arctan(t/l)
print "alpha =", alpha
beta = delta - alpha
print "beta =", beta
gamma = pi/2 - 2*delta
print "gamma =", gamma
delta = arctan(a/(l-t))
print "delta =", delta
print "alpha_deg =", 180*alpha/pi
print "beta_deg =", 180*beta/pi
print "gamma_deg =", 180*gamma/pi
print "delta_deg =", 180*delta/pi
c = (l-t)*sin(alpha)
print "c =", c
d = (l-t)*cos(alpha)
print "d =", d
e = e = (l-t)*sin(delta)
print "e =", e
f = (l-t)*cos(delta)
print "f =", f
q = q = sqrt(2)/2
print "q =", q
print "l_minus_qt =", l - q * t
print "minus_qt =", - q * t

