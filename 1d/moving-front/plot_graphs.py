# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
pylab.title("DOF history")
pylab.xlabel("Physical time")
pylab.ylabel("Number of degrees of freedom")
data = numpy.loadtxt("dof_history_hp.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, "-s", label="hp-FEM")
data = numpy.loadtxt("dof_history_h2.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, "-s", label="h-FEM (p=2)")
data = numpy.loadtxt("dof_history_h1.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, "-s", label="h-FEM (p=1)")

legend()

# finalize
show()
