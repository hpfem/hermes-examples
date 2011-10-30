KNOWN ISSUES
============

There are a few examples that do not work properly yet. 
We will fix them as we find time for that. If you 
happen to be interested in one of these, let us know 
and we will give it a higher priority. If you are able 
to fix the example yourself, we'll buy you a beer :)

navier-stokes/rayleigh-benard
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Does not exhibit the vortices as expected. This might be
just a wrong choice of parameters. 

2d-benchmarks-nist/09-wave-front-kelly-dealii
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Runs till the end and then segfaults.
Runs fine on Windows.

1d/system
~~~~~~~~~

Segfaults.

neutronics/4-group
~~~~~~~~~~~~~~~~~~

For some reason the fourth solution component remains zero. 

neutronics/4-group-adapt
~~~~~~~~~~~~~~~~~~~~~~~~

See line #300 of main.cpp for problem description.

richards/seepage-adapt 
~~~~~~~~~~~~~~~~~~~~~~

Needs conversion to new weak forms.

richards/capillary-barrier-rk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Somehow the solution might not be right.

Richard's examples in general
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most of them do a void first time step.

heat-transfer/heat-and-moisture
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Adaptivity did not work well since the two fields had 
several orders of magnitude difference. So I rescaled 
the equations which improved adaptivity, but I must have 
introduced a bug, since the results are now different from 
before. The rescaling was just to multiply moisture ‘w’ 
with a scaling coefficient.

Also the polynomial degree and time step were raised.

maxwell/maxwell-debye
~~~~~~~~~~~~~~~~~~~~~

This example uses the R-K method with two Hcurl spaces 
and one H1 space. Now matter how small the time step size, 
after the first time the solution jumps abruptly. Some 
time ago we observed something similar with the resonator 
example. I think that this might be a bug in rk_time_step().

flame-propagation/laminar flame
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After conversion to new forms, Newton stopped converging.
Forms were checked many times. The problem may be in 
rk_time_step().
