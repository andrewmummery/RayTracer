# README #
## RayTracer - A Python and C++ package for solving photon geodesics in full General Relativity. ##
## What is in this README? ##
* Description of what the package RayTracer.py can do
* Description of the physics problem that is solved
* Brief description of the algorithm used
* Instructions of how to set-up and use this program

## What do I do?  ##
![](https://github.com/andrewmummery/RayTracer/blob/main/ExampleAnimation.gif)

![](https://github.com/andrewmummery/RayTracer/blob/main/ExampleNewtonianFigure.png)
![](https://github.com/andrewmummery/RayTracer/blob/main/ExampleGRFigure.png)

## The Physics ##
In brief, we solve the standard equations of motion for a photon moving through a general relativistic spacetime:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{{\rm&space;d^2}x^\mu}{{\rm&space;d}\tau^2}&space;&plus;&space;\Gamma^\mu_{\sigma&space;\nu}&space;\frac{{\rm&space;d}x^\sigma}{{\rm&space;d}\tau}\frac{{\rm&space;d}x^\nu}{{\rm&space;d}\tau}&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{{\rm&space;d^2}x^\mu}{{\rm&space;d}\tau^2}&space;&plus;&space;\Gamma^\mu_{\sigma&space;\nu}&space;\frac{{\rm&space;d}x^\sigma}{{\rm&space;d}\tau}\frac{{\rm&space;d}x^\nu}{{\rm&space;d}\tau}&space;=&space;0" title="\frac{{\rm d^2}x^\mu}{{\rm d}\tau^2} + \Gamma^\mu_{\sigma \nu} \frac{{\rm d}x^\sigma}{{\rm d}\tau}\frac{{\rm d}x^\nu}{{\rm d}\tau} = 0" /></a>

where our co-ordinates are <a href="https://www.codecogs.com/eqnedit.php?latex=x^\mu&space;=&space;(t,&space;r,&space;\theta,&space;\varphi)," target="_blank"><img src="https://latex.codecogs.com/gif.latex?x^\mu&space;=&space;(t,&space;r,&space;\theta,&space;\varphi)," title="x^\mu = (t, r, \theta, \varphi)," /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\Gamma^\mu_{\sigma&space;\nu}&space;\equiv&space;\frac{1}{2}&space;g^{\mu&space;\alpha}&space;\left[\frac{\partial&space;g_{\alpha&space;\nu}&space;}{\partial&space;x^\sigma}&space;&plus;&space;\frac{\partial&space;g_{\alpha&space;\sigma}&space;}{\partial&space;x^\nu}&space;-&space;\frac{\partial&space;g_{\sigma&space;\nu}&space;}{\partial&space;x^\alpha}\right&space;]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Gamma^\mu_{\sigma&space;\nu}&space;\equiv&space;\frac{1}{2}&space;g^{\mu&space;\alpha}&space;\left[\frac{\partial&space;g_{\alpha&space;\nu}&space;}{\partial&space;x^\sigma}&space;&plus;&space;\frac{\partial&space;g_{\alpha&space;\sigma}&space;}{\partial&space;x^\nu}&space;-&space;\frac{\partial&space;g_{\sigma&space;\nu}&space;}{\partial&space;x^\alpha}\right&space;]" title="\Gamma^\mu_{\sigma \nu} \equiv \frac{1}{2} g^{\mu \alpha} \left[\frac{\partial g_{\alpha \nu} }{\partial x^\sigma} + \frac{\partial g_{\alpha \sigma} }{\partial x^\nu} - \frac{\partial g_{\sigma \nu} }{\partial x^\alpha}\right ]" /></a>

and <a href="https://www.codecogs.com/eqnedit.php?latex=g_{\mu\nu}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?g_{\mu\nu}" title="g_{\mu\nu}" /></a> is the Kerr metric.

The idea is to start with a photon in the observers image plane, specified by two (cartesian) co-ordinates <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha" title="\alpha" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\beta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /></a>, at a large distance <a href="https://www.codecogs.com/eqnedit.php?latex=D" target="_blank"><img src="https://latex.codecogs.com/gif.latex?D" title="D" /></a> away from a blackhole of (dimensionless) spin <a href="https://www.codecogs.com/eqnedit.php?latex=a" target="_blank"><img src="https://latex.codecogs.com/gif.latex?a" title="a" /></a>, observed at an inclination angle <a href="https://www.codecogs.com/eqnedit.php?latex=\theta_{\rm&space;obs}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_{\rm&space;obs}" title="\theta_{\rm obs}" /></a>. The co-ordinate settup is displayed in the below figure:  

![](https://github.com/andrewmummery/RayTracer/blob/main/Schematic_of_co-ordinates.png)

Fortunately, in the Kerr metric we have two constants of motion, the angular momentum and energy of the photon. This of course results from the time and phi independence of the Kerr metric. The two constants of motion can be calculated at the initial conditions in the image plane. By defining:

<a href="https://www.codecogs.com/eqnedit.php?latex=L&space;\equiv&space;p_{\varphi}&space;=&space;g_{\varphi\varphi}&space;\frac{{\rm&space;d}\varphi}{{\rm&space;d}\tau}&space;&plus;&space;g_{\varphi&space;t}\frac{{\rm&space;d}t}{{\rm&space;d}\tau}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?L&space;\equiv&space;p_{\varphi}&space;=&space;g_{\varphi\varphi}&space;\frac{{\rm&space;d}\varphi}{{\rm&space;d}\tau}&space;&plus;&space;g_{\varphi&space;t}\frac{{\rm&space;d}t}{{\rm&space;d}\tau}" title="L \equiv p_\varphi = g_{\varphi\varphi} \frac{{\rm d}\varphi}{{\rm d}\tau} + g_{\varphi t}\frac{{\rm d}t}{{\rm d}\tau}" /></a>

and 

<a href="https://www.codecogs.com/eqnedit.php?latex=E&space;\equiv&space;-p_t&space;=&space;-&space;g_{tt}&space;\frac{{\rm&space;d}t}{{\rm&space;d}\tau}&space;-&space;g_{t&space;\varphi}\frac{{\rm&space;d}\varphi}{{\rm&space;d}\tau}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?E&space;\equiv&space;-p_t&space;=&space;-&space;g_{tt}&space;\frac{{\rm&space;d}t}{{\rm&space;d}\tau}&space;-&space;g_{t&space;\varphi}\frac{{\rm&space;d}\varphi}{{\rm&space;d}\tau}" title="E \equiv -p_t = - g_{tt} \frac{{\rm d}t}{{\rm d}\tau} - g_{t \varphi}\frac{{\rm d}\varphi}{{\rm d}\tau}" /></a>

We can write two simplified equations of motion for the <a href="https://www.codecogs.com/eqnedit.php?latex=\varphi" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\varphi" title="\varphi" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=t" target="_blank"><img src="https://latex.codecogs.com/gif.latex?t" title="t" /></a> co-ordinates:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{{\rm&space;d}&space;t&space;}{{\rm&space;d}\tau'}&space;=&space;-&space;\frac{l&space;g_{&space;\varphi&space;t}&space;&plus;&space;g_{&space;\varphi&space;\varphi}}{g_{\varphi&space;\varphi}&space;g_{tt}&space;-&space;g_{\varphi&space;t}^2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{{\rm&space;d}&space;t&space;}{{\rm&space;d}\tau'}&space;=&space;-&space;\frac{l&space;g_{&space;\varphi&space;t}&space;&plus;&space;g_{&space;\varphi&space;\varphi}}{g_{\varphi&space;\varphi}&space;g_{tt}&space;-&space;g_{\varphi&space;t}^2}" title="\frac{{\rm d} t }{{\rm d}\tau'} = - \frac{l g_{ \varphi t} + g_{ \varphi \varphi}}{g_{\varphi \varphi} g_{tt} - g_{\varphi t}^2}" /></a>

and 

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{{\rm&space;d}&space;\varphi&space;}{{\rm&space;d}\tau'}&space;=&space;\frac{l&space;g_{&space;t&space;t}&space;&plus;&space;g_{&space;t&space;\varphi}}{g_{\varphi&space;\varphi}&space;g_{tt}&space;-&space;g_{\varphi&space;t}^2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{{\rm&space;d}&space;\varphi&space;}{{\rm&space;d}\tau'}&space;=&space;\frac{l&space;g_{&space;t&space;t}&space;&plus;&space;g_{&space;t&space;\varphi}}{g_{\varphi&space;\varphi}&space;g_{tt}&space;-&space;g_{\varphi&space;t}^2}" title="\frac{{\rm d} \varphi }{{\rm d}\tau'} = \frac{l g_{ t t} + g_{ t \varphi}}{g_{\varphi \varphi} g_{tt} - g_{\varphi t}^2}" /></a>

where we normalise our proper time and angular momentum values by the photons energy: <a href="https://www.codecogs.com/eqnedit.php?latex=\tau'&space;=&space;E&space;\tau,&space;\,\,\,&space;l&space;=&space;L/E" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\tau'&space;=&space;E&space;\tau,&space;\,\,\,&space;l&space;=&space;L/E" title="\tau' = E \tau, \,\,\, l = L/E" /></a>. 

Our final two equations of motion are then:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{{\rm&space;d^2}&space;r}{{\rm&space;d}{\tau'}^2}&space;=&space;-&space;\Gamma^r_{\sigma&space;\nu}&space;\frac{{\rm&space;d}x^\sigma}{{\rm&space;d}\tau'}\frac{{\rm&space;d}x^\nu}{{\rm&space;d}\tau'}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{{\rm&space;d^2}&space;r}{{\rm&space;d}{\tau'}^2}&space;=&space;-&space;\Gamma^r_{\sigma&space;\nu}&space;\frac{{\rm&space;d}x^\sigma}{{\rm&space;d}\tau'}\frac{{\rm&space;d}x^\nu}{{\rm&space;d}\tau'}" title="\frac{{\rm d^2} r}{{\rm d}{\tau'}^2} = - \Gamma^r_{\sigma \nu} \frac{{\rm d}x^\sigma}{{\rm d}\tau'}\frac{{\rm d}x^\nu}{{\rm d}\tau'}" /></a>

and 

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{{\rm&space;d^2}&space;\theta&space;}{{\rm&space;d}{\tau'}^2}&space;=&space;-&space;\Gamma^\theta_{\sigma&space;\nu}&space;\frac{{\rm&space;d}x^\sigma}{{\rm&space;d}\tau'}\frac{{\rm&space;d}x^\nu}{{\rm&space;d}\tau'}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{{\rm&space;d^2}&space;\theta&space;}{{\rm&space;d}{\tau'}^2}&space;=&space;-&space;\Gamma^\theta_{\sigma&space;\nu}&space;\frac{{\rm&space;d}x^\sigma}{{\rm&space;d}\tau'}\frac{{\rm&space;d}x^\nu}{{\rm&space;d}\tau'}" title="\frac{{\rm d^2} \theta }{{\rm d}{\tau'}^2} = - \Gamma^\theta_{\sigma \nu} \frac{{\rm d}x^\sigma}{{\rm d}\tau'}\frac{{\rm d}x^\nu}{{\rm d}\tau'}" /></a>



## Brief description of the algorithm ##
The actual numerical solution of the photon equations is done in the file ray_tracing.cpp. The algorithm works in the following way.

#### Step 1 ####
Read in the parameters <a href="https://www.codecogs.com/eqnedit.php?latex=(\alpha,&space;\beta,&space;\theta_{\rm&space;obs},&space;a)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(\alpha,&space;\beta,&space;\theta_{\rm&space;obs},&space;a)" title="(\alpha, \beta, \theta_{\rm obs}, a)" /></a>  and calculate the initial values (image plane) of the photons 4-position and 4-momentum. 

#### Step 2 ####
Evolve the photons trajectory back towards the blackhole-disc plane by solving the four dynamical equations described above. (We actually split the dynamical equations into six coupled first order equations which are solved by the RK4 method, for more details see Appendix A of https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.5655M/abstract)

#### Step 3 ####
The path of the photon is terminated when it either: 
* Hits the event horizon of the blackhole
* Hits an accretion disc in the equatorial plane of the blackhole (the presence of the disc can be switched on and off)
* Passes between the disc inner edge and the blackholes event horizon and escapes to a large distance (~ 30 gravitational radii) from the blackhole

#### Step 4 ####
Load the trajectories <a href="https://www.codecogs.com/eqnedit.php?latex=r(\tau'),&space; &space;\theta(\tau'),&space; &space;\varphi&space;(\tau'),&space;&space;t(\tau')" target="_blank"><img src="https://latex.codecogs.com/gif.latex?r(\tau'),&space;&space;\theta(\tau'),&space;&space;\varphi&space;(\tau'),&space;&space;t(\tau')" title="r(\tau'), \, \theta(\tau'), \, \varphi (\tau'), \, t(\tau')" /></a> into Python, where they can be analysed. These trajectories are then converted to cartesian co-ordinates <a href="https://www.codecogs.com/eqnedit.php?latex=x(\tau'),&space;y(\tau'),&space;z(\tau')" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x(\tau'),&space;y(\tau'),&space;z(\tau')" title="x(\tau'), y(\tau'), z(\tau')" /></a> internally in the Python script, and are then plotted. 



## How do I get set up? ##

* Pull the codes into their own folder
* Run examples.py for a series of examples of what the code package can do. Note: you will have to manually edit examples.py to display different examples, the script defaults to 'example_zero' 
* Make your own images/animations by creating Python scripts in the same folder as RayTracer.py, having copied the first 20 lines of examples.py

