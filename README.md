# README #
## RayTracer - A Python and C++ package for solving photon geodesics in full General Relativity. ##
### What is in this README? ###
* Description of what the package RayTracer.py can do
* Description of the physics problem that is solved
* Brief description of the algorithm used
* Instructions of how to set-up and use this program

### What do I do?  ###
![](https://github.com/andrewmummery/RayTracer/blob/main/ExampleAnimation.gif)

![](https://github.com/andrewmummery/RayTracer/blob/main/ExampleNewtonianFigure.png)
![](https://github.com/andrewmummery/RayTracer/blob/main/ExampleGRFigure.png)

### The Physics ###
In brief, we solve the standard equations of motion for a photon moving through a general relativistic spacetime:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{{\rm&space;d^2}x^\mu}{{\rm&space;d}\tau^2}&space;&plus;&space;\Gamma^\mu_{\sigma&space;\nu}&space;\frac{{\rm&space;d}x^\sigma}{{\rm&space;d}\tau}\frac{{\rm&space;d}x^\nu}{{\rm&space;d}\tau}&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{{\rm&space;d^2}x^\mu}{{\rm&space;d}\tau^2}&space;&plus;&space;\Gamma^\mu_{\sigma&space;\nu}&space;\frac{{\rm&space;d}x^\sigma}{{\rm&space;d}\tau}\frac{{\rm&space;d}x^\nu}{{\rm&space;d}\tau}&space;=&space;0" title="\frac{{\rm d^2}x^\mu}{{\rm d}\tau^2} + \Gamma^\mu_{\sigma \nu} \frac{{\rm d}x^\sigma}{{\rm d}\tau}\frac{{\rm d}x^\nu}{{\rm d}\tau} = 0" /></a>

where our co-ordinates are <a href="https://www.codecogs.com/eqnedit.php?latex=x^\mu&space;=&space;(t,&space;r,&space;\theta,&space;\varphi)," target="_blank"><img src="https://latex.codecogs.com/gif.latex?x^\mu&space;=&space;(t,&space;r,&space;\theta,&space;\varphi)," title="x^\mu = (t, r, \theta, \varphi)," /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\Gamma^\mu_{\sigma&space;\nu}&space;\equiv&space;\frac{1}{2}&space;g^{\mu&space;\alpha}&space;\left[\frac{\partial&space;g_{\alpha&space;\nu}&space;}{\partial&space;x^\sigma}&space;&plus;&space;\frac{\partial&space;g_{\alpha&space;\sigma}&space;}{\partial&space;x^\nu}&space;-&space;\frac{\partial&space;g_{\sigma&space;\nu}&space;}{\partial&space;x^\alpha}\right&space;]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Gamma^\mu_{\sigma&space;\nu}&space;\equiv&space;\frac{1}{2}&space;g^{\mu&space;\alpha}&space;\left[\frac{\partial&space;g_{\alpha&space;\nu}&space;}{\partial&space;x^\sigma}&space;&plus;&space;\frac{\partial&space;g_{\alpha&space;\sigma}&space;}{\partial&space;x^\nu}&space;-&space;\frac{\partial&space;g_{\sigma&space;\nu}&space;}{\partial&space;x^\alpha}\right&space;]" title="\Gamma^\mu_{\sigma \nu} \equiv \frac{1}{2} g^{\mu \alpha} \left[\frac{\partial g_{\alpha \nu} }{\partial x^\sigma} + \frac{\partial g_{\alpha \sigma} }{\partial x^\nu} - \frac{\partial g_{\sigma \nu} }{\partial x^\alpha}\right ]" /></a>

and <a href="https://www.codecogs.com/eqnedit.php?latex=g_{\mu\nu}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?g_{\mu\nu}" title="g_{\mu\nu}" /></a> is the Kerr metric.

The idea is to start with a photon in the observers image plane, specified by two (cartesian) co-ordinates <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha" title="\alpha" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\beta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /></a>, at a large distance <a href="https://www.codecogs.com/eqnedit.php?latex=D" target="_blank"><img src="https://latex.codecogs.com/gif.latex?D" title="D" /></a> away from a blackhole of (dimensionless) spin <a href="https://www.codecogs.com/eqnedit.php?latex=a" target="_blank"><img src="https://latex.codecogs.com/gif.latex?a" title="a" /></a>, observed at an inclination angle <a href="https://www.codecogs.com/eqnedit.php?latex=\theta_{\rm&space;obs}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_{\rm&space;obs}" title="\theta_{\rm obs}" /></a>. The co-ordinate settup is displayed in the below figure:  

![](https://github.com/andrewmummery/RayTracer/blob/main/Schematic_of_co-ordinates.png)

The path of these photons are then ran backwards towards the blackhole-disc plane. The path of the photon is terminated when it either: 
* Hits the event horizon of the blackhole
* Hits an accretion disc in the equatorial plane of the blackhole (the presence of the disc can be switched on and off)
* Passes between the disc inner edge and the blackholes event horizon and escapes to a large distance (~ 30 gravitational radii) from the blackhole



### How do I get set up? ###

**Method one:**
* Pull the codes into their own folder
* Run examples.py for a series of examples of what the code package can do. Note: you will have to manually edit examples.py to display different examples, the script defaults to 'example_zero' 
* Make your own images/animations by creating Python scripts in the same folder as RayTracer.py

