"""
    Some simple examples of what RayTracer.py can do. 
"""
## Import for the ray tracing modules
import RayTracer

## 3rd party imports
## numpy for maths and any linear algebra
import numpy as np
## matplotlib for displaying plots
import matplotlib.pyplot as plt

def example_zero():
    a = 0
    theta0 = 85
    a0 = 3
    b0 = 4.25063
    
    plot_fig, plot_ax, anim_fig, anim_ax, anim = RayTracer.plot_and_animate_rays_from_parameters(spins=a, thetas=theta0, alphas=a0, betas=b0, technique='NoDisc')
    
    plt.show()
    
    
def example_one():
    a = 0
    theta0 = 85
    a0_s = 3
    
    b0_s = [4.23,4.24,4.245,4.2475,4.25063,4.255,4.2575,4.26,4.27,4.3,4.35,4.4,4.5,4.6,4.7]
    
    plot_fig, plot_ax, anim_fig, anim_ax, anim = RayTracer.plot_and_animate_rays_from_parameters(spins=a, thetas=theta0, alphas=a0_s, betas=b0_s, technique='NoDisc')
    
    plt.show()
    

def example_two_a():
    a = 1
    theta0 = 60
    
    a0_s = [2.24 + k/3000 for k in range(100)]
    b0_s = [2.24 + k/3000 for k in range(100)]
    
    plot_fig, plot_ax, anim_fig, anim_ax, anim = RayTracer.plot_and_animate_rays_from_parameters(spins=a, thetas=theta0, alphas=a0_s, betas=b0_s, technique='NoDisc')
    
    plt.show()


def example_two_b():
    a = 1
    theta0 = 60

    a0_s = [2.3 + k/700 for k in range(100)]
    b0_s = [2.3 + k/700 for k in range(100)]

    plot_fig, plot_ax, anim_fig, anim_ax, anim = RayTracer.plot_and_animate_rays_from_parameters(spins=a, thetas=theta0, alphas=a0_s, betas=b0_s)

    plt.show()        


def example_two_c():
    a = 1.0
    theta0 = 60
    a0_s = [2.24 + k/3000 for k in range(100)]
    b0_s = [2.24 + k/3000 for k in range(100)]

    a0_s_2 = [2.3 + k/700 for k in range(100)]
    b0_s_2 = [2.3 + k/700 for k in range(100)]

    for k in range(100):
        a0_s.append(a0_s_2[k])
        b0_s.append(b0_s_2[k])

    plot_fig, plot_ax, anim_fig, anim_ax, anim = RayTracer.plot_and_animate_rays_from_parameters(spins=a, thetas=theta0, alphas=a0_s, betas=b0_s)

    plt.show()


def example_three():
    a = 0
    theta0 = [5 + k*80/20 for k in range(21)]
    a0_s = 4
    b0_s = 4

    plot_fig, plot_ax, anim_fig, anim_ax, anim = RayTracer.plot_and_animate_rays_from_parameters(spins=a, thetas=theta0, alphas=a0_s, betas=b0_s)

    plt.show()
 
         
def example_four():
    a = [-0.998, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.998]
    theta0 = 60
    a0_s = 4.375
    b0_s = 4.35

    plot_fig, plot_ax, anim_fig, anim_ax, anim = RayTracer.plot_and_animate_rays_from_parameters(spins=a, thetas=theta0, alphas=a0_s, betas=b0_s)

    plt.show()
    
    
def lambda_from_xy_1(x, y):
    """ This has no physical motivation. 
    """
    lambda_0 = 550
    eps = 0.2
    x0, y0 = 5, 5
    
    lambda_tru = lambda_0 * (1 + eps * np.cos(x/x0) * np.sin(y/y0))
    return lambda_tru

def lambda_from_xy_2(x, y):
    ''' 4 "cold spots" on the disc. '''
    lambda_0 = 450
    d_lambda = 200
    x0, y0 = 9, 9
    
    lambda_tru = lambda_0
    lambda_tru += d_lambda * np.exp(-((x-x0)/5)**2) * np.exp(-((y-y0)/5)**2) 
    lambda_tru += d_lambda * np.exp(-((x+x0)/5)**2) * np.exp(-((y-y0)/5)**2) 
    lambda_tru += d_lambda * np.exp(-((x-x0)/5)**2) * np.exp(-((y+y0)/5)**2) 
    lambda_tru += d_lambda * np.exp(-((x+x0)/5)**2) * np.exp(-((y+y0)/5)**2) 
    return lambda_tru

def lambda_from_xy_3(x, y):
    ''' A realistic temperature decay with increasing disc radius. 
        Based on approximate temperature profiles from real disc physics. 
    '''
    r = np.sqrt(x**2 + y**2 - a**2)
    lambda_0 = 500
    
    lambda_tru = lambda_0 * (r/10)**(0.75)
    return lambda_tru

def lambda_from_xy_4(x, y):
    ''' A hot spiral in the disc. '''
    lambda_0 = 650
    r = np.sqrt(x**2 + y**2 - a**2)    
    d_lambda = 250
    spiral_a = 0
    spiral_b = 18/np.pi
    phi = np.arctan(y/x)
    if type(phi) != type(np.zeros(1)):
        if phi < 0:
            phi = phi + np.pi
    else:
        for i, phi_ in enumerate(phi):
            if phi_ < 0:
                phi[i] = phi_ + np.pi
    spiral_par = r - (spiral_a + spiral_b * phi)
    lambda_tru = lambda_0 - d_lambda * np.exp(-abs(spiral_par)/2.5)
    return lambda_tru

example_zero()
# example_one()
# example_two_a()
# example_two_b()
# example_two_c()
# example_three()
# example_four()

a = 0.9
theta0 = 80

# RayTracer.camera_image(a, theta0, set_unobservable_grey=True)
# RayTracer.camera_image(a, theta0, set_unobservable_grey=True, set_intensity=False, rest_wavelength=500)
# RayTracer.camera_image(a, theta0, set_unobservable_grey=True, set_intensity=False, rest_wavelength=600)
# RayTracer.camera_image(a, theta0, wavelength_function=lambda_from_xy_2, set_unobservable_grey=False, set_intensity=False)
# RayTracer.camera_image(a, theta0, wavelength_function=lambda_from_xy_3, set_unobservable_grey=False, set_intensity=False)
# RayTracer.camera_image(a, theta0, wavelength_function=lambda_from_xy_4, set_unobservable_grey=False, set_intensity=False)


plt.show()


# End.
