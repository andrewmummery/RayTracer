### For compiling/running the c++ file
import os
import platform ## for checking if im windows or linux
import sys ## for printing progress on large loops (camera image function).

### 3rd party imports ###

## Maths and linear algebra
import numpy as np
from scipy.interpolate import interp1d##Interpolation useful for making smooth anomations

### Plotting and animating
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d as ax3d
from matplotlib import animation
import matplotlib.colors


def printProgressBar(Q,size,preText='',postText=''):
    """
        Prints a progress bar to show how far through a loop you are. 
    """
    n_bar =30 #size of progress bar
    q = Q/size
    sys.stdout.write('\r')
    sys.stdout.write(f" {preText} [{'=' * int(n_bar * q):{n_bar}s}] {postText} ")
    sys.stdout.flush()
    return None

def compile_cpp(cpp_file='ray_tracing.cpp', exe_file='ray_tracing.o'):
    ''' Compiles the C++ file.
    '''
    os.system("echo Compiling " + cpp_file)
    os.system('g++ ' + cpp_file + ' -o ' + exe_file)
    print('Done.')
    return exe_file

def run_cpp_ray_trace(params, technique='Simple',exe_file='ray_tracing.o', cpp_file='ray_tracing.cpp', force_new_compile=False, get_times=False, get_j=False):
    ''' I run the C++ script ray_tracing.cpp, either by first compiling ray_tracing.o
        or by simply running ray_tracing.o if it already exists.
        
        I take as paramaters 
        
        params[0] = alpha_0 = X co-ordinate in image plane
        params[1] = beta_0 = Y co-ordinate in image plane
        params[2] = theta_0 = image plane inclination angle (degrees)
        params[3] = a = Normalised black hole spin -1 < a < 1. 
        
        technique is the name of logic used for ray tracing.
    '''
    if force_new_compile:
        exe_file = compile_cpp(cpp_file, exe_file)
    elif os.path.isfile(exe_file):
        pass
    else:
        exe_file = compile_cpp(cpp_file, exe_file)
    
    if platform.system() == 'Darwin':
        exe_file = './' + exe_file
    elif platform.system() == 'Windows':
        exe_file = '.\\' + exe_file
    else:
        exe_file = './' + exe_file#I am probably on Linux
    
    if len(params) != 4:
        print('######################### ERROR #########################')
        print('Ray trace got %d params but needs 4'%(len(params))) 
        print('#########################################################')
        raise ValueError('Ray trace got %d params but needs 4'%(len(params)))
    
    
    params.append(technique)

    allowed_techniques = ['Simple', 'Disc', 'NoDisc']
    if params[-1] not in allowed_techniques:
        print('################################### WARNING ###################################')
        print('Ray trace got %s technique but must be one of %s , %s , %s .'%(technique, *allowed_techniques)) 
        print('Defaulting to Simple technique.')
        print('###############################################################################')
        params[-1] = 'Simple'
        
    params = [str(param) for param in params]
    
    pipe = os.popen(exe_file + ' ' + ' '.join(params))
    ## Note that os.popen() is not technically the recommended technique for doing this.
    ## See RunCpp.py for a better way. 
    
    if get_times:
        if get_j:
            x, y, z, t, j = process_cpp_output(pipe,a=float(params[-2]),get_times=get_times, get_j = get_j)
            return x, y, z, t, j
        x, y, z, t = process_cpp_output(pipe,a=float(params[-2]),get_times=get_times, get_j = get_j)
        return x, y, z, t
    if get_j:
        x, y, z, j = process_cpp_output(pipe,a=float(params[-2]),get_times=get_times, get_j = get_j)
        return x, y, z, j
    x, y, z = process_cpp_output(pipe, a=float(params[-2]))
    return x, y, z
    
def piped_line_to_array(line):
    ''' Takes an output line and returns an array.
        The C++ code is written to return the r(t), 
        theta(t) and phi(t) arrays as a single string 
        with values seperated by commas. This transforms
        that back into the original array.
    '''
    arr = []
    l = 0
    for k in range(len(line)):
        if line[k] == ',':
            arr.append(float(line[l:k]))
            l=k+1
    return np.array(arr)

def process_cpp_output(pipe, a, get_times=False, get_j = False):
    ''' Transforms the output from the C++ file into 
        useful numpy arrays. Returns x(t), y(t), z(t).
        Can return the time as well if get_times=True.
    '''
    line_r = pipe.readline()# C++ code written so that 
    line_theta = pipe.readline()# the first line of output
    line_phi = pipe.readline()# is r, then theta, then phi, then t. 
    line_t = pipe.readline()
    line_j = pipe.readline()
    
    r = piped_line_to_array(line_r)
    theta = piped_line_to_array(line_theta)
    phi = piped_line_to_array(line_phi) 
    
    if get_times:
        t = piped_line_to_array(line_t)
    
    if get_j:
        j = float(line_j)
    
    
    x, y, z = cartesian_from_sphereical_polar(r, theta, phi, a)
    
    if get_times:
        if get_j:
            return x, y, z, t, j
        return x, y, z, t
    if get_j:
        return x, y, z, j
    return x, y, z

def cartesian_from_sphereical_polar(r, theta, phi, a):
    ''' Transforms sphereical polar co-ordinates 
        (r, theta, phi) to cartesian co-ordinates 
        (x, y, z).
    '''
    x = np.sqrt(a*a + r*r)*np.sin(theta)*np.cos(phi)
    y = np.sqrt(a*a + r*r)*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return x, y, z

def sphereical_polar_from_cartesian(x, y, z, a):
    ''' Transforms cartesian co-ordinates 
        (x, y, z) to sphereical polar
        co-ordinates (r, theta, phi).
    '''
    phi = np.arctan(y/x)
    R = np.sqrt(x**2 + y**2 + z**2)
    r = np.sqrt(1/2 * (R**2 - a**2 + np.sqrt((R**2 - a**2)**2 + 4*a**2*z**2)))
    theta = np.arccos(z/r)
    return r, theta, phi


def get_isco(a):
	""" Gets the ISCO radius for a given BH spin parameter"""
	Z_1 = 1 + (1-a**2)**(1/3) * ((1+a)**(1/3) + (1-a)**(1/3)) 
	Z_2 = np.sqrt(3*a**2 + Z_1**2)
	return (3 + Z_2 - np.sign(a) * np.sqrt((3-Z_1)*(3 + Z_1 + 2 * Z_2)))

def get_event_horizon(a):
    ''' Gets the event horizon for a given BH spin parameter '''
    return 1 + np.sqrt(1 - a*a) # black hole event horizon

def make_canvas(a, fig_width=9, fig_height=6, view_theta=60,view_phi=-130,axis='off'):
    '''
    Returns a figure and axis with a black hole and
    disc drawn on
    '''
    fig = plt.figure(figsize=(fig_width, fig_height))
    ax = fig.add_subplot(111,projection='3d')
    ax.set_xlim(-15,15)
    ax.set_ylim(-15,15)
    ax.set_zlim(-9,9)
    ax.view_init(90-view_theta,view_phi)
    
    ax = plot_black_hole(a, ax)
    ax = plot_disc(a, ax)
    
    plt.axis(axis)
    
    return fig, ax

def make_canvas_no_disc(a, fig_width=9, fig_height=6, view_theta=60,view_phi=-130,axis='off'):
    '''
    Returns a figure and axis with a black hole and
    disc drawn on
    '''
    fig = plt.figure(figsize=(fig_width, fig_height))
    ax = fig.add_subplot(111,projection='3d')
    ax.set_xlim(-15,15)
    ax.set_ylim(-15,15)
    ax.set_zlim(-9,9)
    ax.view_init(90-view_theta,view_phi)
    
    ax = plot_black_hole(a, ax)
    
    plt.axis(axis)
    
    
    return fig, ax


def plot_black_hole(a, ax):
    """
    plot the black hole
    """
    rH = get_event_horizon(a)
    
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:20j]
    xs = np.sqrt(rH*rH + a*a)*np.cos(u)*np.sin(v)
    ys = np.sqrt(rH*rH + a*a)*np.sin(u)*np.sin(v)
    zs = rH*np.cos(v)

    ax.plot_surface(xs, ys, zs, color="k",alpha=0.5)
    plt.axis('off')    
    
    return ax
    

def plot_disc(a, ax):
    """
    plot the disc -- p = outer edge radius, q = inner edge radius
    """
    rI = get_isco(a)
    
    N = 100

    thetad = np.linspace(0, 2.*np.pi, N)
    phid = np.linspace(0, 2.*np.pi, N)
    thetad, phid = np.meshgrid(thetad, phid)
    p, q = 20,rI
    c, b = (p+q)/2, (p-q)/2
    xd = (c + b*np.cos(thetad)) * np.cos(phid)
    yd = (c + b*np.cos(thetad)) * np.sin(phid)
    zd = 0.000000001*b * np.sin(thetad)

    ax.contourf(xd,yd,zd,[0.0000000001,2],zdir='z',cmap=cm.autumn,alpha = 0.5)
    
    return ax
    
def plot_ray(ax, x, y, z, color=None):
    ''' Plots the photon path. 
    '''
    ax.plot(x, y, z, color=color)
    return


def animate_rays(xs, ys, zs, ts, a, cmap='jet', n_frame=None, disc=False, burst_mode=False, lw=2.0, ls='-', interval = 1, blit = False, repeat = True,
     fig_width=9, fig_height=6, view_theta=60,view_phi=-130):
    ''' Animates a set of ray paths described by x_i(t_i), y_i(t_i), z_i(t_i). 
        Takes as input:
            xs = [x1(t1), x2(t2), x3(t3), ....]
            ys = [y1(t1), y2(t2), y3(t3), ....]
            zs = .... etc. 
    '''
    
    if n_frame is None:
        n_frame = int(max([len(x) for x in xs])/3)
        
    n_rays = len(xs)
    
    if cmap is None:
        colors = [None for k in range(n_rays)]
    else:
        cm = plt.cm.get_cmap(cmap)
        colors = cm(np.linspace(0,1,n_rays))

    def get_t_i(x, y, z, t, r0 = 30):
        ''' Returns the time when r ~ r0, 
            and the array index when that happens.
        '''
        r,_,_ = sphereical_polar_from_cartesian(x, y, z, a)
        ind = np.argmin(r>r0) # gets first instance of r < r0 in array.
        t_i = t[ind]
        return t_i, ind
    
    def get_t_max(ts):
        ''' Returns the largest proper time 
            of all the photon paths
        '''
        tm = max([t[-1] for t in ts])
        return tm
        
    x_animate = [[] for k in range(n_rays)]
    y_animate = [[] for k in range(n_rays)]
    z_animate = [[] for k in range(n_rays)]
    
    lamb_animate = np.linspace(0,1,n_frame)
    
    t_j = get_t_max(ts)
    
    for k in range(n_rays):
        t_i, ind = get_t_i(xs[k],ys[k],zs[k],ts[k])
        t_m_i = max(ts[k])
        if t_m_i != t_j:
            x_interp = interp1d((np.append(ts[k][ind:],t_j)-t_i)/(t_j-t_i), np.append(xs[k][ind:], xs[k][-1]))
            y_interp = interp1d((np.append(ts[k][ind:],t_j)-t_i)/(t_j-t_i), np.append(ys[k][ind:], ys[k][-1]))
            z_interp = interp1d((np.append(ts[k][ind:],t_j)-t_i)/(t_j-t_i), np.append(zs[k][ind:], zs[k][-1]))
        else:
            x_interp = interp1d((ts[k][ind:]-t_i)/(t_j-t_i), xs[k][ind:])
            y_interp = interp1d((ts[k][ind:]-t_i)/(t_j-t_i), ys[k][ind:])
            z_interp = interp1d((ts[k][ind:]-t_i)/(t_j-t_i), zs[k][ind:])

        x_animate[k] = x_interp(lamb_animate)
        y_animate[k] = y_interp(lamb_animate)
        z_animate[k] = z_interp(lamb_animate)


    
    if disc:
        fig, ax = make_canvas(a, fig_width=fig_width, fig_height=fig_height, view_theta=view_theta,view_phi=view_phi)
    else:
        fig, ax = make_canvas_no_disc(a, fig_width=fig_width, fig_height=fig_height, view_theta=view_theta,view_phi=view_phi)
    
    line, = ax.plot([], [])
    plotlays, plotcols = [n_rays], colors

    lines = []
    for index in range(n_rays):
        lobj = ax.plot([],[],lw=lw,color=plotcols[index],linestyle=ls)[0]
        lines.append(lobj)

    def init():
        for line in lines:
            line.set_data(np.asarray([]),np.asarray([]))
            line.set_3d_properties(np.asarray([]))
        return lines

    def animate(i):
        xlist = [[] for k in range(n_rays)]
        ylist = [[] for k in range(n_rays)]
        zlist = [[] for k in range(n_rays)]
        if burst_mode:
            for j in range(n_rays):
                xlist[j] = x_animate[j][max(0,i-int(n_frame/10)):i]
                ylist[j] = y_animate[j][max(0,i-int(n_frame/10)):i]
                zlist[j] = z_animate[j][max(0,i-int(n_frame/10)):i]
        else:
            for j in range(n_rays):
                xlist[j] = x_animate[j][:i]
                ylist[j] = y_animate[j][:i]
                zlist[j] = z_animate[j][:i]

        for lnum,line in enumerate(lines):
            line.set_data(np.asarray(xlist[lnum]), np.asarray(ylist[lnum]))
            line.set_3d_properties(np.asarray(zlist[lnum]))
        return lines

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames = n_frame, interval = interval, blit = blit, repeat = repeat)
    
    return fig, ax, anim


def wavelength_to_rgb(wavelength, gamma=0.8, alpha = 0.1, set_unobservable_grey=True):
    ''' Modified from http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
    This converts a given wavelength of light to an 
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html

    Additionally allows any wavelengths out of visible range to be 
    set to grey with alpha value = alpha.  
    '''
        
    wavelength = float(wavelength)
    if wavelength >= 380 and wavelength <= 750:
        A = 1.0
    else:
        if set_unobservable_grey:
            A = alpha
        else:
            A = 1.0
    if wavelength < 380:
        if set_unobservable_grey:
            wavelength = 379.0
        else:
            wavelength = 380.0
    
    if wavelength >750:
        if set_unobservable_grey:
            wavelength = 751.0
        else:
            wavelength = 750.0
    
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    return (R,G,B,A)

def get_spectral_color_map(set_unobservable_grey=True):
    ''' Returns a spectral color map
    '''
    clim=(350,780)
    norm = plt.Normalize(*clim)
    wl = np.arange(clim[0],clim[1]+1,2)
    colorlist = list(zip(norm(wl),[wavelength_to_rgb(w, set_unobservable_grey=set_unobservable_grey) for w in wl]))
    spectralmap = matplotlib.colors.LinearSegmentedColormap.from_list("spectrum", colorlist)
    return spectralmap

def camera_image(a, theta0, N_r=50, N_phi=200, r_out=20, rest_wavelength=550, wavelength_function=None, set_unobservable_grey=False, set_intensity=False, print_progress=False):
    ''' Takes a camera image. 
    '''
    
    save_string = ['N_r=',str(N_r),'N_phi=',str(N_phi),'a=',str(a),'theta=',str(theta0),'r_out=',str(r_out)]
    save_string = ''.join(save_string)
        
    files_exist = os.path.isfile('im_x'+save_string+'.npy')
    files_exist = files_exist and os.path.isfile('im_y'+save_string+'.npy')
    files_exist = files_exist and os.path.isfile('phys_y'+save_string+'.npy')
    files_exist = files_exist and os.path.isfile('phys_y'+save_string+'.npy')
    files_exist = files_exist and os.path.isfile('red_shift'+save_string+'.npy')
    
    if files_exist:
        print('Loading physical parameters...')
        ri = get_isco(a)
        rH = get_event_horizon(a)
        
        im_x = np.load('im_x'+save_string+'.npy')
        im_y = np.load('im_y'+save_string+'.npy')
        phys_x = np.load('phys_x'+save_string+'.npy')
        phys_y = np.load('phys_y'+save_string+'.npy')
        red_shifts = np.load('red_shift'+save_string+'.npy')
    else:
        print('Computing physical parameters of %d photons... \n This may take a couple of minutes. '%(N_r*N_phi))
        def U0(r,a):
            return (1+a*np.power(r,-3/2))/np.sqrt(1-3/r+2*a*np.power(r,-3/2))
        def Uphi(r,a):
            return 1/(np.power(r,3/2)*np.sqrt(1-3/r+2*a*np.power(r,-3/2)))
        def g(r,a,j,ri,r_out=20):
            if r < ri:
                return 10
            elif r > r_out:
                return 0
            return 1/(U0(r,a)) * 1/(1 + j*(Uphi(r,a)/U0(r,a))) 
    
    
        ri = get_isco(a)
        rH = get_event_horizon(a)
    
        im_r = [ri + (r_out + 5 - ri) * k/(N_r-1) for k in range(N_r)]
        im_phi = [np.pi/2 + 2*np.pi * k/(N_phi-1) for k in range(N_phi)]

        im_x = []
        im_y = []
        red_shifts = []
        phys_x = []
        phys_y = []
    
        
        f_min = 0.01
        f_max = 2.00
        n_rays_traced = 0
        for i in range(N_r):
            for j in range(N_phi):
                n_rays_traced += 1
                alpha = im_r[i] * np.cos(im_phi[j])
                beta = im_r[i] * np.sin(im_phi[j]) * np.cos(theta0*np.pi/180)
                if beta > 0:
                    beta *= np.power(1/np.cos(theta0*np.pi/180), .7)
                im_x.append(alpha)
                im_y.append(beta)
                x, y, z, j = run_cpp_ray_trace([alpha, beta, theta0, a], get_times=False, get_j = True)
                rf, _, _ = sphereical_polar_from_cartesian(x[-1], y[-1], z[-1], a)## gets the final radial co-ordinate
                if rf < rH:
                    f = 10
                elif z[-1] > 0.5:
                    f = 10
                else:
                    f = g(rf, a, j, ri, r_out)
                red_shifts.append(f)
                phys_x.append(x[-1])
                phys_y.append(y[-1])
                if print_progress:
                    if n_rays_traced % 100 == 0:
                        printProgressBar(Q=n_rays_traced,size=N_r*N_phi,preText='Ray tracing....',postText='Photons traced: '+str(n_rays_traced)+'/'+str(N_r*N_phi))
        
        if print_progress:            
            printProgressBar(Q=N_r*N_phi,size=N_r*N_phi,preText='Ray tracing....',postText='Photons traced: '+str(N_r*N_phi)+'/'+str(N_r*N_phi))
            print('')
        
        red_shifts = [red_shifts[k] if not (red_shifts[k] < f_min or red_shifts[k] > f_max) else np.nan for k in range(N_r*N_phi)]
    
        save_string = ['N_r=',str(N_r),'N_phi=',str(N_phi),'a=',str(a),'theta=',str(theta0),'r_out=',str(r_out)]
        save_string = ''.join(save_string)
        np.save('im_x'+save_string+'.npy',im_x)
        np.save('im_y'+save_string+'.npy',im_y)
        np.save('phys_x'+save_string+'.npy',phys_x)
        np.save('phys_y'+save_string+'.npy',phys_y)
        np.save('red_shift'+save_string+'.npy',red_shifts)
        # End else. 
    
    print('Done. ')
    ### Plotting from here. 
    if wavelength_function is not None:
        wavelengths = [wavelength_function(phys_x[k], phys_y[k])/red_shifts[k] if not np.isnan(red_shifts[k]) else np.nan for k in range(N_r*N_phi)]
    else:
        wavelengths = [rest_wavelength/red_shifts[k] if not np.isnan(red_shifts[k]) else np.nan for k in range(N_r*N_phi)]
    
    imx_lim_u = max([im_x[k] if not np.isnan(wavelengths[k]) else 0 for k in range(N_r*N_phi)])
    imx_lim_l = min([im_x[k] if not np.isnan(wavelengths[k]) else 0 for k in range(N_r*N_phi)])
    imy_lim_u = max([im_y[k] if not np.isnan(wavelengths[k]) else 0 for k in range(N_r*N_phi)])
    imy_lim_l = min([im_y[k] if not np.isnan(wavelengths[k]) else 0 for k in range(N_r*N_phi)])
        
    spectralmap = get_spectral_color_map(set_unobservable_grey=set_unobservable_grey)
    
    fig = plt.figure(figsize=(9,6))
    ax = fig.add_subplot(111)
    
    if set_intensity:
        colors = spectralmap((np.array([wavelengths[k] if not np.isnan(wavelengths[k]) else 0 for k in range(N_r*N_phi)]) - 350)/(780-350))
        colors[:,3] = [red_shifts[k]**3 if not np.isnan(red_shifts[k]) else 0.01 for k in range(N_r*N_phi)]
        colors[:,3] = [.999 * (colors[k,3] - min(colors[:,3])) /(max(colors[:,3]) - min(colors[:,3]))  for k in range(N_r*N_phi)]
        sc = ax.scatter(im_x,im_y,s=100,color=colors, edgecolors=None)
    else:
        sc = ax.scatter(im_x,im_y,s=100,c=wavelengths,cmap = spectralmap, vmin=349, vmax=781, edgecolors=None)
            
    ax.set_xlabel(r'$X_{im}$',fontsize=20)
    ax.set_ylabel(r'$Y_{im}$',fontsize=20,rotation = 0)

    ax.set_xlim(imx_lim_l-1.5, imx_lim_u+1.5)
    ax.set_ylim(imy_lim_l-1.0, imy_lim_u+1.0)

    fig1 = plt.figure(figsize=(9, 6))
    ax1 = fig1.add_subplot(111)

    ax1.set_xlabel(r'$X_{im}$',fontsize=20)
    ax1.set_ylabel(r'$Y_{im}$',fontsize=20,rotation = 0)
    
    ax1.set_xlim(imx_lim_l-1.5, imx_lim_u+1.5)
    ax1.set_ylim(imy_lim_l-1.0, imy_lim_u+1.0)
    
    if wavelength_function is None:
        plot_phi = np.linspace(0,2*np.pi,100)

        disc_x_in = ri * np.cos(plot_phi)
        disc_y_in = ri * np.cos(theta0*np.pi/180) * np.sin(plot_phi)

        disc_x_out = r_out * np.cos(plot_phi)
        disc_y_out = r_out * np.cos(theta0*np.pi/180) * np.sin(plot_phi)

        ax1.plot(disc_x_in, disc_y_in, color=wavelength_to_rgb(rest_wavelength),lw=0)
        ax1.plot(disc_x_out, disc_y_out, color=wavelength_to_rgb(rest_wavelength),lw=0)
        ax1.fill(np.append(disc_x_in, disc_x_out[::-1]), np.append(disc_y_in, disc_y_out[::-1]), color=wavelength_to_rgb(rest_wavelength))

        circle1 = plt.Circle((0, 0), rH, color='k')

        ax1.add_artist(circle1)

        if rH > ri * np.cos(theta0*np.pi/180):
            n = 40
            for k in range(n):
                plot_phi = np.linspace(np.pi/2,3*np.pi/2,100)
                r = ri + (rH/np.cos(theta0*np.pi/180) - ri) * k/(n-1)
                disc_x_line = r * np.sin(plot_phi)
                disc_y_line = r * np.cos(theta0*np.pi/180) * np.cos(plot_phi)
                ax1.plot(disc_x_line, disc_y_line, color=wavelength_to_rgb(rest_wavelength),lw=1)
    else:
        n = 40
        for k in range(n):
            plot_phi = np.linspace(0,2*np.pi,200)
            r = ri + (r_out - ri) * k/(n-1)
            disc_x_line = r * np.cos(plot_phi)
            disc_y_line = r * np.cos(theta0*np.pi/180) * np.sin(plot_phi)
            ax1.scatter(disc_x_line, disc_y_line,s=100, c=wavelength_function(disc_x_line, disc_y_line/np.cos(theta0*np.pi/180)),cmap=spectralmap,vmin=350,vmax=780)
    
        if rH > ri * np.cos(theta0*np.pi/180):
            for k in range(n):#### Can do this better with ax1.fill(), like above, with semicircle function etc.
                plot_phi = np.linspace(0,np.pi,200)
                r = rH * k/(n-1)
                disc_x_line = r * np.cos(plot_phi)
                disc_y_line = r * np.sin(plot_phi)
                ax1.scatter(disc_x_line, disc_y_line,s=1, color='k')
        else:
            for k in range(n):#### Can do this better with ax1.fill(), like above, with semicircle function etc.
                plot_phi = np.linspace(0,2*np.pi,400)
                r = rH * k/(n-1)
                disc_x_line = r * np.cos(plot_phi)
                disc_y_line = r * np.sin(plot_phi)
                ax1.scatter(disc_x_line, disc_y_line,s=1, color='k')
        
    
    ax.set_title('General Relativistic Universe', fontsize=18)
    ax1.set_title("'Newtonian' Universe", fontsize=18)
    ### End of camera image.
    return fig, ax, fig1, ax1


def plot_rays_from_parameters(spins, thetas, alphas, betas, technique='NoDisc', fig_width=9, fig_height=6, view_theta=60,view_phi=-130):
    
    ## Cast all parameters into lists
    if type(spins) != type([]):
        spins = [spins]
    if type(thetas) != type([]):
        thetas = [thetas]
    if type(alphas) != type([]):
        alphas = [alphas]
    if type(betas) != type([]):
        betas = [betas]
    
    ## Check they have the right lengths
    max_len = max([len(spins), len(thetas), len(alphas), len(betas)])
    len_test = (len(spins) == max_len or len(spins) == 1)
    len_test = len_test and (len(thetas) == max_len or len(thetas) == 1)
    len_test = len_test and (len(alphas) == max_len or len(alphas) == 1)
    len_test = len_test and (len(betas) == max_len or len(betas) == 1)
    
    if not len_test:
        raise ValueError('Initial ray parameters must be lists of length 1 or equal length. \n You have provided lists with lengths %d, %d, %d, %d.'
        %(len(spins), len(thetas), len(alphas), len(betas)))
    
    if len(spins) == 1:
        spins = [spins[0] for k in range(max_len)]
    if len(thetas) == 1:
        thetas = [thetas[0] for k in range(max_len)]
    if len(alphas) == 1:
        alphas = [alphas[0] for k in range(max_len)]
    if len(betas) == 1:
        betas = [betas[0] for k in range(max_len)]
    
    n_ray = max_len
    xs, ys, zs, ts = [[] for k in range(n_ray)], [[] for k in range(n_ray)], [[] for k in range(n_ray)], [[] for k in range(n_ray)]

    for i, a0, b0, theta0, a in zip(range(n_ray), alphas,  betas, thetas, spins):
        x, y, z, t = run_cpp_ray_trace([a0, b0, theta0, a], technique, get_times=True)
        xs[i], ys[i], zs[i], ts[i] = x, y, z, t
    
    if technique == 'NoDisc':
        fig, ax = make_canvas_no_disc(spins[0], fig_width=fig_width, fig_height=fig_height, view_theta=view_theta,view_phi=view_phi)
    else:
        fig, ax = make_canvas(spins[0], fig_width=fig_width, fig_height=fig_height, view_theta=view_theta,view_phi=view_phi)
    
    cm = plt.cm.get_cmap('jet')
    colors = cm(np.linspace(0,1,n_ray))

    for i, x, y, z in zip(range(n_ray), xs, ys, zs):
        plot_ray(ax, x, y, z,color=colors[i])

    return fig, ax

def animate_rays_from_parameters(spins, thetas, alphas, betas, technique='NoDisc', cmap='jet', n_frame=None, burst_mode=False,
                lw=2.0, ls='-', interval = 1, blit = False, repeat = True, fig_width=9, fig_height=6, view_theta=60,view_phi=-130):
    
    ## Cast all parameters into lists
    if type(spins) != type([]):
        spins = [spins]
    if type(thetas) != type([]):
        thetas = [thetas]
    if type(alphas) != type([]):
        alphas = [alphas]
    if type(betas) != type([]):
        betas = [betas]
    
    ## Check they have the right lengths
    max_len = max([len(spins), len(thetas), len(alphas), len(betas)])
    len_test = (len(spins) == max_len or len(spins) == 1)
    len_test = len_test and (len(thetas) == max_len or len(thetas) == 1)
    len_test = len_test and (len(alphas) == max_len or len(alphas) == 1)
    len_test = len_test and (len(betas) == max_len or len(betas) == 1)
    
    if not len_test:
        raise ValueError('Initial ray parameters must be lists of length 1 or equal length. \n You have provided lists with lengths %d, %d, %d, %d.'
        %(len(spins), len(thetas), len(alphas), len(betas)))
    
    if len(spins) == 1:
        spins = [spins[0] for k in range(max_len)]
    if len(thetas) == 1:
        thetas = [thetas[0] for k in range(max_len)]
    if len(alphas) == 1:
        alphas = [alphas[0] for k in range(max_len)]
    if len(betas) == 1:
        betas = [betas[0] for k in range(max_len)]
    
    n_ray = max_len
    xs, ys, zs, ts = [[] for k in range(n_ray)], [[] for k in range(n_ray)], [[] for k in range(n_ray)], [[] for k in range(n_ray)]

    for i, a0, b0, theta0, a in zip(range(n_ray), alphas,  betas, thetas, spins):
        x, y, z, t = run_cpp_ray_trace([a0, b0, theta0, a],technique, get_times=True)
        xs[i], ys[i], zs[i], ts[i] = x, y, z, t
    
    if technique == 'NoDisc':
        anim_fig, anim_ax, anim = animate_rays(xs, ys, zs, ts, spins[0], disc=False, cmap=cmap, n_frame=n_frame, burst_mode=burst_mode,
     lw=lw, ls=ls, interval = interval, blit = blit, repeat = repeat, fig_width=fig_width, fig_height=fig_height, view_theta=view_theta,view_phi=view_phi)
    else:
        anim_fig, anim_ax, anim = animate_rays(xs, ys, zs, ts, spins[0], disc=True, cmap=cmap, n_frame=n_frame, burst_mode=burst_mode,
     lw=lw, ls=ls, interval = interval, blit = blit, repeat = repeat, fig_width=fig_width, fig_height=fig_height, view_theta=view_theta,view_phi=view_phi)

    return anim_fig, anim_ax, anim

def plot_and_animate_rays_from_parameters(spins, thetas, alphas, betas, technique='NoDisc', cmap='jet', n_frame=None, burst_mode=False,
                lw=2.0, ls='-', interval = 1, blit = False, repeat = True, fig_width=9, fig_height=6, view_theta=60,view_phi=-130):
    
    ## Cast all parameters into lists
    if type(spins) != type([]):
        spins = [spins]
    if type(thetas) != type([]):
        thetas = [thetas]
    if type(alphas) != type([]):
        alphas = [alphas]
    if type(betas) != type([]):
        betas = [betas]
    
    ## Check they have the right lengths
    max_len = max([len(spins), len(thetas), len(alphas), len(betas)])
    len_test = (len(spins) == max_len or len(spins) == 1)
    len_test = len_test and (len(thetas) == max_len or len(thetas) == 1)
    len_test = len_test and (len(alphas) == max_len or len(alphas) == 1)
    len_test = len_test and (len(betas) == max_len or len(betas) == 1)
    
    if not len_test:
        raise ValueError('Initial ray parameters must be lists of length 1 or equal length. \n You have provided lists with lengths %d, %d, %d, %d.'
        %(len(spins), len(thetas), len(alphas), len(betas)))
    
    if len(spins) == 1:
        spins = [spins[0] for k in range(max_len)]
    if len(thetas) == 1:
        thetas = [thetas[0] for k in range(max_len)]
    if len(alphas) == 1:
        alphas = [alphas[0] for k in range(max_len)]
    if len(betas) == 1:
        betas = [betas[0] for k in range(max_len)]
    
    n_ray = max_len
    xs, ys, zs, ts = [[] for k in range(n_ray)], [[] for k in range(n_ray)], [[] for k in range(n_ray)], [[] for k in range(n_ray)]

    for i, a0, b0, theta0, a in zip(range(n_ray), alphas,  betas, thetas, spins):
        x, y, z, t = run_cpp_ray_trace([a0, b0, theta0, a],technique, get_times=True)
        xs[i], ys[i], zs[i], ts[i] = x, y, z, t
    
    if technique == 'NoDisc':
        plot_fig, plot_ax = make_canvas_no_disc(spins[0], fig_width=fig_width, fig_height=fig_height, view_theta=view_theta,view_phi=view_phi)
    else:
        plot_fig, plot_ax = make_canvas(spins[0], fig_width=fig_width, fig_height=fig_height, view_theta=view_theta,view_phi=view_phi)
    
    cm = plt.cm.get_cmap('jet')
    colors = cm(np.linspace(0,1,n_ray))

    for i, x, y, z in zip(range(n_ray), xs, ys, zs):
        plot_ray(plot_ax, x, y, z,color=colors[i])
    
    if technique == 'NoDisc':
        anim_fig, anim_ax, anim = animate_rays(xs, ys, zs, ts, spins[0], disc=False, cmap=cmap, n_frame=n_frame, burst_mode=burst_mode,
     lw=lw, ls=ls, interval = interval, blit = blit, repeat = repeat, fig_width=fig_width, fig_height=fig_height, view_theta=view_theta,view_phi=view_phi)
    else:
        anim_fig, anim_ax, anim = animate_rays(xs, ys, zs, ts, spins[0], disc=True, cmap=cmap, n_frame=n_frame, burst_mode=burst_mode,
     lw=lw, ls=ls, interval = interval, blit = blit, repeat = repeat, fig_width=fig_width, fig_height=fig_height, view_theta=view_theta,view_phi=view_phi)
        
    return plot_fig, plot_ax, anim_fig, anim_ax, anim


#########################   END.    #########################



