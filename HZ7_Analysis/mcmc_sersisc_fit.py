"""
MCMC Fitting routing for Sersic Profiles
"""

import numpy as np
import emcee
import pylab as plt
import corner
from sersic_fitter import SersicFitter, sersic

def model(theta, radius):
    amplitude, effective_radius, sersic_index, offset = theta 
    model = amplitude * np.exp(-(2*sersic_index -1./3) * ((radius/effective_radius)**(1.0/sersic_index) - 1.0)) + offset
    return model 

def lnlike(theta, x, y, y_err):
    LnLike = -0.5 * np.sum(((y - model(theta, x)) / (y_err))**2 )
    return LnLike

def lnprior(theta):
    amplitude, effective_radius, sersic_index, offset = theta 
    if (0 < amplitude < 100) and ( 0 < effective_radius < 100) and (0.3 <sersic_index < 10) and (0 < offset < 100):
        return 0.0
    return -np.inf

def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not lp == 0:
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)

def fit_sersic_mcmc(x, y, y_err):
    nwalkers = 1000
    niter = 5000
    data = (x, y, y_err)
    initial = np.array([np.max(y)/4, np.max(x)/4, 0.75, np.min(y)])
    ndim = len(initial)
    p0 = [initial + 0.01 * np.random.randn(ndim) for i in range(nwalkers)]
    sampler, pos, prob, state = main(p0, nwalkers, niter, ndim, lnprob, data)
    return sampler, pos, prob, state

def main(p0, nwalkers, niter, ndim, lnprob, data):
    print('burning in...')
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=data)
    p0, _, _ = sampler.run_mcmc(p0, 100, progress= True)
    sampler.reset()

    print('running fit...')
    pos, prob, state = sampler.run_mcmc(p0, niter, progress=True)

    return sampler, pos, prob, state 

def sample_walkers(x_for_plot, nsamples,flattened_chain):
    models = []
    draw = np.floor(np.random.uniform(0,len(flattened_chain),size=nsamples)).astype(int)
    thetas = flattened_chain[draw]
    for i in thetas:
        mod = model(i, x_for_plot)
        models.append(mod)
    spread = np.std(models,axis=0)
    med_model = np.median(models,axis=0)
    return med_model,spread

def sersic_fit(x_data, y_data, y_err):
    x_plotting = np.linspace(x_data[0], x_data[-1], 1000)
    sampler, pos, prob, state = fit_sersic_mcmc(x_data, y_data, y_err)
    new_samples = sampler.flatchain
    samples = sampler.get_chain(discard = 100,thin=10, flat = True)

    best_params = []
    uncertainties = []
    for i in range(4):
        mcmc = np.percentile(samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        uncertainties.append(q)
        best_params.append(mcmc[1])
        print (mcmc[1], q[0], q[1])
    best_params = np.array(best_params)
    uncertainties = np.array(uncertainties)
    return best_params, uncertainties

    '''theta_max = np.percentile(samples[:, i], [16, 50, 84])
    print(theta_max)
    for theta in samples[np.random.randint(len(samples), size = 100)]:
        plt.plot(x_plotting, model(theta, x_plotting), color='b', alpha=0.1)
    plt.ylim(0, np.max(y_data) + 5)
    plt.errorbar(x_data, y_data, y_err, fmt = 'ro')
    plt.show()

    plt.errorbar(x_data, y_data, y_err, fmt = 'ro', label = 'data points')
    plt.plot(x_plotting, sersic(x_plotting, *best_params),color='r', lw=3, alpha=0.5, label = 'MCMC Best Fit')
    plt.legend()
    plt.show()

    labels = ['amp', 'reff', 'n', 'offset']
    fig = corner.corner(samples, show_titles=True, labels = labels, plot_datapoints = True, quantiles = [0.16, 0.5, 0.84])
    plt.show()'''