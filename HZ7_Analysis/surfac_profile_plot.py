"""Surface Density Plotting"""


from fits_images import FitsImage
from sersic_fitter import SersicFitter

def apply_sersic_fits(fits_object: FitsImage, center, radii):
    """Fits sersic profiles to surface profile data for a given fits image."""
    x_vals, y_vals, y_uncertainties, areas = fits_object.generate_radial_profile(center, radii)
    fit_free, fit_n1 = SersicFitter(x_vals, y_vals / areas, y_uncertainties / areas).fit_sersics()
    
    return fit_free, fit_n1
