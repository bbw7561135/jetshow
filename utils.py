import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity


def fit_kde(data, bandwidths_to_search=None):
    """
    Fit 1D density with gaussian KDE.

    :param data:
        Numpy array.
    :param bandwidths_to_search: (optional)
        Iterable of bandwidths to check with 5-fold CV. If ``None`` then use
        ``np.logspace(-3, 1, 10)``. (default: ``None``)
    :return:
        Instance of ``sklearn.neighbors.KernelDensity`` class with best density
        estimate chosen by 5-fold CV.
    """
    if bandwidths_to_search is None:
        bandwidths_to_search = np.logspace(-3, 1, 10)
    params = {'bandwidth': bandwidths_to_search}
    grid = GridSearchCV(KernelDensity(), params, cv=5)
    grid.fit(data[:, np.newaxis])

    return grid.best_estimator_


def generate_grid(x_min, x_max, kde, N):
    from scipy import optimize, integrate
    xs = list()

    def f(x):
        try:
            result = np.exp(density.score_samples(x.reshape(1, -1)))
        except AttributeError:
            result = np.exp(density.score_samples(np.array(x).reshape(1, -1)))
        return result

    while True:
        x = optimize.fsolve(lambda x: integrate.quad(f, x_min, x)[0]-1./N, 0)[0]
        x_min = x
        xs.append(x)
        if x > x_max:
            break

    return np.array(xs)