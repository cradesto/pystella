from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel as C

kernel = C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e1)) + WhiteKernel(noise_level=0.05)
model = GaussianProcessRegressor(alpha=1e-10, copy_X_train=True,
                                 kernel=kernel)
model.fit(x_train, y_train)


x = np.linspace(0, 10, 100).reshape(-1, 1)
y, sigma = model.predict(x, return_std=True)
pyplot.figure()
plot_with_sigma(x_train, y_train, x, y, sigma, 'GP sklearn')
pyplot.show()