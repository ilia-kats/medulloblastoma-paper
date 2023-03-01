import tensorflow as tf
import gpflow

class FullRandomEffectsModel(tf.Module):
    def __init__(
        self,
        Y1,  # ct1 density vector
        Y2,  # ct2 density vector
        X,  # spot coordinates
        l_value,  # numerical value of l parameter to start with (then will converge)
    ):
        self.Y1 = np.squeeze(Y1)
        self.Y2 = np.squeeze(Y2)

        self.Y1 = (self.Y1 - np.mean(self.Y1)) / np.std(self.Y1)
        self.Y2 = (self.Y2 - np.mean(self.Y2)) / np.std(self.Y2)

        self.Y = tf.concat((self.Y1, self.Y2), axis=0)
        self.coords = tf.cast(X, dtype=tf.float64)
        self.e = tf.Variable(tf.ones(shape=(2, 2), dtype=tf.float64))
        self.sigma = gpflow.Parameter(
            [1, 1], dtype=tf.float64, transform=gpflow.utilities.positive(lower=1e-6)
        )

        self.kernel = gpflow.kernels.SquaredExponential()
        gpflow.utilities.set_trainable(self.kernel.variance, False)
        self.kernel.lengthscales = gpflow.Parameter(
            value=l_value, transform=gpflow.utilities.positive(lower=1e-1)
        )

        self.e_null = gpflow.Parameter(
            [1, 1], dtype=tf.float64, transform=gpflow.utilities.positive(lower=1e-6)
        )
        self.sigma_null = gpflow.Parameter(
            [1, 1], dtype=tf.float64, transform=gpflow.utilities.positive(lower=1e-6)
        )
        self.kernel_null = gpflow.kernels.SquaredExponential()
        gpflow.utilities.set_trainable(self.kernel_null.variance, False)
        self.kernel_null.lengthscales = gpflow.Parameter(
            value=l_value, transform=gpflow.utilities.positive(lower=1e-1)
        )

    def covmat(self, null=False):
        bigeye = tf.linalg.LinearOperatorIdentity(self.Y1.shape[0], dtype=tf.float64)

        if not null:
            spatialvar = tf.linalg.LinearOperatorLowerTriangular(self.e)
            spatialvar = spatialvar @ spatialvar.adjoint()
            kern = self.kernel

            noisevar = tf.linalg.LinearOperatorDiag(self.sigma)
        else:
            spatialvar = tf.linalg.LinearOperatorDiag(self.e_null)
            kern = self.kernel_null
            noisevar = tf.linalg.LinearOperatorDiag(self.sigma_null)

        spatialvar = tf.linalg.LinearOperatorKronecker(
            (spatialvar, tf.linalg.LinearOperatorFullMatrix(kern(self.coords)))
        )
        noisevar = tf.linalg.LinearOperatorKronecker((noisevar, bigeye))
        return spatialvar.to_dense() + noisevar.to_dense()

    def loglik(self, null=False):
        var = self.covmat(null)
        chol = tf.linalg.cholesky(var)
        ldet = tf.reduce_sum(2 * tf.math.log(tf.linalg.diag_part(chol)))

        residual = self.Y
        quad = tf.reduce_sum(
            tf.square(tf.linalg.cholesky_solve(chol, residual[:, tf.newaxis]))
        )

        return -0.5 * ldet - 0.5 * quad

    def trainable_vars(self, what="all"):
        if what == "all":
            return self.trainable_variables
        elif what == "alt":
            return (
                (self.e,)
                + self.sigma.trainable_variables
                + self.kernel.trainable_variables
            )
        elif what == "null":
            return (
                self.e_null.trainable_variables
                + self.sigma_null.trainable_variables
                + self.kernel_null.trainable_variables
            )
