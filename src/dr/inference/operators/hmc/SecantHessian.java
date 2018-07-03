package dr.inference.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.ReadableVector;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
class SecantHessian implements HessianWrtParameterProvider {
    private final int dim;
    private final GradientWrtParameterProvider gradientProvider;
    private final int secantSize;

    private final Secant[] queue;
    private int secantIndex;
    private int secantUpdateCount;
    private double[][] secantHessian;  //TODO: instead store and update the Cholesky decomposition H = LL^T

    SecantHessian(GradientWrtParameterProvider gradientProvider, int secantSize) {
        this.gradientProvider = gradientProvider;
        this.secantSize = secantSize;
        this.dim = gradientProvider.getDimension();

        this.secantHessian = new double[dim][dim]; //TODO: use WrappedVector
        for (int i = 0; i < dim; i++) {
            secantHessian[i][i] = 1.0;
        }

        this.queue = new Secant[secantSize];
        this.secantIndex = 0;
        this.secantUpdateCount = 0;
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {
        throw new RuntimeException("DiagonalHessian not permitted for seacant Hessian approximation, exiting ...");
    }

    @Override
    public double[][] getHessianLogDensity() {
        return secantHessian.clone();
    }

    @Override
    public Likelihood getLikelihood() {
        return gradientProvider.getLikelihood();
    }

    @Override
    public Parameter getParameter() {
        return gradientProvider.getParameter();
    }

    @Override
    public int getDimension() {
        return gradientProvider.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return gradientProvider.getGradientLogDensity();
    }

    private class Secant { // TODO Inner class because we may change the storage depending on efficiency

        ReadableVector gradient;
        ReadableVector position;

        Secant(ReadableVector gradient, ReadableVector position) {
            this.gradient = gradient;
            this.position = position;
        }

        public double getPosition(int index) {
            return position.get(index);
        }

        public double getGradient(int index) {
            return gradient.get(index);
        }
    }

    public void storeSecant(ReadableVector gradient, ReadableVector position) {
        queue[secantIndex] = new Secant(gradient, position);

        if (secantUpdateCount == 0) {
            initializeSecantHessian(queue[secantIndex]);
        } else {
            final int lastSecantIndex = (secantIndex + secantSize - 1) % queue.length;
            updateSecantHessian(queue[secantIndex], queue[lastSecantIndex]);
        }

        secantIndex = (secantIndex + 1) % queue.length;
        ++secantUpdateCount;
    }

    private void initializeSecantHessian(Secant secant) {
            /*
            Equation 6.20 on page 143 of Numerical Optimization 2nd by Nocedal and Wright
             */
        double numerator = 0.0;
        double denominator = 0.0;
        for (int i = 0; i < dim; i++) {
            numerator += secant.getGradient(i) * secant.getPosition(i);
            denominator += secant.getGradient(i) * secant.getGradient(i);
        }

        final double fraction = (denominator == 0 ? 1.0 : numerator / denominator);

        for (int i = 0; i < dim; i++) {
            secantHessian[i][i] = fraction;
        }
    }

    private void updateSecantHessian(Secant lastSecant, Secant currentSecant) {
            /*
            Implementation of formula 11.17 on Page 281 of Ken Lange's Optimization book.
             */
        double[] d = new double[dim];
        double[] g = new double[dim];
        double[] Hd = new double[dim];
        double b = 0.0;
        double c = 0.0;

        for (int i = 0; i < dim; i++) {
            d[i] = currentSecant.getPosition(i) - lastSecant.getPosition(i);
            g[i] = currentSecant.getGradient(i) - lastSecant.getGradient(i);
            b += d[i] * g[i];
        }
        b = 1.0 / b;

        for (int i = 0; i < dim; i++) {
            double sum = 0.0;
            for (int j =0; j < dim; j++) {
                sum += secantHessian[i][j] * d[j];
            }
            Hd[i] = sum;
            c += d[i] * Hd[i];
        }
        c = -1.0 / c;

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                secantHessian[i][j] += b * g[i] * g[j] + c * Hd[i] * Hd[j];
            }
        }
    }
}
