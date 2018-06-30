package dr.inference.operators.hmc;

import dr.inference.hmc.HessianWrtParameterProvider;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;

/**
 * @author Marc A. Suchard
 */
public interface MassPreconditioner {

    // It looks like my interface is slowly reverting to something close to what Xiang had previously

     WrappedVector drawInitialMomentum();

     double getKineticEnergy(ReadableVector momentum);

     double getVelocity(int index, ReadableVector momentum);

     class NoPreconditioning implements MassPreconditioner {

         final int dim;

         NoPreconditioning(int dim) { this.dim = dim; }

         @Override
         public WrappedVector drawInitialMomentum() {

             double[] momentum = new double[dim];

             for (int i = 0; i < dim; ++i) {
                 momentum[i] = MathUtils.nextGaussian();
             }

             return new WrappedVector.Raw(momentum);
         }

         @Override
         public double getKineticEnergy(ReadableVector momentum) {
             double total = 0.0;
             for (int i = 0; i < dim; i++) {
                 total += momentum.get(i) * momentum.get(i) / 2.0;
             }
             return total;
         }

         @Override
         public double getVelocity(int i, ReadableVector momentum) {
             return momentum.get(i);
         }
     }

     class DiagonalPreconditioning implements MassPreconditioner {

         final int dim;
         WrappedVector mass = null;
         WrappedVector massInverse = null;
         final HessianWrtParameterProvider hessianWrtParameterProvider;

         DiagonalPreconditioning(int dim, HessianWrtParameterProvider hessianWrtParameterProvider) {
             this.dim = dim;
             this.hessianWrtParameterProvider = hessianWrtParameterProvider;
         }

         @Override
         public WrappedVector drawInitialMomentum() {

             if (mass == null || massInverse == null) {
                 updateMass();
             }

             double[] momentum = new double[dim];

             for (int i = 0; i < dim; i++) {
                 momentum[i] = MathUtils.nextGaussian() * Math.sqrt(mass.get(i));
             }

             return new WrappedVector.Raw(momentum);
         }

         private void updateMass() {
             double[] diagonalHessian = hessianWrtParameterProvider.getDiagonalHessianLogDensity();
             double[] boundedMassInverse = boundMassInverse(diagonalHessian);
             for (int i = 0; i < dim; i++) {
                 massInverse.set(i, boundedMassInverse[i]);
                 mass.set(i, 1.0 / massInverse.get(i));
             }
         }

         private double[] boundMassInverse(double[] diagonalHessian) {

             double sum = 0.0;
             final double lowerBound = 1E-2; //TODO bad magic numbers
             final double upperBound = 1E2;
             double[] boundedMassInverse = new double[dim];

             for (int i = 0; i < dim; i++) {
                 boundedMassInverse[i] = -1.0 / diagonalHessian[i];
                 if (boundedMassInverse[i] < lowerBound) {
                     boundedMassInverse[i] = lowerBound;
                 } else if (boundedMassInverse[i] > upperBound) {
                     boundedMassInverse[i] = upperBound;
                 }
                 sum += 1.0 / boundedMassInverse[i];
             }
             final double mean = sum / dim;
             for (int i = 0; i < dim; i++) {
                 boundedMassInverse[i] = boundedMassInverse[i] * mean;
             }
             return boundedMassInverse;
         }

         @Override
         public double getKineticEnergy(ReadableVector momentum) {
             double total = 0.0;
             for (int i = 0; i < dim; i++) {
                 total += momentum.get(i) * momentum.get(i) / (2.0 * mass.get(i));
             }
             return total;
         }

         @Override
         public double getVelocity(int i, ReadableVector momentum) {
             return momentum.get(i) / mass.get(i);
         }
     }

     class FullPreconditioning implements MassPreconditioner {

         final int dim;
         WrappedVector mass;
         WrappedVector massInverse;
         HessianWrtParameterProvider hessianWrtParameterProvider;

         FullPreconditioning(int dim, HessianWrtParameterProvider hessianWrtParameterProvider) {
             this.dim = dim;
             this.hessianWrtParameterProvider = hessianWrtParameterProvider;
         }

         @Override
         public WrappedVector drawInitialMomentum() {
             return null;
         }

         @Override
         public double getKineticEnergy(ReadableVector momentum) {
             return 0.0;
         }

         @Override
         public double getVelocity(int i, ReadableVector momentum) {
             return 0.0; // M^{-1}_{i} momentum
         }
     }
}
