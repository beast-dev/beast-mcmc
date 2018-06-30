package dr.inference.operators.hmc;

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
             return 0.0;
         }

         @Override
         public double getVelocity(int i, ReadableVector momentum) {
             return momentum.get(i);
         }
     }

     class DiagonalPreconditioning implements MassPreconditioner {

         WrappedVector mass;

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
             return momentum.get(i) / mass.get(i);
         }
     }

     class FullPreconditioning implements MassPreconditioner {

         WrappedVector mass;

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
