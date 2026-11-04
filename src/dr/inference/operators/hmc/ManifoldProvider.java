package dr.inference.operators.hmc;

import dr.inference.model.Parameter;
import dr.math.geodesics.Manifold;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.math.matrixAlgebra.WritableVector;

import java.util.ArrayList;


public interface ManifoldProvider {


    void projectTangent(double[] momentum, double[] position);

    void updatePositionAndMomentum(double[] position, WrappedVector momentum, double functionalStepSize);

    int getDimension();


    public class BlockManifoldProvider implements ManifoldProvider {

        private final ArrayList<ManifoldProvider> providers;
        private final double[][] positionBuffers;
        private final double[][] momentumBuffers;


        public BlockManifoldProvider(ArrayList<ManifoldProvider> manifoldProviders) {
            this.providers = manifoldProviders;


            int n = manifoldProviders.size();

            positionBuffers = new double[n][];
            momentumBuffers = new double[n][];
            for (int i = 0; i < n; i++) {
                int ni = manifoldProviders.get(i).getDimension();
                positionBuffers[i] = new double[ni];
                momentumBuffers[i] = new double[ni];
            }

        }


        @Override
        public void projectTangent(double[] momentum, double[] position) {
            int start = 0;
            for (int i = 0; i < providers.size(); i++) {
                ManifoldProvider provider = providers.get(i);
                int dim = provider.getDimension();

                System.arraycopy(momentum, start, momentumBuffers[i], 0, dim);
                System.arraycopy(position, start, positionBuffers[i], 0, dim);

                provider.projectTangent(momentumBuffers[i], positionBuffers[i]);

                System.arraycopy(momentumBuffers[i], 0, momentum, start, dim);
                System.arraycopy(positionBuffers[i], 0, position, start, dim);
                start = start + dim;
            }

        }

        @Override
        public void updatePositionAndMomentum(double[] position, WrappedVector momentum, double functionalStepSize) {
            int start = 0;
            for (int i = 0; i < providers.size(); i++) {
                ManifoldProvider provider = providers.get(i);
                int dim = provider.getDimension();


                WrappedVector momentumI = new WrappedVector.View(momentum, start, dim);

                System.arraycopy(position, start, positionBuffers[i], 0, positionBuffers[i].length);

                provider.updatePositionAndMomentum(positionBuffers[i], momentumI, functionalStepSize);

                System.arraycopy(positionBuffers[i], 0, position, start, positionBuffers[i].length);
                start = start + dim;
            }

        }

        @Override
        public int getDimension() {
            int dim = 0;
            for (ManifoldProvider provider : providers) {
                dim += provider.getDimension();
            }
            return dim;
        }
    }


    public class BasicManifoldProvider implements ManifoldProvider {

        private final Manifold manifold;
        private final double[] momentumBuffer;
        private final double[] positionBuffer;
        private final Parameter mask;

        public BasicManifoldProvider(Manifold manifold, int dim, Parameter mask) {
            this.manifold = manifold;
            this.momentumBuffer = new double[dim];
            this.positionBuffer = new double[dim];
            this.mask = mask;

        }


        //TODO: refactor to WrappedVectors and avoid all this messy mask stuff
        private void copyMask(double[] src, double[] dest) {
            int destPos = 0;
            for (int srcPos = 0; srcPos < src.length; srcPos++) {
                double val = mask.getParameterValue(srcPos);
                if (val == 1) {
                    dest[destPos] = src[srcPos];
                    destPos++;
                }
            }
        }

        private void copyMask(ReadableVector src, double[] dest) {
            int destPos = 0;
            for (int srcPos = 0; srcPos < src.getDim(); srcPos++) {
                double val = mask.getParameterValue(srcPos);
                if (val == 1) {
                    dest[destPos] = src.get(srcPos);
                    destPos++;
                }
            }
        }

        private void copyMaskReverse(double[] src, double[] dest) {
            int srcPos = 0;
            for (int destPos = 0; destPos < dest.length; destPos++) {
                double val = mask.getParameterValue(destPos);
                if (val == 1) {
                    dest[destPos] = src[srcPos];
                    srcPos++;
                }
            }
        }

        private void copyMaskReverse(double[] src, WritableVector dest) {
            int srcPos = 0;
            for (int destPos = 0; destPos < dest.getDim(); destPos++) {
                double val = mask.getParameterValue(destPos);
                if (val == 1) {
                    dest.set(destPos, src[srcPos]);
                    srcPos++;
                }
            }
        }

        @Override
        public void projectTangent(double[] momentum, double[] position) {
            if (mask == null) {
                manifold.projectTangent(momentum, position);
            } else {
                copyMask(momentum, momentumBuffer);
                copyMask(position, positionBuffer);
                manifold.projectTangent(momentumBuffer, positionBuffer);
                copyMaskReverse(momentumBuffer, momentum);
                copyMaskReverse(positionBuffer, position);
            }
        }

        @Override
        public void updatePositionAndMomentum(double[] position, WrappedVector momentum, double functionalStepSize) {
            if (mask != null) { //TODO: don't to this if else stuff
                copyMask(position, positionBuffer);
                copyMask(momentum, momentumBuffer);
                manifold.geodesic(positionBuffer, momentumBuffer, functionalStepSize);
                copyMaskReverse(positionBuffer, position);
                copyMaskReverse(momentumBuffer, momentum);
            } else {


                int dim = position.length;
                for (int i = 0; i < momentumBuffer.length; i++) { //TODO: don't jump back and forth between WrappedVector and double[]
                    momentumBuffer[i] = momentum.get(i);
                }
                manifold.geodesic(position, momentumBuffer, functionalStepSize);
                for (int i = 0; i < dim; i++) {
                    momentum.set(i, momentumBuffer[i]);
                }
            }

        }

        @Override
        public int getDimension() {
            return momentumBuffer.length;
        }
    }
}


