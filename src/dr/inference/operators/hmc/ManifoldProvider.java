package dr.inference.operators.hmc;

import dr.inference.model.Parameter;
import dr.math.geodesics.Manifold;
import dr.math.matrixAlgebra.EJMLUtils;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.math.matrixAlgebra.WritableVector;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.ArrayList;


public interface ManifoldProvider {

//    public final ArrayList<Manifold> manifolds;

//    protected ManifoldProvider(ArrayList<Manifold> manifolds) {
//        this.manifolds = manifolds;
//    }
//
//
//    public int getManifoldCount() {
//        return manifolds.size();
//    }

//    public Manifold getManifold(int i) {
//        return manifolds.get(i);
//    }

//    public abstract double[] extractManifoldData(int i, double[] data);
//
//    public double[] extractManifoldData(int i, WrappedVector data) {
//        throw new RuntimeException("Not yet implemented");
//    }
//
//    public abstract void injectManifoldData(int i, double[] manifoldData, double[] sinkData);
//
//    public void injectManifoldData(int i, double[] manifoldData, WrappedVector sinkData) {
//        throw new RuntimeException("Not yet implemented");
//    }

    void projectTangent(double[] momentum, double[] position);

    void updatePositionAndMomentum(double[] position, WrappedVector momentum, double functionalStepSize);

    int getDimension();


    public class BlockManifoldProvider implements ManifoldProvider {

        private final ArrayList<ManifoldProvider> providers;
        private final double[][] positionBuffers;
        private final double[][] momentumBuffers;
//
//    private final int[] starts;
//    private final int dim;

        public BlockManifoldProvider(ArrayList<ManifoldProvider> manifoldProviders) {
            this.providers = manifoldProviders;
//        this.manifolds = manifolds;
//        this.starts = starts;
//        this.dim = dim;
//        int n = starts.length;
//        if (n != manifolds.size()) {
//            throw new RuntimeException("Bad structure");
//        }

            int n = manifoldProviders.size();
//
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

//                for (int j = 0; j < momentumBuffers[i].length; j++) {
//                    momentumBuffers[i][j] = momentum.get(j + start);
//                }
                WrappedVector momentumI = new WrappedVector.View(momentum, start, dim);

                System.arraycopy(position, start, positionBuffers[i], 0, positionBuffers[i].length);

                provider.updatePositionAndMomentum(positionBuffers[i], momentumI, functionalStepSize);
//                manifolds.get(i).geodesic(positionBuffers[i], momentumBuffers[i], functionalStepSize);


//                for (int j = 0; j < momentumBuffers[i].length; j++) {
//                    momentum.set(start + j, momentumBuffers[i][j]);
//                }
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


//    public class BlockStiefelManifoldProvider implements ManifoldProvider { //TODO: make this just a block manifold provider, all the Stiefel manifold code should be in a separate class
//
//        private final int rowDim;
//        private final int colDim;
//
//        private final ArrayList<ArrayList<Integer>> orthogonalityStructure;
//        private final ArrayList<ArrayList<Integer>> orthogonalityBlockRows;
//
//        public BlockStiefelManifoldProvider(ArrayList<Manifold> manifolds,
//                                            int rowDim,
//                                            int colDim,
//                                            ArrayList<ArrayList<Integer>> orthogonalityStructure,
//                                            ArrayList<ArrayList<Integer>> orthogonalityBlockRows) {
//
//            super(manifolds);
//            this.rowDim = rowDim;
//            this.colDim = colDim;
//
//            this.orthogonalityStructure = orthogonalityStructure;
//            this.orthogonalityBlockRows = orthogonalityBlockRows;
//
//        }
//
//
////        @Override
////        public double[] extractManifoldData(int i, double[] data) {
////            throw new RuntimeException("Not yet implemented");
//////            return new double[0];
////        }
////
////        @Override
////        public void injectManifoldData(int i, double[] manifoldData, double[] sinkData) {
////            throw new RuntimeException("Not yet implemented");
////        }
//
//        @Override
//        public void projectTangent(double[] momentum, double[] position) {
//            for (int block = 0; block < getManifoldCount(); block++) {
//                DenseMatrix64F positionMatrix = setOrthogonalSubMatrix(position, block);
//                DenseMatrix64F momentumMatrix = setOrthogonalSubMatrix(momentum, block);
//
//                int nCols = orthogonalityStructure.get(block).size();
//                int nRows = orthogonalityBlockRows.get(block).size();
////            positionMatrix.setData(position);
////            momentumMatrix.setData(momentum);
//
//                DenseMatrix64F innerProduct = new DenseMatrix64F(nCols, nCols);
//
//                CommonOps.multTransB(positionMatrix, momentumMatrix, innerProduct);
//                EJMLUtils.addWithTransposed(innerProduct);
//
//                DenseMatrix64F projection = new DenseMatrix64F(nCols, nRows);
//
//                CommonOps.mult(0.5, innerProduct, positionMatrix, projection);
//                CommonOps.subtractEquals(momentumMatrix, projection);
//
//                unwrapSubMatrix(momentumMatrix, block, momentum);
//            }
//        }
//
//        @Override
//        public void updatePositionAndMomentum(double[] position, WrappedVector momentum,
//                                              double functionalStepSize) {
//
//            for (int block = 0; block < getManifoldCount(); block++) {
//
////                int nCols = orthogonalityStructure.get(block).size();
////                int nRows = orthogonalityBlockRows.get(block).size();
//
////            positionMatrix.setData(position);
////                DenseMatrix64F positionMatrix = setOrthogonalSubMatrix(position, block);
////                DenseMatrix64F momentumMatrix = setOrthogonalSubMatrix(momentum.getBuffer(), momentum.getOffset(), block);
////            System.arraycopy(momentum.getBuffer(), momentum.getOffset(), momentumMatrix.data, 0, momentum.getDim());
//
//                double[] blockPosition = setOrthogonalSubMatrix(position, block).getData(); // TODO: un-hack the code so it doesn't constantly convert between double[] and DenseMatrix64F
//                double[] blockMomentum = setOrthogonalSubMatrix(momentum.getBuffer(), momentum.getOffset(), block).getData(); // TODO: un-hack the code so it doesn't constantly convert between double[] and DenseMatrix64F
//
//                Manifold manifold = getManifold(block);
//
//                manifold.geodesic(blockPosition, blockMomentum, functionalStepSize);
//
////                injectManifoldData(block, blockPosition, position);
////                injectManifoldData(block, blockMomentum, momentum);
//
//                int colDimi = orthogonalityStructure.get(block).size();
//                int rowDimi = orthogonalityBlockRows.get(block).size();
//                unwrapSubMatrix(DenseMatrix64F.wrap(colDimi, rowDimi, blockPosition), block, position); // TODO: un-hack the code so it doesn't constantly convert between double[] and DenseMatrix64F
//                unwrapSubMatrix(DenseMatrix64F.wrap(colDimi, rowDimi, blockMomentum), block, momentum.getBuffer(), momentum.getOffset()); // TODO: un-hack the code so it doesn't constantly convert between double[] and DenseMatrix64F
////                unwrapSubMatrix(momentumMatrix, block, momentum.getBuffer(), momentum.getOffset());
////            System.arraycopy(positionMatrix.data, 0, position, 0, position.length);
////            System.arraycopy(momentumMatrix.data, 0, momentum.getBuffer(), momentum.getOffset(), momentum.getDim());
//            }
//        }
//
//
////            ArrayList<Integer> subColList = new ArrayList<>();
////            for (int i : subColumns) {
////                subColList.add(i);
////            }
////
////            //check that orthogonalityStructure is consistent with the subRows
////            ArrayList<Integer> alreadyOrthogonal = new ArrayList<>();
////
////            for (int i = 0; i < newOrthogonalColumns.size(); i++) {
////                for (int j = 0; j < newOrthogonalColumns.get(i).length; j++) {
////                    if (!subColList.contains(newOrthogonalColumns.get(i)[j])) { //TODO: check that we're doing this by row (or allow to do by row or column)
////                        throw new RuntimeException("Cannot enforce orthogonality structure.");
////                    }
////                    if (alreadyOrthogonal.contains(newOrthogonalColumns.get(i)[j])) {
////                        throw new RuntimeException("Orthogonal blocks must be non-overlapping");
////                    }
////                    alreadyOrthogonal.add(newOrthogonalColumns.get(i)[j]);
////                }
////                orthogonalityStructure.add(newOrthogonalColumns.get(i));
////            }
////
////            for (int i = 0; i < subColumns.length; i++) {
////                if (!alreadyOrthogonal.contains(subColumns[i])) {
////                    orthogonalityStructure.add(new int[]{subColumns[i]});
////                }
////            }
//
//
//        private DenseMatrix64F setOrthogonalSubMatrix(double[] src, int srcOffset, int block) {
//
//            ArrayList<Integer> blockCols = orthogonalityStructure.get(block);
//            ArrayList<Integer> blockRows = orthogonalityBlockRows.get(block);
//            int nCols = blockCols.size();
//            int nRows = blockRows.size();
//
//            DenseMatrix64F dest = new DenseMatrix64F(nCols, nRows);
//
//            for (int row = 0; row < nRows; row++) {
//                for (int col = 0; col < nCols; col++) {
//                    int ind = rowDim * blockCols.get(col) + blockRows.get(row) + srcOffset;
//                    dest.set(col, row, src[ind]);
//                }
//            }
//
//            return dest;
//        }
//
//        private DenseMatrix64F setOrthogonalSubMatrix(double[] src, int block) {
//            return setOrthogonalSubMatrix(src, 0, block);
//        }
//
//        private void unwrapSubMatrix(DenseMatrix64F src, int block, double[] dest, int destOffset) {
//            int nRowsOriginal = rowDim;
//            ArrayList<Integer> blockCols = orthogonalityStructure.get(block);
//            ArrayList<Integer> blockRows = orthogonalityBlockRows.get(block);
//
//            for (int row = 0; row < blockRows.size(); row++) {
//                for (int col = 0; col < blockCols.size(); col++) {
//                    int ind = nRowsOriginal * blockCols.get(col) + blockRows.get(row) + destOffset;
//                    dest[ind] = src.get(col, row);
//                }
//            }
//        }
//
//        private void unwrapSubMatrix(DenseMatrix64F src, int block, double[] dest) {
//            unwrapSubMatrix(src, block, dest, 0);
//        }
//
//    }

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

//        @Override
//        public double[] extractManifoldData(int i, double[] data) {
//            throw new RuntimeException("not yet implemented");
//        }
//
//        @Override
//        public void injectManifoldData(int i, double[] manifoldData, double[] sinkData) {
//            throw new RuntimeException("not yet implemented");
//        }

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


