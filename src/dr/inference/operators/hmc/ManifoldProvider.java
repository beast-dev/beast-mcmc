package dr.inference.operators.hmc;

import dr.math.geodesics.Manifold;
import dr.math.matrixAlgebra.EJMLUtils;
import dr.math.matrixAlgebra.WrappedVector;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.ArrayList;
import java.util.Arrays;

public abstract class ManifoldProvider {

    public final ArrayList<Manifold> manifolds;

    protected ManifoldProvider(ArrayList<Manifold> manifolds) {
        this.manifolds = manifolds;
    }


    public int getManifoldCount() {
        return manifolds.size();
    }

    public Manifold getManifold(int i) {
        return manifolds.get(i);
    }

    public abstract double[] extractManifoldData(int i, double[] data);

    public double[] extractManifoldData(int i, WrappedVector data) {
        throw new RuntimeException("Not yet implemented");
    }

    public abstract void injectManifoldData(int i, double[] manifoldData, double[] sinkData);

    public void injectManifoldData(int i, double[] manifoldData, WrappedVector sinkData) {
        throw new RuntimeException("Not yet implemented");
    }

    public abstract void projectTangent(double[] momentum, double[] position);

    public abstract void updatePositionAndMomentum(double[] position, WrappedVector momentum, double functionalStepSize);


    public static class BlockStiefelManifoldProvider extends ManifoldProvider { //TODO: make this just a block manifold provider, all the Stiefel manifold code should be in a separate class

        private final int rowDim;
        private final int colDim;

        private final ArrayList<ArrayList<Integer>> orthogonalityStructure;
        private final ArrayList<ArrayList<Integer>> orthogonalityBlockRows;

        public BlockStiefelManifoldProvider(ArrayList<Manifold> manifolds,
                                            int rowDim,
                                            int colDim,
                                            ArrayList<ArrayList<Integer>> orthogonalityStructure,
                                            ArrayList<ArrayList<Integer>> orthogonalityBlockRows) {

            super(manifolds);
            this.rowDim = rowDim;
            this.colDim = colDim;

            this.orthogonalityStructure = orthogonalityStructure;
            this.orthogonalityBlockRows = orthogonalityBlockRows;

        }


        @Override
        public double[] extractManifoldData(int i, double[] data) {
            throw new RuntimeException("Not yet implemented");
//            return new double[0];
        }

        @Override
        public void injectManifoldData(int i, double[] manifoldData, double[] sinkData) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public void projectTangent(double[] momentum, double[] position) {
            for (int block = 0; block < getManifoldCount(); block++) {
                DenseMatrix64F positionMatrix = setOrthogonalSubMatrix(position, block);
                DenseMatrix64F momentumMatrix = setOrthogonalSubMatrix(momentum, block);

                int nCols = orthogonalityStructure.get(block).size();
                int nRows = orthogonalityBlockRows.get(block).size();
//            positionMatrix.setData(position);
//            momentumMatrix.setData(momentum);

                DenseMatrix64F innerProduct = new DenseMatrix64F(nCols, nCols);

                CommonOps.multTransB(positionMatrix, momentumMatrix, innerProduct);
                EJMLUtils.addWithTransposed(innerProduct);

                DenseMatrix64F projection = new DenseMatrix64F(nCols, nRows);

                CommonOps.mult(0.5, innerProduct, positionMatrix, projection);
                CommonOps.subtractEquals(momentumMatrix, projection);

                unwrapSubMatrix(momentumMatrix, block, momentum);
            }
        }

        @Override
        public void updatePositionAndMomentum(double[] position, WrappedVector momentum,
                                              double functionalStepSize) {

            for (int block = 0; block < getManifoldCount(); block++) {

//                int nCols = orthogonalityStructure.get(block).size();
//                int nRows = orthogonalityBlockRows.get(block).size();

//            positionMatrix.setData(position);
//                DenseMatrix64F positionMatrix = setOrthogonalSubMatrix(position, block);
//                DenseMatrix64F momentumMatrix = setOrthogonalSubMatrix(momentum.getBuffer(), momentum.getOffset(), block);
//            System.arraycopy(momentum.getBuffer(), momentum.getOffset(), momentumMatrix.data, 0, momentum.getDim());

                double[] blockPosition = setOrthogonalSubMatrix(position, block).getData(); // TODO: un-hack the code so it doesn't constantly convert between double[] and DenseMatrix64F
                double[] blockMomentum = setOrthogonalSubMatrix(momentum.getBuffer(), momentum.getOffset(), block).getData(); // TODO: un-hack the code so it doesn't constantly convert between double[] and DenseMatrix64F

                Manifold manifold = getManifold(block);

                manifold.geodesic(blockPosition, blockMomentum, functionalStepSize);

//                injectManifoldData(block, blockPosition, position);
//                injectManifoldData(block, blockMomentum, momentum);

                int colDimi = orthogonalityStructure.get(block).size();
                int rowDimi = orthogonalityBlockRows.get(block).size();
                unwrapSubMatrix(DenseMatrix64F.wrap(colDimi, rowDimi, blockPosition), block, position); // TODO: un-hack the code so it doesn't constantly convert between double[] and DenseMatrix64F
                unwrapSubMatrix(DenseMatrix64F.wrap(colDimi, rowDimi, blockMomentum), block, momentum.getBuffer(), momentum.getOffset()); // TODO: un-hack the code so it doesn't constantly convert between double[] and DenseMatrix64F
//                unwrapSubMatrix(momentumMatrix, block, momentum.getBuffer(), momentum.getOffset());
//            System.arraycopy(positionMatrix.data, 0, position, 0, position.length);
//            System.arraycopy(momentumMatrix.data, 0, momentum.getBuffer(), momentum.getOffset(), momentum.getDim());
            }
        }


//            ArrayList<Integer> subColList = new ArrayList<>();
//            for (int i : subColumns) {
//                subColList.add(i);
//            }
//
//            //check that orthogonalityStructure is consistent with the subRows
//            ArrayList<Integer> alreadyOrthogonal = new ArrayList<>();
//
//            for (int i = 0; i < newOrthogonalColumns.size(); i++) {
//                for (int j = 0; j < newOrthogonalColumns.get(i).length; j++) {
//                    if (!subColList.contains(newOrthogonalColumns.get(i)[j])) { //TODO: check that we're doing this by row (or allow to do by row or column)
//                        throw new RuntimeException("Cannot enforce orthogonality structure.");
//                    }
//                    if (alreadyOrthogonal.contains(newOrthogonalColumns.get(i)[j])) {
//                        throw new RuntimeException("Orthogonal blocks must be non-overlapping");
//                    }
//                    alreadyOrthogonal.add(newOrthogonalColumns.get(i)[j]);
//                }
//                orthogonalityStructure.add(newOrthogonalColumns.get(i));
//            }
//
//            for (int i = 0; i < subColumns.length; i++) {
//                if (!alreadyOrthogonal.contains(subColumns[i])) {
//                    orthogonalityStructure.add(new int[]{subColumns[i]});
//                }
//            }


        private DenseMatrix64F setOrthogonalSubMatrix(double[] src, int srcOffset, int block) {

            ArrayList<Integer> blockCols = orthogonalityStructure.get(block);
            ArrayList<Integer> blockRows = orthogonalityBlockRows.get(block);
            int nCols = blockCols.size();
            int nRows = blockRows.size();

            DenseMatrix64F dest = new DenseMatrix64F(nCols, nRows);

            for (int row = 0; row < nRows; row++) {
                for (int col = 0; col < nCols; col++) {
                    int ind = rowDim * blockCols.get(col) + blockRows.get(row) + srcOffset;
                    dest.set(col, row, src[ind]);
                }
            }

            return dest;
        }

        private DenseMatrix64F setOrthogonalSubMatrix(double[] src, int block) {
            return setOrthogonalSubMatrix(src, 0, block);
        }

        private void unwrapSubMatrix(DenseMatrix64F src, int block, double[] dest, int destOffset) {
            int nRowsOriginal = rowDim;
            ArrayList<Integer> blockCols = orthogonalityStructure.get(block);
            ArrayList<Integer> blockRows = orthogonalityBlockRows.get(block);

            for (int row = 0; row < blockRows.size(); row++) {
                for (int col = 0; col < blockCols.size(); col++) {
                    int ind = nRowsOriginal * blockCols.get(col) + blockRows.get(row) + destOffset;
                    dest[ind] = src.get(col, row);
                }
            }
        }

        private void unwrapSubMatrix(DenseMatrix64F src, int block, double[] dest) {
            unwrapSubMatrix(src, block, dest, 0);
        }

    }

    public static class BasicManifoldProvider extends ManifoldProvider {

        private final Manifold manifold;
        private final double[] momentumBuffer;

        public BasicManifoldProvider(Manifold manifold, int dim) {
            super(new ArrayList(Arrays.asList(manifold)));
            this.manifold = manifold;
            this.momentumBuffer = new double[dim];

        }

        @Override
        public double[] extractManifoldData(int i, double[] data) {
            throw new RuntimeException("not yet implemented");
        }

        @Override
        public void injectManifoldData(int i, double[] manifoldData, double[] sinkData) {
            throw new RuntimeException("not yet implemented");
        }

        @Override
        public void projectTangent(double[] momentum, double[] position) {
            manifold.projectTangent(momentum, position);
        }

        @Override
        public void updatePositionAndMomentum(double[] position, WrappedVector momentum, double functionalStepSize) {
            int dim = position.length;
            for (int i = 0; i < position.length; i++) { //TODO: don't jump back and forth between WrappedVector and double[]
                momentumBuffer[i] = momentum.get(i);
            }
            manifold.geodesic(position, momentumBuffer, functionalStepSize);
            for (int i = 0; i < dim; i++) {
                momentum.set(i, momentumBuffer[i]);
            }
        }
    }


}

