package dr.inference.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inferencexml.operators.hmc.GeodesicHamiltonianMonteCarloOperatorParser;
import dr.math.matrixAlgebra.EJMLUtils;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SkewSymmetricMatrixExponential;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;
import dr.xml.Reportable;
import org.ejml.alg.dense.decomposition.TriangularSolver;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.CholeskyDecomposition;
import org.ejml.ops.CommonOps;

import java.util.ArrayList;
import java.util.Collections;

public class GeodesicHamiltonianMonteCarloOperator extends HamiltonianMonteCarloOperator implements Reportable {

    public GeodesicHamiltonianMonteCarloOperator(AdaptationMode mode, double weight,
                                                 GradientWrtParameterProvider gradientProvider, Parameter parameter,
                                                 Transform transform, Parameter maskParameter, Options runtimeOptions,
                                                 MassPreconditioner preconditioner) {
        super(mode, weight, gradientProvider, parameter, transform, maskParameter, runtimeOptions, preconditioner);
        this.leapFrogEngine = new GeodesicLeapFrogEngine(parameter, getDefaultInstabilityHandler(), preconditioning, mask);
    }

    @Override
    public String getOperatorName() {
        return "GeodesicHMC(" + parameter.getParameterName() + ")";
    }

    @Override
    public String getReport() {

        MatrixParameterInterface matParam = (MatrixParameterInterface) parameter;
        int k = matParam.getColumnDimension();
        int p = matParam.getRowDimension();

        StringBuilder sb = new StringBuilder("operator: " + GeodesicHamiltonianMonteCarloOperatorParser.OPERATOR_NAME);
        sb.append("\n");
        sb.append("\toriginal position:\n");
        Matrix originalPosition = new Matrix(matParam.getParameterAsMatrix());
        sb.append(originalPosition.toString(2));

        double[] momentum = new double[parameter.getDimension()];

        for (int i = 0; i < momentum.length; i++) { //Need some deterministic way to assign momentum for test
            momentum[i] = i;
        }

        Matrix originalMomentum = new Matrix(p, k);
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < k; j++) {
                originalMomentum.set(i, j, momentum[i + j * p]);
            }
        }
        sb.append("\toriginal momentum (unprojected):\n");
        sb.append(originalMomentum.toString(2));

        WrappedVector wrappedMomentum = new WrappedVector.Raw(momentum);
        double hastings;
        try {
            hastings = leapFrogGivenMomentum(wrappedMomentum);
        } catch (NumericInstabilityException e) {
            e.printStackTrace();
            throw new RuntimeException("HMC failed");
        }

        Matrix finalPosition = new Matrix(matParam.getParameterAsMatrix());
        sb.append("\n");
        sb.append("\tfinal position:\n");
        sb.append(finalPosition.toString(2));
        sb.append("\n");
        sb.append("\thastings ratio: " + hastings + "\n\n");


        return sb.toString();
    }

    public void setOrthogonalityStructure(ArrayList<ArrayList<Integer>> oldOrthogonalityStructure) {
        ((GeodesicLeapFrogEngine) leapFrogEngine).setOrthogonalityStructure(oldOrthogonalityStructure);
    }


    public static class GeodesicLeapFrogEngine extends HamiltonianMonteCarloOperator.LeapFrogEngine.Default {

        private final MatrixParameterInterface matrixParameter;
        //        private final DenseMatrix64F positionMatrix;
//        private final DenseMatrix64F innerProduct;
//        private final DenseMatrix64F innerProduct2;
//        private final DenseMatrix64F projection;
//        private final DenseMatrix64F momentumMatrix;
//        private final int nRows;
//        private final int nCols;

        //        private final int[] subRows;
//        private final int[] subColumns;
        private final ArrayList<ArrayList<Integer>> orthogonalityStructure;
        private final ArrayList<ArrayList<Integer>> orthogonalityBlockRows;


        GeodesicLeapFrogEngine(Parameter parameter, HamiltonianMonteCarloOperator.InstabilityHandler instabilityHandler,
                               MassPreconditioner preconditioning, double[] mask) {
            super(parameter, instabilityHandler, preconditioning, mask);
            this.matrixParameter = (MatrixParameterInterface) parameter;


//            this.subRows = parseSubRowsFromMask();
//            this.subColumns = parseSubColumnsFromMask();


            this.orthogonalityStructure = new ArrayList<>();
            this.orthogonalityBlockRows = new ArrayList<>();

            if (mask == null) {
                ArrayList<Integer> rows = new ArrayList<>();
                for (int i = 0; i < matrixParameter.getRowDimension(); i++) {
                    rows.add(i);
                }
                ArrayList<Integer> cols = new ArrayList<>();
                for (int i = 0; i < matrixParameter.getColumnDimension(); i++) {
                    cols.add(i);
                }
                orthogonalityStructure.add(cols);
                orthogonalityBlockRows.add(rows);
            } else {
                parseStructureFromMask(mask);
            }
//            orthogonalityStructure.add(subColumns);

//            this.nRows = subRows.length;

//            if (mask != null) checkMask(subRows, subColumns);
//            this.nCols = subColumns.length;
//            this.positionMatrix = new DenseMatrix64F(nCols, nRows);
//            this.innerProduct = new DenseMatrix64F(nCols, nCols);
//            this.innerProduct2 = new DenseMatrix64F(nCols, nCols);
//            this.projection = new DenseMatrix64F(nCols, nRows);
//            this.momentumMatrix = new DenseMatrix64F(nCols, nRows);
        }

        private void parseStructureFromMask(double[] mask) {
            int nRows = matrixParameter.getRowDimension();
            int nCols = matrixParameter.getColumnDimension();

            ArrayList<Integer> colRows = new ArrayList<>();

            for (int i = 0; i < nCols; i++) {
                colRows.clear();
                int offset = i * nRows;
                for (int j = 0; j < nRows; j++) {
                    if (mask[offset + j] == 1) {
                        colRows.add(j);
                    }
                }

                if (!colRows.isEmpty()) {
                    int matchingInd = findMatchingArray(orthogonalityBlockRows, colRows);
                    if (matchingInd == -1) {
                        ArrayList<Integer> newBlock = new ArrayList<>();
                        newBlock.add(i);
                        orthogonalityStructure.add(newBlock);
                        orthogonalityBlockRows.add(new ArrayList<>(colRows));
                    } else {
                        orthogonalityStructure.get(matchingInd).add(i);
                    }
                }
            }
        }

        private int findMatchingArray(ArrayList<ArrayList<Integer>> listOfLists, ArrayList<Integer> list) {
            int nLists = listOfLists.size();
            for (int i = 0; i < nLists; i++) {
                ArrayList<Integer> subList = listOfLists.get(i);
                boolean matching = true;
                if (list.size() == subList.size()) {
                    for (int j = 0; j < list.size(); j++) {
                        if (list.get(j) != subList.get(j)) {
                            matching = false;
                            break;
                        }
                    }

                    if (matching) {
                        return i;
                    }
                }
            }

            return -1;
        }

        private int findSubArray(ArrayList<ArrayList<Integer>> listOfLists, ArrayList<Integer> list, ArrayList<Integer> remainingList) { //assumes both are sorted
            int nLists = listOfLists.size();
            for (int i = 0; i < nLists; i++) {
                ArrayList<Integer> subList = listOfLists.get(i);
                remainingList.clear();

                if (list.size() <= subList.size()) {
                    int currentInd = 0;
                    for (int j = 0; j < subList.size(); j++) {

                        if (currentInd < list.size() && subList.get(j) == list.get(currentInd)) {
                            currentInd += 1;
                        } else {
                            remainingList.add(subList.get(j));
                        }
                    }

                    if (currentInd == list.size()) {
                        return i;
                    }

                }
            }

            return -1;
        }


        public void setOrthogonalityStructure(ArrayList<ArrayList<Integer>> newOrthogonalColumns) {

            for (int i = 0; i < newOrthogonalColumns.size(); i++) {
                ArrayList<Integer> remainingList = new ArrayList<>();
                ArrayList<Integer> cols = newOrthogonalColumns.get(i);
                Collections.sort(cols);
                int matchingCol = findSubArray(orthogonalityStructure, cols, remainingList);
                if (matchingCol == -1) {
                    throw new RuntimeException("Orthogonality structure incompatible with itself or mask.");
                }

                ArrayList<Integer> existingCols = orthogonalityStructure.get(matchingCol);


                if (remainingList.size() > 0) {
                    orthogonalityStructure.set(matchingCol, remainingList);
                    orthogonalityStructure.add(cols);
                    orthogonalityBlockRows.add(orthogonalityBlockRows.get(matchingCol));

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

        }

//        private int[] parseSubColumnsFromMask() {
//
//            int originalRows = matrixParameter.getRowDimension();
//            int originalColumns = matrixParameter.getColumnDimension();
//
//            ArrayList<Integer> subArray = new ArrayList<Integer>();
//
//            for (int col = 0; col < originalColumns; col++) {
//                int offset = col * originalRows;
//                for (int row = 0; row < originalRows; row++) {
//                    int ind = offset + row;
//                    if (mask == null || mask[ind] == 1.0) {
//                        subArray.add(col);
//                        break;
//                    }
//                }
//            }
//
//            int[] subColumns = new int[subArray.size()];
//            for (int i = 0; i < subColumns.length; i++) {
//                subColumns[i] = subArray.get(i);
//            }
//
//            return subColumns;
//        }
//
//        private int[] parseSubRowsFromMask() {
//            int originalRows = matrixParameter.getRowDimension();
//            int originalColumns = matrixParameter.getColumnDimension();
//
//            ArrayList<Integer> subArray = new ArrayList<Integer>();
//
//            for (int row = 0; row < originalRows; row++) {
//                for (int col = 0; col < originalColumns; col++) {
//                    int ind = col * originalRows + row;
//                    if (mask == null || mask[ind] == 1.0) {
//                        subArray.add(row);
//                        break;
//                    }
//                }
//            }
//
//            int[] subRows = new int[subArray.size()];
//            for (int i = 0; i < subRows.length; i++) {
//                subRows[i] = subArray.get(i);
//            }
//
//            return subRows;
//        }

//        private void checkMask(int[] rows, int[] cols) {
//            int originalRows = matrixParameter.getRowDimension();
//            int originalColumns = matrixParameter.getColumnDimension();
//
//            int subRowInd = 0;
//            int subColInd = 0;
//
//            Boolean isSubRow;
//            Boolean isSubCol;
//
//            for (int row = 0; row < originalRows; row++) {
//                if (subRowInd < rows.length && row == rows[subRowInd]) {
//                    isSubRow = true;
//                    subRowInd++;
//                } else {
//                    isSubRow = false;
//                }
//
//                subColInd = 0;
//
//                for (int col = 0; col < originalColumns; col++) {
//                    if (subColInd < cols.length && col == cols[subColInd]) {
//                        isSubCol = true;
//                        subColInd++;
//                    } else {
//                        isSubCol = false;
//                    }
//
//                    int ind = originalRows * col + row;
//
//                    if (isSubCol && isSubRow) {
//                        if (mask[ind] != 1.0) {
//                            throw new RuntimeException("mask is incompatible with " +
//                                    GeodesicHamiltonianMonteCarloOperatorParser.OPERATOR_NAME +
//                                    ". All elements in sub-matrix must be set to 1.");
//                        }
//                    } else {
//                        if (mask[ind] != 0.0) {
//                            throw new RuntimeException("mask is incompatible with " +
//                                    GeodesicHamiltonianMonteCarloOperatorParser.OPERATOR_NAME +
//                                    ". All elements outside of sub-matrix must be set to 0.");
//                        }
//                    }
//
//                }
//            }
//        }

        private DenseMatrix64F setOrthogonalSubMatrix(double[] src, int srcOffset, int block) {

            int nRowsOriginal = matrixParameter.getRowDimension();
            ArrayList<Integer> blockCols = orthogonalityStructure.get(block);
            ArrayList<Integer> blockRows = orthogonalityBlockRows.get(block);
            int nCols = blockCols.size();
            int nRows = blockRows.size();

            DenseMatrix64F dest = new DenseMatrix64F(nCols, nRows);

            for (int row = 0; row < nRows; row++) {
                for (int col = 0; col < nCols; col++) {
                    int ind = nRowsOriginal * blockCols.get(col) + blockRows.get(row) + srcOffset;
                    dest.set(col, row, src[ind]);
                }
            }

            return dest;
        }

        private DenseMatrix64F setOrthogonalSubMatrix(double[] src, int block) {
            return setOrthogonalSubMatrix(src, 0, block);
        }

        private void unwrapSubMatrix(DenseMatrix64F src, int block, double[] dest, int destOffset) {
            int nRowsOriginal = matrixParameter.getRowDimension();
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

        @Override
        public void updateMomentum(double[] position, double[] momentum, double[] gradient,
                                   double functionalStepSize) throws HamiltonianMonteCarloOperator.NumericInstabilityException {
            super.updateMomentum(position, momentum, gradient, functionalStepSize);
            projectMomentum(momentum, position);

        }

        @Override
        public void updatePosition(double[] position, WrappedVector momentum,
                                   double functionalStepSize) throws HamiltonianMonteCarloOperator.NumericInstabilityException {

            for (int block = 0; block < orthogonalityStructure.size(); block++) {

                int nCols = orthogonalityStructure.get(block).size();
                int nRows = orthogonalityBlockRows.get(block).size();

//            positionMatrix.setData(position);
                DenseMatrix64F positionMatrix = setOrthogonalSubMatrix(position, block);
                DenseMatrix64F momentumMatrix = setOrthogonalSubMatrix(momentum.getBuffer(), momentum.getOffset(), block);
//            System.arraycopy(momentum.getBuffer(), momentum.getOffset(), momentumMatrix.data, 0, momentum.getDim());

                DenseMatrix64F innerProduct = new DenseMatrix64F(nCols, nCols);
                DenseMatrix64F innerProduct2 = new DenseMatrix64F(nCols, nCols);

                CommonOps.multTransB(positionMatrix, momentumMatrix, innerProduct);
                CommonOps.multTransB(momentumMatrix, momentumMatrix, innerProduct2);

                double[][] VtV = new double[2 * nCols][2 * nCols];

                for (int i = 0; i < nCols; i++) {
                    VtV[i + nCols][i] = 1;
                    for (int j = 0; j < nCols; j++) {
                        VtV[i][j] = innerProduct.get(i, j);
                        VtV[i + nCols][j + nCols] = innerProduct.get(i, j);
                        VtV[i][j + nCols] = -innerProduct2.get(j, i);
                    }
                }

                double[] expBuffer = new double[nCols * nCols];
                CommonOps.scale(-functionalStepSize, innerProduct);
                SkewSymmetricMatrixExponential matExp1 = new SkewSymmetricMatrixExponential(nCols);
                matExp1.exponentiate(innerProduct.data, expBuffer);

                double[] expBuffer2 = new double[nCols * nCols * 4];
                SkewSymmetricMatrixExponential matExp2 = new SkewSymmetricMatrixExponential(nCols * 2); //TODO: better matrix exponential
                DenseMatrix64F VtVmat = new DenseMatrix64F(VtV);
                CommonOps.scale(functionalStepSize, VtVmat);
                matExp2.exponentiate(VtVmat.data, expBuffer2);

                DenseMatrix64F X = new DenseMatrix64F(nCols * 2, nCols * 2);
                DenseMatrix64F Y = new DenseMatrix64F(nCols * 2, nCols * 2);

                for (int i = 0; i < nCols; i++) {
                    for (int j = 0; j < nCols; j++) {
                        X.set(i, j, expBuffer[i * nCols + j]);
                        X.set(i + nCols, j + nCols, expBuffer[i * nCols + j]);
                    }
                }
                Y.setData(expBuffer2);

                DenseMatrix64F Z = new DenseMatrix64F(nCols * 2, nCols * 2);

                CommonOps.mult(Y, X, Z);

                DenseMatrix64F PM = new DenseMatrix64F(nCols * 2, nRows);
                for (int i = 0; i < nRows; i++) {
                    for (int j = 0; j < nCols; j++) {
                        PM.set(j, i, positionMatrix.get(j, i));
                        PM.set(j + nCols, i, momentumMatrix.get(j, i));
                    }
                }

                DenseMatrix64F W = new DenseMatrix64F(2 * nCols, nRows);
                CommonOps.transpose(Z);
                CommonOps.mult(Z, PM, W);

                for (int i = 0; i < nRows; i++) {
                    for (int j = 0; j < nCols; j++) {
                        positionMatrix.set(j, i, W.get(j, i));
                        momentumMatrix.set(j, i, W.get(j + nCols, i));
                    }
                }

                //TODO: only run chunk below occasionally
                CommonOps.multTransB(positionMatrix, positionMatrix, innerProduct);
                CholeskyDecomposition cholesky = DecompositionFactory.chol(nCols, true);
                cholesky.decompose(innerProduct);
                TriangularSolver.invertLower(innerProduct.data, nCols);

                DenseMatrix64F projection = new DenseMatrix64F(nCols, nRows);

                CommonOps.mult(innerProduct, positionMatrix, projection);

                double sse = 0;
                for (int i = 0; i < positionMatrix.data.length; i++) {
                    double diff = projection.data[i] - positionMatrix.data[i];
                    sse += diff * diff;
                }

                if (sse / position.length > 1e-2) { //TODO: actually figure out if I want this
                    System.err.println("unstable"); //TODO: REMOVE
                    throw new NumericInstabilityException();
                }


                System.arraycopy(projection.data, 0, positionMatrix.data, 0, positionMatrix.data.length);

                unwrapSubMatrix(positionMatrix, block, position);
                unwrapSubMatrix(momentumMatrix, block, momentum.getBuffer(), momentum.getOffset());
//            System.arraycopy(positionMatrix.data, 0, position, 0, position.length);
//            System.arraycopy(momentumMatrix.data, 0, momentum.getBuffer(), momentum.getOffset(), momentum.getDim());
            }

            matrixParameter.setAllParameterValuesQuietly(position, 0);
            matrixParameter.fireParameterChangedEvent();
        }

        @Override
        public void projectMomentum(double[] momentum, double[] position) {
            for (int block = 0; block < orthogonalityStructure.size(); block++) {
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
    }
}
