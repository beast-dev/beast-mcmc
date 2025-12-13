package dr.evomodel.continuous;

import dr.inference.model.BlockDiagonalCosSinMatrixParameter;
import dr.inference.model.MatrixParameterInterface;

/**
 * Strategy for block-diagonal strength-of-selection matrices:
 * A = R D R^{-1}, D given in compressed form [diag|upper|lower].
 */
final class BlockDiagonalStrategy implements MultivariateElasticModel.SelectionMatrixStrategy {

    @Override
    public BasisRepresentation computeBasis(MatrixParameterInterface param, int dim) {
        if (!(param instanceof BlockDiagonalCosSinMatrixParameter)) {
            throw new IllegalArgumentException(
                    "BlockDiagonalStrategy requires BlockDiagonalCosSinMatrixParameter");
        }

        BlockDiagonalCosSinMatrixParameter blockParam = (BlockDiagonalCosSinMatrixParameter) param;

        // D compressed: [diag (dim), upper (dim-1), lower (dim-1)]
        double[] valuesD = new double[blockParam.getCompressedDDimension()];

        // R and Rinv row-major
        double[] valuesR = new double[dim*dim];
        double[] valuesRinv = new double[dim*dim];

        blockParam.fillBlockDiagonalElements(valuesD); //TODO this is a useless allocation, fix later
        blockParam.fillRAndRinv(valuesR, valuesRinv); //TODO this is a useless allocation, fix later

        BlockStructure bs = new BlockStructure(
                blockParam.getBlockStarts(),
                blockParam.getBlockSizes()
        );
        int nParameters = blockParam.getNumberOfParameters();

        return new BasisRepresentation(dim, nParameters, valuesD, valuesR, valuesRinv, bs);
    }


    @Override public boolean isDiagonal() { return false; }
    @Override public boolean isSymmetric() { return false; }
    @Override public boolean isBlockDiagonal() { return true; }
}
