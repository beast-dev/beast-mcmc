package dr.evomodel.treedatalikelihood.continuous;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.missingData.MissingOps.wrap;

public final class OULegacyReportHelper {

    private OULegacyReportHelper() {
    }

    public static double[][] getJointVarianceFull(final int dim,
                                                  final int ntaxa,
                                                  final double priorSampleSize,
                                                  final double[][] treeSharedLengths,
                                                  final double[][] traitVariance,
                                                  final double[] eigenValues,
                                                  final double[] eigenVectors) {
        final DenseMatrix64F v = wrap(eigenVectors, 0, dim, dim);
        final DenseMatrix64F vInv = new DenseMatrix64F(dim, dim);
        CommonOps.invert(v, vInv);

        final DenseMatrix64F transTraitVariance = new DenseMatrix64F(traitVariance);
        final DenseMatrix64F tmp = new DenseMatrix64F(dim, dim);
        CommonOps.mult(vInv, transTraitVariance, tmp);
        CommonOps.multTransB(tmp, vInv, transTraitVariance);

        final double[][] invEigVals = new double[dim][dim];
        for (int p = 0; p < dim; ++p) {
            for (int q = 0; q < dim; ++q) {
                invEigVals[p][q] = factorFunction(eigenValues[p] + eigenValues[q], 1.0);
            }
        }

        final DenseMatrix64F varTemp = new DenseMatrix64F(dim, dim);
        final double[][] jointVariance = new double[dim * ntaxa][dim * ntaxa];

        for (int i = 0; i < ntaxa; ++i) {
            for (int j = 0; j < ntaxa; ++j) {
                final double ti = treeSharedLengths[i][i];
                final double tj = treeSharedLengths[j][j];
                final double tij = treeSharedLengths[i][j];
                for (int p = 0; p < dim; ++p) {
                    for (int q = 0; q < dim; ++q) {
                        final double ep = eigenValues[p];
                        final double eq = eigenValues[q];
                        varTemp.set(
                                p,
                                q,
                                Math.exp(-ep * ti) * Math.exp(-eq * tj)
                                        * (invEigVals[p][q] * (Math.exp((ep + eq) * tij) - 1.0) + 1.0 / priorSampleSize)
                                        * transTraitVariance.get(p, q));
                    }
                }
                CommonOps.mult(v, varTemp, tmp);
                CommonOps.multTransB(tmp, v, varTemp);
                for (int p = 0; p < dim; ++p) {
                    for (int q = 0; q < dim; ++q) {
                        jointVariance[i * dim + p][j * dim + q] = varTemp.get(p, q);
                    }
                }
            }
        }
        return jointVariance;
    }

    public static double[][] getJointVarianceDiagonal(final int dim,
                                                      final int ntaxa,
                                                      final double priorSampleSize,
                                                      final double[][] treeSharedLengths,
                                                      final double[][] traitVariance,
                                                      final double[] eigenValues) {
        final DenseMatrix64F varTemp = new DenseMatrix64F(dim, dim);
        final double[][] jointVariance = new double[dim * ntaxa][dim * ntaxa];

        for (int i = 0; i < ntaxa; ++i) {
            for (int j = 0; j < ntaxa; ++j) {
                final double ti = treeSharedLengths[i][i];
                final double tj = treeSharedLengths[j][j];
                final double tij = treeSharedLengths[i][j];
                for (int p = 0; p < dim; ++p) {
                    for (int q = 0; q < dim; ++q) {
                        final double ep = eigenValues[p];
                        final double eq = eigenValues[q];
                        final double var;
                        if (ep + eq == 0.0) {
                            var = (tij + 1.0 / priorSampleSize) * traitVariance[p][q];
                        } else {
                            var = Math.exp(-ep * ti) * Math.exp(-eq * tj)
                                    * (Math.expm1((ep + eq) * tij) / (ep + eq) + 1.0 / priorSampleSize)
                                    * traitVariance[p][q];
                        }
                        varTemp.set(p, q, var);
                    }
                }
                for (int p = 0; p < dim; ++p) {
                    for (int q = 0; q < dim; ++q) {
                        jointVariance[i * dim + p][j * dim + q] = varTemp.get(p, q);
                    }
                }
            }
        }
        return jointVariance;
    }

    public static void getMeanTipVariancesFull(final int dim,
                                               final double priorSampleSize,
                                               final double[] treeLengths,
                                               final DenseMatrix64F traitVariance,
                                               final DenseMatrix64F varSum,
                                               final double[] eigenValues,
                                               final double[] eigenVectors) {
        final DenseMatrix64F v = wrap(eigenVectors, 0, dim, dim);
        final DenseMatrix64F vInv = new DenseMatrix64F(dim, dim);
        CommonOps.invert(v, vInv);

        final DenseMatrix64F transTraitVariance = new DenseMatrix64F(traitVariance);
        final DenseMatrix64F tmp = new DenseMatrix64F(dim, dim);
        CommonOps.mult(vInv, transTraitVariance, tmp);
        CommonOps.multTransB(tmp, vInv, transTraitVariance);

        getMeanTipVariancesDiagonal(dim, priorSampleSize, treeLengths, transTraitVariance, varSum, eigenValues);

        CommonOps.mult(v, varSum, tmp);
        CommonOps.multTransB(tmp, v, varSum);
    }

    public static void getMeanTipVariancesDiagonal(final int dim,
                                                   final double priorSampleSize,
                                                   final double[] treeLengths,
                                                   final DenseMatrix64F traitVariance,
                                                   final DenseMatrix64F varSum,
                                                   final double[] eigenValues) {
        for (int i = 0; i < treeLengths.length; ++i) {
            final double ti = treeLengths[i];
            for (int p = 0; p < dim; ++p) {
                final double ep = eigenValues[p];
                for (int q = 0; q < dim; ++q) {
                    final double eq = eigenValues[q];
                    final double var;
                    if (ep + eq == 0.0) {
                        var = (ti + 1.0 / priorSampleSize) * traitVariance.get(p, q);
                    } else {
                        var = Math.exp(-(ep + eq) * ti)
                                * (Math.expm1((ep + eq) * ti) / (ep + eq) + 1.0 / priorSampleSize)
                                * traitVariance.get(p, q);
                    }
                    varSum.set(p, q, varSum.get(p, q) + var);
                }
            }
        }
        CommonOps.scale(1.0 / treeLengths.length, varSum);
    }

    private static double factorFunction(final double x, final double l) {
        if (x == 0.0) {
            return l;
        }
        return -Math.expm1(-x * l) / x;
    }
}
