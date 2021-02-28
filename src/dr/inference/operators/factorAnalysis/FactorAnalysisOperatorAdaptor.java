package dr.inference.operators.factorAnalysis;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.LatentFactorModel;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;
import static dr.evomodelxml.treedatalikelihood.ContinuousDataLikelihoodParser.FACTOR_NAME;

/**
 * @author Marc A. Suchard
 * @author Gabriel Hassler
 */
public interface FactorAnalysisOperatorAdaptor {

    int getNumberOfTaxa();

    int getNumberOfTraits();

    int getNumberOfFactors();

    double getFactorValue(int factor, int taxon);

    double getDataValue(int trait, int taxon);

    double getLoadingsValue(int dim);

    double getColumnPrecision(int index);

    void setLoadingsForTraitQuietly(int trait, double[] value);

    void reflectLoadingsForFactor(int factor);

    void fireLoadingsChanged();

    void drawFactors();

    boolean isNotMissing(int trait, int taxon);

    Parameter[] getFactorDependentParameters();

    abstract class Abstract implements FactorAnalysisOperatorAdaptor, Reportable {

        private final MatrixParameterInterface loadings;

        Abstract(MatrixParameterInterface loadings) {
            this.loadings = loadings;
        }

        @Override
        public void setLoadingsForTraitQuietly(int trait, double[] value) {
            for (int j = 0; j < value.length; j++) {
                loadings.setParameterValueQuietly(trait, j, value[j]);
            }
        }

        @Override
        public void reflectLoadingsForFactor(int factor) {
            double pivot = loadings.getParameterValue(factor, factor);
            if (pivot < 0.0) {
                for (int trait = factor; trait < loadings.getRowDimension(); ++trait) {
                    loadings.setParameterValueQuietly(trait, factor,
                            -1.0 * loadings.getParameterValue(trait, factor));
                }
            }
        }

        @Override
        public void fireLoadingsChanged() {
            loadings.fireParameterChangedEvent();
        }

        @Override
        public double getLoadingsValue(int dim) {
            return loadings.getParameterValue(dim);
        }

        @Override
        public String getReport() {
            int repeats = 1000000;
            int nFac = getNumberOfFactors();
            int nTaxa = getNumberOfTaxa();
            int dim = nFac * nTaxa;

            double[] sums = new double[dim];
            double[][] sumSquares = new double[dim][dim];


            for (int i = 0; i < repeats; i++) {

                fireLoadingsChanged();
                drawFactors();
                for (int j = 0; j < nTaxa; j++) {
                    for (int k = 0; k < nFac; k++) {
                        double x = getFactorValue(k, j);
                        sums[k * nTaxa + j] += x;

                        for (int l = 0; l < nTaxa; l++) {
                            for (int m = 0; m < nFac; m++) {
                                double y = getFactorValue(m, l);
                                sumSquares[k * nTaxa + j][m * nTaxa + l] += x * y;
                            }
                        }
                    }
                }
            }


            double[] mean = new double[dim];
            double[][] cov = new double[dim][dim];
            for (int i = 0; i < dim; i++) {
                mean[i] = sums[i] / repeats;
                for (int j = 0; j < dim; j++) {
                    sumSquares[i][j] /= repeats;
                }
            }
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    cov[i][j] = sumSquares[i][j] - mean[i] * mean[j];
                }
            }

            StringBuilder sb = new StringBuilder();
            sb.append(this.getClass() + "Report:\n");
            sb.append("Factor mean:\n");
            sb.append(new Vector(mean));
            sb.append("\n\n");
            sb.append("Factor covariance:\n");
            sb.append(new Matrix(cov));
            sb.append("\n\n");
            return sb.toString();
        }
    }

    class SampledFactors extends Abstract {

        private final LatentFactorModel LFM;

        public SampledFactors(LatentFactorModel LFM) {
            super(LFM.getLoadings());
            this.LFM = LFM;
        }

        @Override
        public int getNumberOfTaxa() {
            return LFM.getFactors().getColumnDimension();
        }

        @Override
        public int getNumberOfTraits() {
            return LFM.getLoadings().getRowDimension();
        }

        @Override
        public int getNumberOfFactors() {

            assert (LFM.getFactors().getRowDimension() == LFM.getLoadings().getColumnDimension());

            return LFM.getFactors().getRowDimension();
        }

        @Override
        public double getFactorValue(int factor, int taxon) {
            return LFM.getFactors().getParameterValue(factor, taxon);
        }

        @Override
        public double getDataValue(int trait, int taxon) {

            assert (taxon < getNumberOfTaxa());
            assert (trait < getNumberOfTraits());

            assert (LFM.getScaledData().getRowDimension() == getNumberOfTraits());
            assert (LFM.getScaledData().getColumnDimension() == getNumberOfTaxa());

            return LFM.getScaledData().getParameterValue(trait, taxon);
        }

        @Override
        public double getColumnPrecision(int index) {
            return LFM.getColumnPrecision().getParameterValue(index, index);
        }

        @Override
        public void drawFactors() {
            // Do nothing
        }

        @Override
        public boolean isNotMissing(int trait, int taxon) {
            Parameter missing = LFM.getMissingIndicator();
            int index = taxon * getNumberOfTraits() + trait;

            return missing == null || missing.getParameterValue(index) != 1.0;
        }

        @Override
        public Parameter[] getFactorDependentParameters() {
            return new Parameter[]{LFM.getFactors()};
        }
    }

    class IntegratedFactors extends Abstract {

        private final IntegratedFactorAnalysisLikelihood factorLikelihood;
        private final TreeDataLikelihood treeLikelihood;

        private final Parameter precision;
        private final CompoundParameter data;

        private final TreeTrait factorTrait;
        private double[] factors;

        public IntegratedFactors(IntegratedFactorAnalysisLikelihood factorLikelihood,
                                 TreeDataLikelihood treeLikelihood) {
            super(factorLikelihood.getLoadings());
            this.factorLikelihood = factorLikelihood;
            this.treeLikelihood = treeLikelihood;

            this.precision = factorLikelihood.getPrecision();
            this.data = factorLikelihood.getParameter();

            factorTrait = treeLikelihood.getTreeTrait(factorLikelihood.getTipTraitName());

            assert (factorTrait != null);
        }

        @Override
        public int getNumberOfTaxa() {
            return factorLikelihood.getNumberOfTaxa();
        }

        @Override
        public int getNumberOfTraits() {
            return factorLikelihood.getNumberOfTraits();
        }

        @Override
        public int getNumberOfFactors() {
            return factorLikelihood.getNumberOfFactors();
        }

        @Override
        public double getFactorValue(int factor, int taxon) {
            return factors[taxon * getNumberOfFactors() + factor];
        }

        @Override
        public double getDataValue(int trait, int taxon) {
            return data.getParameterValue(trait, taxon);
        }

        @Override
        public double getColumnPrecision(int index) {
            return precision.getParameterValue(index);
        }

        @Override
        public void drawFactors() {
            factors = (double[]) factorTrait.getTrait(treeLikelihood.getTree(), null);

            if (DEBUG) {
                System.err.println("factors: " + new Vector(factors));
            }
        }

        @Override
        public boolean isNotMissing(int trait, int taxon) {
            int index = taxon * getNumberOfTraits() + trait;
            return !factorLikelihood.getDataMissingIndicators()[index];
        }

        @Override
        public Parameter[] getFactorDependentParameters() {
            return new Parameter[]{
                    factorLikelihood.getLoadings(),
                    factorLikelihood.getParameter(),
                    factorLikelihood.getPrecision()
            };
        }

        private static final boolean DEBUG = false;
    }
}
