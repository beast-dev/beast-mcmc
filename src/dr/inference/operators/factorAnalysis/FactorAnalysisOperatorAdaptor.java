package dr.inference.operators.factorAnalysis;

import dr.inference.model.LatentFactorModel;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public interface FactorAnalysisOperatorAdaptor {

    int getNumberOfTaxa();

    int getNumberOfTraits();

    int getNumberOfFactors();

    double getFactorValue(int factor, int taxon);

    double getDataValue(int trait, int taxon);

    double getColumnPrecision(int index);

    void setLoadingsForTraitQuietly(int trait, double[] value);

    void fireLoadingsChanged();
    
    boolean isNotMissing(int trait, int taxon);

    class Original implements  FactorAnalysisOperatorAdaptor {

        private final LatentFactorModel LFM;

        public Original(LatentFactorModel LFM) {
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
        public void setLoadingsForTraitQuietly(int trait, double[] value) {
            MatrixParameterInterface changing = LFM.getLoadings();
            for (int j = 0; j < value.length; j++) {
                changing.setParameterValueQuietly(trait, j, value[j]);
            }
        }

        @Override
        public void fireLoadingsChanged() {
            LFM.getLoadings().fireParameterChangedEvent();
        }

        @Override
        public double getColumnPrecision(int index) {
            return LFM.getColumnPrecision().getParameterValue(index, index);
        }

        @Override
        public boolean isNotMissing(int trait, int taxon) {
            Parameter missing = LFM.getMissingIndicator();
            int index = taxon * getNumberOfTraits() + trait;

            return missing == null || missing.getParameterValue(index) != 1.0;
        }
    }
}
