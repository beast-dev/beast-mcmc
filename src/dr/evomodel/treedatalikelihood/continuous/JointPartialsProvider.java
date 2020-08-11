package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.CompoundParameter;

import java.util.List;

public class JointPartialsProvider implements ContinuousTraitPartialsProvider {

    private final String name;
    private final ContinuousTraitPartialsProvider[] providers;
    private final int traitDim;
    private final int dataDim;
    private static final PrecisionType precisionType = PrecisionType.FULL; //TODO: base on child precisionTypes

    public JointPartialsProvider(String name, ContinuousTraitPartialsProvider[] providers) {
        this.name = name;
        this.providers = providers;

        int traitDim = 0;
        int dataDim = 0;
        for (ContinuousTraitPartialsProvider provider : providers) {
            traitDim += provider.getTraitDimension();
            dataDim += provider.getDataDimension();
        }

        this.traitDim = traitDim;
        this.dataDim = dataDim;

    }


    @Override
    public boolean bufferTips() {
        return true; //TODO: not sure what this does, but it's set to `true` for every implementation
    }

    @Override
    public int getTraitCount() {
        return providers[0].getTraitCount(); //TODO: make sure all have same trait count in parser
    }

    @Override
    public int getTraitDimension() {
        return traitDim; //TODO
    }

    @Override
    public int getDataDimension() {
        return dataDim; //TODO: maybe throw error here? Used for model extension, mse stuff and it might be worth putting conditions if JointPartialsProvider
    }

    @Override
    public PrecisionType getPrecisionType() {
        return precisionType;
    }

    @Override
    public double[] getTipPartial(int taxonIndex, boolean fullyObserved) {
        double[] partial = new double[precisionType.getPartialsDimension(traitDim)];

        int meanOffset = precisionType.getMeanOffset(traitDim);
        int precOffset = precisionType.getPrecisionOffset(traitDim);
        int varOffset = precisionType.getVarianceOffset(traitDim);
        int effDimDim = precisionType.getEffectiveDimensionOffset(traitDim);

        for (ContinuousTraitPartialsProvider provider : providers) {
            double[] subPartial = getTipPartial(taxonIndex, fullyObserved);
            int subDim = provider.getTraitDimension();
            int subPrecDim = precisionType.getPrecisionLength(subDim);
            int subVarDim = precisionType.getVarianceLength(subDim);
            int matrixIncrement = getMatrixIncrement(traitDim, subDim);

            System.arraycopy(subPartial, precisionType.getMeanOffset(subDim), partial, meanOffset, subDim);
            meanOffset += subDim;

            System.arraycopy(subPartial, precisionType.getPrecisionOffset(subDim), partial, precOffset, subPrecDim);
            precOffset += matrixIncrement;

            System.arraycopy(subPartial, precisionType.getVarianceOffset(subDim), partial, varOffset, subVarDim);
            varOffset += matrixIncrement;

            if (precisionType.hasEffectiveDimension()) {
                partial[effDimDim] += subPartial[precisionType.getEffectiveDimensionOffset(subDim)];
            }
        }
        //Assume conditional independence (for now)
        return partial;
    }

    private int getMatrixIncrement(int fullDim, int subDim) { //TODO
        return 0;
    }


    @Override
    public List<Integer> getMissingIndices() {
        return null; //TODO: how to merge missing indices
    }

    @Override
    public boolean[] getMissingIndicator() {
        return new boolean[0]; //TODO: see above
    }

    @Override
    public CompoundParameter getParameter() {
        throw new RuntimeException("not implemented"); //TODO: This is going to be the real problem, I think
    }

    @Override
    public String getModelName() {
        return name;
    }


}
