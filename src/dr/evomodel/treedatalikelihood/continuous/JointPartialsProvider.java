package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.CompoundParameter;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.xml.*;

import java.util.List;

import static dr.math.matrixAlgebra.WrappedMatrix.Utils.transferSymmetricBlockDiagonal;
import static dr.math.matrixAlgebra.WrappedMatrix.Utils.wrapBlockDiagonalMatrix;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class JointPartialsProvider implements ContinuousTraitPartialsProvider {

    private final String name;
    private final ContinuousTraitPartialsProvider[] providers;
    private final int traitDim;
    private final int dataDim;

    private final List<Integer> missingIndices;
    private final boolean[] missingIndicators;

    private static final PrecisionType precisionType = PrecisionType.FULL; //TODO: base on child precisionTypes (make sure they're all the same)

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

        this.missingIndicators = setupMissingIndicators();
        this.missingIndices = ContinuousTraitPartialsProvider.indicatorToIndices(missingIndicators);

    }


    private boolean[] setupMissingIndicators() {
        int nTaxa = providers[0].getParameter().getParameterCount();
        boolean[] indicators = new boolean[dataDim * nTaxa];
        boolean[][] subIndicators = new boolean[providers.length][0];
        for (int i = 0; i < providers.length; i++) {
            subIndicators[i] = providers[i].getMissingIndicator();
        }
        for (int taxonI = 0; taxonI < nTaxa; taxonI++) {
            int offset = taxonI * dataDim;

            for (int providerI = 0; providerI < providers.length; providerI++) {
                int srcDim = providers[providerI].getDataDimension();
                int srcOffset = taxonI * srcDim;
                System.arraycopy(subIndicators[providerI], srcOffset, indicators, offset, srcDim);
                offset += srcDim;
            }
        }

        return indicators;
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
        return traitDim;
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

        WrappedMatrix.Indexed precWrap = wrapBlockDiagonalMatrix(partial, precOffset, 0, traitDim); //TODO: this only works for precisionType.FULL, make general
        WrappedMatrix.Indexed varWrap = wrapBlockDiagonalMatrix(partial, varOffset, 0, traitDim); //TODO: see above

        int currentMatrixOffset = 0;


        for (ContinuousTraitPartialsProvider provider : providers) {
            double[] subPartial = provider.getTipPartial(taxonIndex, fullyObserved);
            int subDim = provider.getTraitDimension();

            WrappedMatrix.Raw subPrec = new WrappedMatrix.Raw(subPartial, precisionType.getPrecisionOffset(subDim), subDim, subDim); //TODO: see above
            transferSymmetricBlockDiagonal(subPrec, precWrap, currentMatrixOffset); //TODO: see above

            WrappedMatrix.Raw subVar = new WrappedMatrix.Raw(subPartial, precisionType.getVarianceOffset(subDim), subDim, subDim); //TODO: see above
            transferSymmetricBlockDiagonal(subVar, varWrap, currentMatrixOffset); //TODO: see above

            currentMatrixOffset += subDim;


            System.arraycopy(subPartial, precisionType.getMeanOffset(subDim), partial, meanOffset, subDim);
            meanOffset += subDim;

            if (precisionType.hasEffectiveDimension()) {
                partial[effDimDim] += subPartial[precisionType.getEffectiveDimensionOffset(subDim)];
            }
        }
        //Assume conditional independence (for now)
        return partial;
    }


    @Override
    public List<Integer> getMissingIndices() {
        return missingIndices; //TODO: how to merge missing indices
    }

    @Override
    public boolean[] getMissingIndicator() {
        return missingIndicators; //TODO: see above
    }

    @Override
    public CompoundParameter getParameter() {
        System.err.println("Warning: This is broken. (JointPartialsProvider.getParameter())");
        return providers[0].getParameter(); //TODO: This is going to be the real problem, I think
    }

    @Override
    public String getModelName() {
        return name;
    }


    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        private static final String PARSER_NAME = "jointPartialsProvider";


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            List<ContinuousTraitPartialsProvider> providersList =
                    xo.getAllChildren(ContinuousTraitPartialsProvider.class);

            ContinuousTraitPartialsProvider[] providers = new ContinuousTraitPartialsProvider[providersList.size()];

            for (int i = 0; i < providersList.size(); i++) {
                providers[i] = providersList.get(i);
            }
            return new JointPartialsProvider(PARSER_NAME, providers);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(ContinuousTraitPartialsProvider.class, 0, Integer.MAX_VALUE)
            };
        }

        @Override
        public String getParserDescription() {
            return "Merges two Gaussian processes.";
        }

        @Override
        public Class getReturnType() {
            return JointPartialsProvider.class;
        }

        @Override
        public String getParserName() {
            return PARSER_NAME;
        }
    };


}
