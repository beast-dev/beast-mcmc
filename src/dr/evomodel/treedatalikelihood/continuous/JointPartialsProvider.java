package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.CompoundParameter;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.math.matrixAlgebra.missingData.MissingOps;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;

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

    private final boolean defaultAllowSingular;
    private final Boolean computeDeterminant; // TODO: Maybe pass as argument?

    private static final PrecisionType precisionType = PrecisionType.FULL; //TODO: base on child precisionTypes (make sure they're all the same)

    private String tipTraitName;

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

        this.defaultAllowSingular = setDefaultAllowSingular();
        this.computeDeterminant = defaultAllowSingular; // TODO: not perfect behavior, should be based on actual value of `allowSingular`
    }


    private boolean[] setupMissingIndicators() {
        int nTaxa = providers[0].getParameter().getParameterCount();
        boolean[] indicators = new boolean[dataDim * nTaxa];
        boolean[][] subIndicators = new boolean[providers.length][0];
        for (int i = 0; i < providers.length; i++) {
            subIndicators[i] = providers[i].getDataMissingIndicators();
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
    public String getTipTraitName() {
        return tipTraitName;
    }

    @Override
    public void setTipTraitName(String name) {
        tipTraitName = name;
        for (int i = 0; i < providers.length; i++) {
            providers[i].setTipTraitName(name + "." + i); // TODO: make static method for making name both here and in PartitionedTreeTraitProvider;
        }
    }

    @Override
    public int getDataDimension() {
        return dataDim; //TODO: maybe throw error here? Used for model extension, mse stuff and it might be worth putting conditions if JointPartialsProvider
    }

    @Override
    public int[] getPartitionDimensions() {
        int[] dims = new int[providers.length];
        for (int i = 0; i < providers.length; i++) {
            dims[i] = providers[i].getTraitDimension();
        }
        return dims;
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
        int detDim = precisionType.getDeterminantOffset(traitDim);

        WrappedMatrix.Indexed precWrap = wrapBlockDiagonalMatrix(partial, precOffset, 0, traitDim); //TODO: this only works for precisionType.FULL, make general
        WrappedMatrix.Indexed varWrap = wrapBlockDiagonalMatrix(partial, varOffset, 0, traitDim); //TODO: see above

        int currentMatrixOffset = 0;


        for (ContinuousTraitPartialsProvider provider : providers) {
            double[] subPartial = provider.getTipPartial(taxonIndex, fullyObserved);
            int subDim = provider.getTraitDimension();

            int precisionOffset = precisionType.getPrecisionOffset(subDim);

            WrappedMatrix.Raw subPrec = new WrappedMatrix.Raw(subPartial, precisionOffset, subDim, subDim); //TODO: see above
            transferSymmetricBlockDiagonal(subPrec, precWrap, currentMatrixOffset); //TODO: see above

            WrappedMatrix.Raw subVar = new WrappedMatrix.Raw(subPartial, precisionType.getVarianceOffset(subDim), subDim, subDim); //TODO: see above
            transferSymmetricBlockDiagonal(subVar, varWrap, currentMatrixOffset); //TODO: see above

            currentMatrixOffset += subDim;


            System.arraycopy(subPartial, precisionType.getMeanOffset(subDim), partial, meanOffset, subDim);
            meanOffset += subDim;

            if (precisionType.hasEffectiveDimension()) {
                partial[effDimDim] += subPartial[precisionType.getEffectiveDimensionOffset(subDim)];
            }

            if (precisionType.hasEffectiveDimension() && computeDeterminant) {

                double subDet = subPartial[precisionType.getDeterminantOffset(subDim)];

                if (!precisionType.isMissingDeterminantValue(subDet)) {

                    DenseMatrix64F prec = MissingOps.wrap(subPartial, precisionOffset, subDim, subDim);
                    DenseMatrix64F var = new DenseMatrix64F(subDim, subDim);
                    subDet = MissingOps.safeInvert2(prec, var, true).getLogDeterminant();
                }

                partial[detDim] += subDet;
            }
        }

        if (!computeDeterminant) {
            precisionType.fillNoDeterminantInPartials(partial, 0, traitDim);
        }

        //Assume conditional independence (for now)
        return partial;
    }


    @Override
    public List<Integer> getMissingIndices() {
        return missingIndices; //TODO: how to merge missing indices
    }

    @Override
    public boolean[] getDataMissingIndicators() {
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

    @Override
    public boolean getDefaultAllowSingular() {
        return defaultAllowSingular;
    }

    private boolean setDefaultAllowSingular() {
        boolean allowSingular = false;
        for (ContinuousTraitPartialsProvider provider : providers) {
            allowSingular = allowSingular || provider.getDefaultAllowSingular();
        }
        return allowSingular;
    }

    @Override
    public boolean suppliesWishartStatistics() {
        boolean suppliesStatistics = true;
        for (ContinuousTraitPartialsProvider provider : providers) {
            suppliesStatistics = suppliesStatistics && provider.suppliesWishartStatistics();
        }
        return suppliesStatistics;
    }

    @Override
    public void addTreeAndRateModel(Tree treeModel, ContinuousRateTransformation rateTransformation) {
        for (ContinuousTraitPartialsProvider provider : providers) {
            provider.addTreeAndRateModel(treeModel, rateTransformation);
        }
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
