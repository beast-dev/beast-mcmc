package dr.evomodelxml.continuous;

import dr.evomodel.continuous.LoadingsShrinkagePrior;
import dr.inference.distribution.shrinkage.BayesianBridgeDistributionModel;
import dr.inference.distribution.shrinkage.BayesianBridgeLikelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inferencexml.distribution.shrinkage.BayesianBridgeDistributionModelParser;
import dr.inferencexml.distribution.shrinkage.BayesianBridgeLikelihoodParser;
import dr.xml.*;


/**
 * @author Gabriel Hassler
 */


public class LoadingsShrinkagePriorParser extends AbstractXMLObjectParser {
    private static final String LOADINGS_SHRINKAGE = "loadingsShrinkagePrior";
    private static final String ROW_PRIORS = "rowPriors";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute("id", LoadingsShrinkagePrior.class.toString());

        MatrixParameterInterface loadings = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
        XMLObject rpxo = xo.getChild(ROW_PRIORS);

        BayesianBridgeDistributionModel[] rowModels = new BayesianBridgeDistributionModel[loadings.getRowDimension()];

        int counter = 0;


        for (int i = 0; i < rpxo.getChildCount(); i++) {
            Object cxo = rpxo.getChild(i);
            if (cxo instanceof BayesianBridgeDistributionModel) {
                rowModels[i] = (BayesianBridgeDistributionModel) cxo;
                counter++;
            } else {
                throw new XMLParseException("The " + ROW_PRIORS + " element should only have "
                        + BayesianBridgeDistributionModelParser.BAYESIAN_BRIDGE_DISTRIBUTION + " child elements.");
            }
        }

        if (counter != loadings.getColumnDimension()) {
            throw new XMLParseException("The number of " + BayesianBridgeDistributionModelParser.BAYESIAN_BRIDGE_DISTRIBUTION +
                    " in " + ROW_PRIORS + " must be the same as the number of rows in the matrix " +
                    loadings.getParameterName() + ".\nThere are currently " + counter + " " +
                    BayesianBridgeDistributionModelParser.BAYESIAN_BRIDGE_DISTRIBUTION + " elements and " + loadings.getRowDimension() +
                    " rows in " + loadings.getParameterName() + ".");
        }


        return new LoadingsShrinkagePrior(name, loadings, rowModels);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(MatrixParameterInterface.class),
                new ElementRule(ROW_PRIORS, new XMLSyntaxRule[]{
                        new ElementRule(BayesianBridgeDistributionModel.class, 1, Integer.MAX_VALUE)
                })
        };
    }

    @Override
    public String getParserDescription() {
        return "Bayesian bridge shrinkage prior on the loadings matrix of a latent factor model";
    }

    @Override
    public Class getReturnType() {
        return LoadingsShrinkagePrior.class;
    }

    @Override
    public String getParserName() {
        return LOADINGS_SHRINKAGE;
    }
}
