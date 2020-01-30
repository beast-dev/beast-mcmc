package dr.evomodelxml.continuous;

import dr.evomodel.continuous.LoadingsShrinkagePrior;
import dr.inference.distribution.shrinkage.BayesianBridgeLikelihood;
import dr.inference.model.MatrixParameterInterface;
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

        String name = (String) xo.getAttribute("id", LoadingsShrinkagePrior.class.toString());

        MatrixParameterInterface loadings = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
        XMLObject rpxo = xo.getChild(ROW_PRIORS);

        BayesianBridgeLikelihood[] rowPriors = new BayesianBridgeLikelihood[loadings.getRowDimension()];

        int counter = 0;


        for (int i = 0; i < loadings.getRowDimension(); i++) {
            Object cxo = rpxo.getChild(i);
            if (cxo instanceof BayesianBridgeLikelihood) {
                rowPriors[i] = (BayesianBridgeLikelihood) cxo;
                counter++;
            } else {
                throw new XMLParseException("The " + ROW_PRIORS + " element should only have "
                        + BayesianBridgeLikelihoodParser.BAYESIAN_BRIDGE + " child elements.");
            }
        }

        if (counter != loadings.getRowDimension()) {
            throw new XMLParseException("The number of " + BayesianBridgeLikelihoodParser.BAYESIAN_BRIDGE +
                    " in " + ROW_PRIORS + " must be the same as the number of rows in the matrix " +
                    loadings.getParameterName() + ". There are currently " + counter +
                    BayesianBridgeLikelihoodParser.BAYESIAN_BRIDGE + " elements and " + loadings.getRowDimension() +
                    " rows in " + loadings.getParameterName() + ".");
        }


        return new LoadingsShrinkagePrior(name, loadings, rowPriors);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(MatrixParameterInterface.class),
                new ElementRule(ROW_PRIORS, new XMLSyntaxRule[]{
                        new ElementRule(BayesianBridgeLikelihood.class, 1, Integer.MAX_VALUE)
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
