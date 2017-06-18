package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.BNPRSamplingLikelihood;
import dr.evomodel.coalescent.DemographicModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Created by mkarcher on 3/10/17.
 */
public class BNPRSamplingLikelihoodParser extends AbstractXMLObjectParser {
    public static final String SAMPLING_LIKELIHOOD = "bnprSamplingLikelihood";
    public static final String MODEL = "model";
    public static final String BETAS = "betas";
    public static final String POPULATION_TREE = "populationTree";
    public static final String EPOCH_WIDTHS = "epochWidths";
    public static final String WIDTHS = "widths";

    @Override
    public String getParserName() {
        return SAMPLING_LIKELIHOOD;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(MODEL);
        DemographicModel demoModel = (DemographicModel) cxo.getChild(DemographicModel.class);

        cxo = xo.getChild(BETAS);
        Parameter betas = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(POPULATION_TREE); // May need to adapt to multiple trees, a la CoalescentLikelihoodParser
        TreeModel tree = (TreeModel) cxo.getChild(TreeModel.class);

        double[] epochWidths = null;
        if (xo.hasChildNamed(EPOCH_WIDTHS)) {
            cxo = xo.getChild(EPOCH_WIDTHS);
            epochWidths = cxo.getDoubleArrayAttribute(WIDTHS);
        }

        return new BNPRSamplingLikelihood(tree, betas, demoModel, epochWidths);
    }

    @Override
    public String getParserDescription() {
        return "This element represents the likelihood of the sampling times given the demographic function.";
    }

    @Override
    public Class getReturnType() {
        return BNPRSamplingLikelihood.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MODEL, new XMLSyntaxRule[]{
                    new ElementRule(DemographicModel.class)
            }, "The demographic model which describes the effective population size over time"),

            new ElementRule(BETAS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, "The log-linear coefficients of effective population size, used to calculate sampling intensity"),

            new ElementRule(POPULATION_TREE, new XMLSyntaxRule[]{
                    new ElementRule(TreeModel.class)
            }, "Tree/sampling times to compute likelihood for"),

            new ElementRule(EPOCH_WIDTHS,
                    new XMLSyntaxRule[]{AttributeRule.newDoubleArrayRule(WIDTHS)}, true)
    };
}
