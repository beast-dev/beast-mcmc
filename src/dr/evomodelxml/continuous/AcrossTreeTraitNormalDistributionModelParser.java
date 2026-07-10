package dr.evomodelxml.continuous;

import dr.evomodel.continuous.AcrossTreeTraitNormalDistributionModel;
import dr.evomodel.continuous.TreeTraitNormalDistributionModel;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class AcrossTreeTraitNormalDistributionModelParser extends AbstractXMLObjectParser {

    private static final String ACROSS_TREE_TRAIT_NORMAL = "acrossTreeTraitNormalDistribution";

    public String getParserName() {
        return ACROSS_TREE_TRAIT_NORMAL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<TreeDataLikelihood> likelihoods = xo.getAllChildren(TreeDataLikelihood.class);
        List<ContinuousDataLikelihoodDelegate> delegates = new ArrayList<>();
        for (TreeDataLikelihood likelihood : likelihoods) {
            DataLikelihoodDelegate delegate = likelihood.getDataLikelihoodDelegate();
            if (delegate instanceof ContinuousDataLikelihoodDelegate) {
                delegates.add((ContinuousDataLikelihoodDelegate) delegate);
            }
        }

        Parameter rho = (Parameter) xo.getChild(Parameter.class);

        return new AcrossTreeTraitNormalDistributionModel(delegates.get(0), delegates.get(1), rho);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeDataLikelihood.class, 2, 2),
            new ElementRule(Parameter.class),
    };

    public String getParserDescription() {
        return "Parses TreeTraitNormalDistributionModel";
    }

    public Class getReturnType() {
        return TreeTraitNormalDistributionModel.class;
    }
}
