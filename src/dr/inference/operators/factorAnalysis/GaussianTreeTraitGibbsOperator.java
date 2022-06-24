package dr.inference.operators.factorAnalysis;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class GaussianTreeTraitGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String TRAIT_GIBBS = "gaussianTreeTraitGibbsOperator";
    private final TreeTrait<double[]> treeTrait;
    private final Parameter traitParameter;
    private final TreeDataLikelihood treeDataLikelihood;

    public GaussianTreeTraitGibbsOperator(TreeDataLikelihood treeDataLikelihood, Parameter parameter, String traitName) {
        this.traitParameter = parameter;
        this.treeDataLikelihood = treeDataLikelihood;
        ContinuousDataLikelihoodDelegate delegate = (ContinuousDataLikelihoodDelegate) treeDataLikelihood.getDataLikelihoodDelegate();
        this.treeTrait = treeDataLikelihood.getTreeTrait(delegate.getDataModel().getTipTraitName());
    }


    @Override
    public String getOperatorName() {
        return TRAIT_GIBBS;
    }

    @Override
    public double doOperation() {
        treeDataLikelihood.fireModelChanged();
        double[] traits = treeTrait.getTrait(treeDataLikelihood.getTree(), null);
        traitParameter.setAllParameterValuesQuietly(traits);
        traitParameter.fireParameterChangedEvent();
        return Double.POSITIVE_INFINITY;
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            String traitName = xo.getStringAttribute(TreeTraitParserUtilities.TRAIT_NAME);
            TreeDataLikelihood likelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            return new GaussianTreeTraitGibbsOperator(likelihood, parameter, traitName);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[0];
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return GaussianTreeTraitGibbsOperator.class;
        }

        @Override
        public String getParserName() {
            return TRAIT_GIBBS;
        }
    };

}
