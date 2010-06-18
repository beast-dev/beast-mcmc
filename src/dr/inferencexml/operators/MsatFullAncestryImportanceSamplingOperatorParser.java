package dr.inferencexml.operators;

import dr.xml.*;
import dr.inference.operators.MsatFullAncestryImportanceSamplingOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.model.Parameter;
import dr.evomodel.tree.MicrosatelliteSamplerTreeModel;
import dr.evomodel.substmodel.MicrosatelliteModel;
import dr.evomodel.branchratemodel.BranchRateModel;

/**
 * @author Chieh-Hsi Wu
 *
 * Parser for MsatFullAncestryGibbsOperator
 */
public class MsatFullAncestryImportanceSamplingOperatorParser extends AbstractXMLObjectParser {
    public String getParserName(){
        return MsatFullAncestryImportanceSamplingOperator.MSAT_FULL_ANCESTRY_IMPORTANCE_SAMPLING_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final Parameter parameter = (Parameter)xo.getChild(Parameter.class);
        final MicrosatelliteSamplerTreeModel msatSamplerTreeModel = (MicrosatelliteSamplerTreeModel)xo.getChild(MicrosatelliteSamplerTreeModel.class);
        final MicrosatelliteModel msatModel = (MicrosatelliteModel)xo.getChild(MicrosatelliteModel.class);
        final BranchRateModel branchRateModel = (BranchRateModel)xo.getChild(BranchRateModel.class);

        return new MsatFullAncestryImportanceSamplingOperator(parameter, msatSamplerTreeModel, msatModel, branchRateModel,weight);
    }

    public String getParserDescription() {
        return "This element represents an operator that samples the full ancestry given a microsatellite pattern and a tree";
    }

    public Class getReturnType(){
        return MsatFullAncestryImportanceSamplingOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(Parameter.class),
            new ElementRule(MicrosatelliteSamplerTreeModel.class),
            new ElementRule(MicrosatelliteModel.class),
            new ElementRule(BranchRateModel.class)
    };
}
