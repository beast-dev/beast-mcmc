package dr.evomodelxml.continuous;

import dr.evomodel.continuous.DenseBandedMultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.inference.model.MatrixParameterInterface;
import dr.math.distributions.GaussianMarkovRandomField;
import dr.xml.*;

public class DenseBandedMultivariateDiffusionModelParser extends AbstractXMLObjectParser {

    private static final String DIFFUSION_PROCESS = "bandedMultivariateDiffusionModel";
    private static final String DIFFUSION_CONSTANT = "precisionMatrix";
    private final static String REPLICATES = "replicates";

    public String getParserName() { return DIFFUSION_PROCESS; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(DIFFUSION_CONSTANT);
        MatrixParameterInterface diffusionParam = (MatrixParameterInterface)
                cxo.getChild(MatrixParameterInterface.class);

        MultivariateDiffusionModel.checkIsPositiveDefinite(diffusionParam.getParameterAsMatrix());

        if (xo.hasAttribute(REPLICATES)) {
            int replicates = xo.getIntegerAttribute(REPLICATES);
            return new DenseBandedMultivariateDiffusionModel(diffusionParam, replicates);
        } else {
            GaussianMarkovRandomField field = (GaussianMarkovRandomField)
                    xo.getChild(GaussianMarkovRandomField.class);
            return new DenseBandedMultivariateDiffusionModel(diffusionParam, field);
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Describes a multivariate normal diffusion process.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(DIFFUSION_CONSTANT,
                    new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)}),
            new XORRule(
                    AttributeRule.newIntegerRule(REPLICATES),
                    new ElementRule(GaussianMarkovRandomField.class)),
    };

    public Class getReturnType() {
        return DenseBandedMultivariateDiffusionModel.class;
    }
}
