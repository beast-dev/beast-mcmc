package dr.evomodelxml.continuous;

import dr.evomodel.continuous.AbstractBandedMultivariateDiffusionModel;
import dr.evomodel.continuous.DenseBandedMultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.SparseBandedMultivariateDiffusionModel;
import dr.inference.model.MatrixParameterInterface;
import dr.math.distributions.GaussianMarkovRandomField;
import dr.xml.*;

public class BandedMultivariateDiffusionModelParser extends AbstractXMLObjectParser {

    private static final String DIFFUSION_PROCESS = "bandedMultivariateDiffusionModel";
    private static final String DIFFUSION_CONSTANT = "precisionMatrix";
    private static final String REPLICATES = "replicates";
    private static final String REPRESENTATION = "representation";

    public String getParserName() { return DIFFUSION_PROCESS; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(DIFFUSION_CONSTANT);
        MatrixParameterInterface diffusionParam = (MatrixParameterInterface)
                cxo.getChild(MatrixParameterInterface.class);

        MultivariateDiffusionModel.checkIsPositiveDefinite(diffusionParam.getParameterAsMatrix());

        boolean isSparse = false;
        String representation = xo.getAttribute(REPRESENTATION, "dense");
        if (representation.equalsIgnoreCase("sparse")) {
            isSparse = true;
        } else if (!representation.equalsIgnoreCase("dense")) {
            throw new XMLParseException("Unknown internal representation");
        }

        if (xo.hasAttribute(REPLICATES)) {
            int replicates = xo.getIntegerAttribute(REPLICATES);
            if (isSparse) {
                return new SparseBandedMultivariateDiffusionModel(diffusionParam, replicates);
            } else {
                return new DenseBandedMultivariateDiffusionModel(diffusionParam, replicates);
            }
        } else {
            GaussianMarkovRandomField field = (GaussianMarkovRandomField)
                    xo.getChild(GaussianMarkovRandomField.class);
            if (isSparse) {
                return new SparseBandedMultivariateDiffusionModel(diffusionParam, field);
            } else {
                return new DenseBandedMultivariateDiffusionModel(diffusionParam, field);
            }
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
            AttributeRule.newStringRule(REPRESENTATION, true),
            new ElementRule(DIFFUSION_CONSTANT,
                    new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)}),
            new XORRule(
                    AttributeRule.newIntegerRule(REPLICATES),
                    new ElementRule(GaussianMarkovRandomField.class)),
    };

    public Class getReturnType() { return AbstractBandedMultivariateDiffusionModel.class; }
}
