package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.ScaledTreeLengthRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 */
public class ScaledTreeLengthRateModelParser extends AbstractXMLObjectParser {
    public static final String MODEL_NAME = "scaledTreeLengthModel";
    public static final String SCALING_FACTOR = "scalingFactor";

    public String getParserName() {
        return MODEL_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        Parameter totalLength = (Parameter) xo.getElementFirstChild(SCALING_FACTOR);
        if (totalLength == null) {
            totalLength = new Parameter.Default(1, 1.0);
        }
        Logger.getLogger("dr.evomodel.branchratemodel").info("\n ---------------------------------\nCreating ScaledTreeLengthRateModel model.");
        Logger.getLogger("dr.evomodel.branchratemodel").info("\tTotal tree length will be scaled to " + totalLength.getParameterValue(0) + ".");
        Logger.getLogger("dr.evomodel.branchratemodel").info("\tIf you publish results using this rate model, please reference Alekseyenko, Lee and Suchard (in submision).\n---------------------------------\n");

        return new ScaledTreeLengthRateModel(tree, totalLength);
    }

    public String getParserDescription() {
        return "This element returns a branch rate model that scales the total length " +
                "of the tree to specified valued (default=1.0).";
    }

    public Class getReturnType() {
        return ScaledTreeLengthRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(SCALING_FACTOR,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(TreeModel.class)
    };

}
