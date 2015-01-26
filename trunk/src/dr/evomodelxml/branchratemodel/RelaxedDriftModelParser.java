package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.RelaxedDriftModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Created by IntelliJ IDEA.
 * User: mandevgill
 * Date: 7/28/14
 * Time: 3:24 PM
 * To change this template use File | Settings | File Templates.
 */
public class RelaxedDriftModelParser extends AbstractXMLObjectParser {

    public static final String RELAXED_DRIFT = "relaxedDriftModel";
    public static final String RATES = "rates";
    public static final String RATE_IND = "rateIndicator";
    public static final String DRIFT_RATES = "driftRates";

    public String getParserName() {
        return RELAXED_DRIFT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        Parameter ratesParameter = (Parameter) xo.getElementFirstChild(RATES);

        Parameter rateIndicatorParameter = (Parameter) xo.getElementFirstChild(RATE_IND);

        Parameter driftRates = null;
        if (xo.hasChildNamed(DRIFT_RATES)) {
            driftRates = (Parameter) xo.getElementFirstChild(DRIFT_RATES);
        }


        Logger.getLogger("dr.evomodel").info("Using relaxed drift model.");


        return new RelaxedDriftModel(tree, rateIndicatorParameter,
                ratesParameter, driftRates);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element returns a relaxed drift model.";
    }

    public Class getReturnType() {
        return RelaxedDriftModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TreeModel.class),
            new ElementRule(RATES, Parameter.class, "The rates parameter", false),
            new ElementRule(RATE_IND, Parameter.class, "The indicator parameter", false),
            new ElementRule(DRIFT_RATES, Parameter.class, "the drift rates parameter", true)
    };
}
