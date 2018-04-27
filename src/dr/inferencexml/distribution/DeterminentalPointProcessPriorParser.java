package dr.inferencexml.distribution;

import dr.inference.model.AdaptableSizeFastMatrixParameter;
import dr.inference.distribution.DeterminentalPointProcessPrior;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Created by max on 4/6/16.
 */
public class DeterminentalPointProcessPriorParser extends AbstractXMLObjectParser {
    public static final String DETERMINENTAL_POINT_PROCESS_PRIOR="determinentalPointProcessPrior";
    public static final String THETA = "theta";
    public static final String NORMALIZING_CONSTANTS = "normalizingConstants";
    public static final String PATH_SAMPLING = "pathSampling";
    public static final String NO_ZEROS = "noZeros";
    public static final String RESET_DATA = "resetData";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getName();
        double theta = xo.getDoubleAttribute(THETA);
        MatrixParameterInterface data = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
        Parameter normalizingConstants = null;
        if(xo.getChild(NORMALIZING_CONSTANTS) != null) {
             normalizingConstants = (Parameter) xo.getChild(NORMALIZING_CONSTANTS).getChild(Parameter.class);
        }
        boolean noZeros = xo.getAttribute(NO_ZEROS, false);
        boolean pathSampling = xo.getAttribute(PATH_SAMPLING, false);
        boolean resetData = xo.getAttribute(RESET_DATA, true);

        return new DeterminentalPointProcessPrior(name, theta, data, normalizingConstants, noZeros, pathSampling, resetData);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
//            new ElementRule(Parameter.class, 0, Integer.MAX_VALUE),
        new ElementRule(MatrixParameterInterface.class),
            AttributeRule.newDoubleRule(THETA),
            AttributeRule.newBooleanRule(NO_ZEROS, true),
            AttributeRule.newBooleanRule(PATH_SAMPLING, true),
            AttributeRule.newBooleanRule(RESET_DATA, true),
            new ElementRule(NORMALIZING_CONSTANTS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
    };



    @Override
    public String getParserDescription() {
        return "Returns a blockUpperTriangularMatrixParameter which is a compoundParameter which forces the last element to be of full length, the second to last element to be of full length-1, etc.";  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Class getReturnType() {
        return DeterminentalPointProcessPrior.class;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getParserName() {
        return DETERMINENTAL_POINT_PROCESS_PRIOR;  //To change body of implemented methods use File | Settings | File Templates.
    }
}

