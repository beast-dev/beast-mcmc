package dr.inferencexml.model;

import dr.inference.model.AdaptableSizeFastMatrixParameter;
import dr.inference.model.DeterminentalPointProcessPrior;
import dr.xml.*;

/**
 * Created by max on 4/6/16.
 */
public class DeterminentalPointProcessPriorParser extends AbstractXMLObjectParser {
    public static final String DETERMINENTAL_POINT_PROCESS_PRIOR="determinentalPointProcessPrior";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        return new DeterminentalPointProcessPrior();
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
//            new ElementRule(Parameter.class, 0, Integer.MAX_VALUE),
        new ElementRule(AdaptableSizeFastMatrixParameter.class),

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

