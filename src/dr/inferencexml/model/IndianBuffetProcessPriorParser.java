package dr.inferencexml.model;

import dr.inference.model.IndianBuffetProcessPrior;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Max Tolkoff
 */
public class IndianBuffetProcessPriorParser extends AbstractXMLObjectParser {
    public static final String INDIAN_BUFFET_PROCESS="IndianBuffetProcess";
    public static final String BETA="beta";
    public static final String ALPHA="alpha";
    public static final String DATA="data";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter alpha=(Parameter) xo.getChild(ALPHA).getChild(0);
        MatrixParameter data=(MatrixParameter) xo.getChild(DATA).getChild(0);
        Parameter beta;
        if(xo.hasChildNamed(BETA))
        {
            beta=(Parameter) xo.getChild(BETA).getChild(0);
        }
        else
        {
            beta=new Parameter.Default(1);
        }
        return new IndianBuffetProcessPrior(alpha, beta, data);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(BETA,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }
            , true),
            new ElementRule(ALPHA,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }
            ),
            new ElementRule(DATA,
                    new XMLSyntaxRule[]{
                            new ElementRule(MatrixParameter.class)
                    }
            )
    };

    @Override
    public String getParserDescription() {
        return "Indian Buffet Process prior on a binary matrix parameter";
    }

    @Override
    public Class getReturnType() {
        return IndianBuffetProcessPrior.class;
    }

    @Override
    public String getParserName() {
        return INDIAN_BUFFET_PROCESS;
    }
}
