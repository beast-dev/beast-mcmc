package dr.inferencexml.distribution;

import dr.inference.distribution.MomentDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Created by max on 5/14/15.
 */
public class MomentDistributionModelParser extends AbstractXMLObjectParser {
    public static final String MOMENT_DISTRIBUTION_MODEL = "momentDistributionModel";
    public static final String MEAN = "mean";
//    public static final String STDEV = "stdev";
    public static final String PREC = "precision";
    public static final String CUTOFF="cutoff";
    public static final String DATA="data";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter mean=(Parameter) xo.getChild(MEAN).getChild(0);
        Parameter prec=(Parameter) xo.getChild(PREC).getChild(0);
        Parameter cutoff=(Parameter) xo.getChild(CUTOFF).getChild(0);
        Parameter data=(Parameter) xo.getChild(DATA).getChild(0);

        return new MomentDistributionModel(mean, prec, cutoff, data);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MEAN,
                    new XMLSyntaxRule[]{
                                            new ElementRule(Parameter.class)
                            }
            ),
                    new ElementRule(CUTOFF,
                            new XMLSyntaxRule[]{
                                            new ElementRule(Parameter.class)
                                    }
                    ),
                    new ElementRule(PREC,
                            new XMLSyntaxRule[]{
                                            new ElementRule(Parameter.class)
                                    }
            ),
            new ElementRule(DATA,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }
            )
    };

    @Override
    public String getParserDescription() {
        return "Returns an internally truncated normal distribution for claculating a moment prior";
    }

    @Override
    public Class getReturnType() {
        return MomentDistributionModel.class;
    }

    @Override
    public String getParserName() {
        return MOMENT_DISTRIBUTION_MODEL;
    }
}
