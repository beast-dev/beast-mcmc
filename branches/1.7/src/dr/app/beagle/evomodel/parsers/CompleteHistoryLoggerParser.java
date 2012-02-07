package dr.app.beagle.evomodel.parsers;

import dr.app.beagle.evomodel.treelikelihood.MarkovJumpsBeagleTreeLikelihood;
import dr.app.beagle.evomodel.utilities.CompleteHistoryLogger;
import dr.inference.loggers.Logger;
import dr.xml.*;

/**
 * @author Marc A. Suchrd
 * @author Philippe Lemey
 */
public class CompleteHistoryLoggerParser extends AbstractXMLObjectParser {

    public static final String NAME = "completeHistoryLogger";
 
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MarkovJumpsBeagleTreeLikelihood treeLikelihood =
                (MarkovJumpsBeagleTreeLikelihood) xo.getChild(MarkovJumpsBeagleTreeLikelihood.class);
        return new CompleteHistoryLogger(treeLikelihood);
    }

    public String getParserName() {
        return NAME;
    }

    public String getParserDescription() {
        return "A logger to record all transitions in the complete history.";
    }

    public Class getReturnType() {
        return Logger.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MarkovJumpsBeagleTreeLikelihood.class),
    };
}
