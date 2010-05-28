package dr.app.beagle.evomodel.parsers;

import dr.xml.*;
import dr.app.beagle.evomodel.utilities.StratifiedTraitLogger;
import dr.app.beagle.evomodel.substmodel.CodonPartitionedRobustCounting;
import dr.app.beagle.evomodel.substmodel.StratifiedTraitOutputFormat;
import dr.evomodel.tree.TreeModel;

/**
 * @author Marc A. Suchard
 */

public class StratifiedTraitLoggerParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "stratifiedTraitLogger";
    public static final String LOG_FORMAT = "format";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        CodonPartitionedRobustCounting trait =
                (CodonPartitionedRobustCounting) xo.getChild(CodonPartitionedRobustCounting.class);

        String logFormatString = xo.getAttribute(LOG_FORMAT,
                StratifiedTraitOutputFormat.SUM_OVER_SITES.getText());
        StratifiedTraitOutputFormat logFormat = StratifiedTraitOutputFormat.parseFromString(logFormatString);
        if (logFormat == null) {
            throw new XMLParseException("Unrecognized log output format '" + logFormatString + "'");
        }

        return new StratifiedTraitLogger(xo.getId(), treeModel, trait, logFormat);
    }

  private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
          AttributeRule.newStringRule(LOG_FORMAT, true),
          new ElementRule(TreeModel.class),
          new ElementRule(CodonPartitionedRobustCounting.class), // TODO Will be trait implementation
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    public String getParserDescription() {
        return "A parser to stratify traits to be logged";
    }

    public Class getReturnType() {
        return StratifiedTraitLogger.class;
    }

    public String getParserName() {
        return PARSER_NAME;
    }
}
