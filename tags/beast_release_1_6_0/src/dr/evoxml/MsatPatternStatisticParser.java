package dr.evoxml;

import dr.xml.*;
import dr.evolution.alignment.Patterns;
import dr.evolution.alignment.MsatPatternStatistic;

/**
 * @author Chieh-Hsi Wu
 * Parser for computing the statistics of msat pattern
 */
public class MsatPatternStatisticParser extends AbstractXMLObjectParser {
    public static final String MSAT_PATTERN_STATISTIC_PARSER = "msatPatternStatistic";
    public static final String MODE = "mode";

    public String getParserName(){
        return MSAT_PATTERN_STATISTIC_PARSER;
    }


    public Object parseXMLObject(XMLObject xo) throws XMLParseException {


        Patterns pats = (Patterns)xo.getChild(Patterns.class);
        if(xo.hasAttribute(MODE)){
            return new MsatPatternStatistic(pats, xo.getStringAttribute(MODE));
        }


        return new MsatPatternStatistic(pats);

    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(Patterns.class),
                AttributeRule.newStringRule(MODE, true)
        };
    }



    public String getParserDescription(){
        return "Returns MsatPatternStatistic object";
    }

    public Class getReturnType(){
        return MsatPatternStatistic.class;
    }
}
