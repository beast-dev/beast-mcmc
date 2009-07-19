package dr.evoxml;

import dr.evolution.datatype.HiddenNucleotides;
import dr.xml.*;

/**
 * @author Alexei Drummond
 */
public class HiddenNucleotideParser extends AbstractXMLObjectParser {

    public static final String HIDDEN_NUCLEOTIDES = "hiddenNucleotides";
    public static final String HIDDEN_CLASS_COUNT = "classCount";

    public String getParserName() {
        return HIDDEN_NUCLEOTIDES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int hiddenClassCount = xo.getIntegerAttribute(HIDDEN_CLASS_COUNT);
        return new HiddenNucleotides(hiddenClassCount);
    }

    public String getParserDescription() {
        return "A nucleotide data type that allows hidden substitution classes";
    }

    public Class getReturnType() {
        return HiddenNucleotides.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{AttributeRule.newIntegerRule(HIDDEN_CLASS_COUNT)};
    }
}