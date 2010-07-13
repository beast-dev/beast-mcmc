package dr.app.beagle.evomodel.parsers;

import dr.app.beagle.evomodel.substmodel.CodonPartitionedRobustCounting;
import dr.app.beagle.evomodel.utilities.DnDsLogger;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.xml.*;

/**
 * @author Philippe Lemey
 * @author Marc A. Suchard
 */
public class DnDsLoggerParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "dNdSLogger";

//    public static final String COND_S = "condS";
//    public static final String COND_N = "condN";
//    public static final String UNCOND_S = "uncondS";
//    public static final String UNCOND_N = "uncondN";
//
//    public static final String SYN = "synonymous";
//    public static final String NON_SYN = "nonsynonymous";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String[] names = DnDsLogger.traitNames;
        TreeTrait[] foundTraits = new TreeTrait[names.length];

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object obj = xo.getChild(i);
            if (obj instanceof CodonPartitionedRobustCounting) {
                CodonPartitionedRobustCounting thisCount = (CodonPartitionedRobustCounting) obj;
                for (int j = 0; j < names.length; j++) {
                    TreeTrait trait = thisCount.getTreeTrait(names[j]);
                    if (trait != null) {
                        foundTraits[j] = trait;
                    }
                }
            }
        }

        for (int i = 0; i < foundTraits.length; i++) {
            if (foundTraits[i] == null) {
                throw new XMLParseException("Unable to find trait '" + names[i] + "'");
            }
        }

        Tree tree = (Tree) xo.getChild(Tree.class);

        // Use AttributeRules for options here

        return new DnDsLogger(xo.getId(), tree, foundTraits);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private static XMLSyntaxRule[] rules = {
            new ElementRule(CodonPartitionedRobustCounting.class, 2, 2),
            new ElementRule(Tree.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return DnDsLogger.class;
    }

    public String getParserName() {
        return PARSER_NAME;
    }
}
