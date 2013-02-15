package dr.evomodelxml.tree;

import dr.evomodel.tree.randomlocalmodel.RLTVLoggerOnTree;
import dr.evomodel.tree.randomlocalmodel.RandomLocalTreeVariable;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class RLTVLoggerOnTreeParser extends AbstractXMLObjectParser {

    private static final String RANDOM_LOCAL_LOGGER = "randomLocalLoggerOnTree";

        public String getParserName() {
            return RANDOM_LOCAL_LOGGER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            RandomLocalTreeVariable treeVariable = (RandomLocalTreeVariable) xo.getChild(RandomLocalTreeVariable.class);
            return new RLTVLoggerOnTree(treeVariable);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A parser to log changed == 0 : 1 in a tree log";
        }

        public Class getReturnType() {
            return RLTVLoggerOnTree.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(RandomLocalTreeVariable.class),
        };
}
