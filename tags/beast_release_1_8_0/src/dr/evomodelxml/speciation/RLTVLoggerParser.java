package dr.evomodelxml.speciation;

import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.randomlocalmodel.RLTVLogger;
import dr.evomodel.tree.randomlocalmodel.RandomLocalTreeVariable;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.xml.*;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Parses an element from an DOM document into a SpeciationModel. Recognises YuleModel.
 */
public class RLTVLoggerParser extends AbstractXMLObjectParser {

    private static final String RANDOM_LOCAL_LOGGER = "randomLocalLogger";
    private static final String FILENAME = "fileName";
    private static final String LOG_EVERY = "logEvery";

        public String getParserName() {
            return RANDOM_LOCAL_LOGGER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            RandomLocalTreeVariable randomLocal =
                    (RandomLocalTreeVariable) xo.getChild(RandomLocalTreeVariable.class);

            String fileName = xo.getStringAttribute(FILENAME);
            int logEvery = xo.getIntegerAttribute(LOG_EVERY);

            TabDelimitedFormatter formatter = null;
            try {
                formatter = new TabDelimitedFormatter(new PrintWriter(new FileWriter(fileName)));
            } catch (IOException e) {
                e.printStackTrace();
            }

            return new RLTVLogger(formatter, logEvery, treeModel, randomLocal);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A speciation model of a Yule process whose rate can evolve down the tree.";
        }

        public Class getReturnType() {
            return RLTVLogger.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(TreeModel.class),
                new ElementRule(RandomLocalTreeVariable.class),
                AttributeRule.newIntegerRule(LOG_EVERY),
                AttributeRule.newStringRule(FILENAME),
        };
}
