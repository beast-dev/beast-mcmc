package dr.evomodel.tree.randomlocalmodel;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.xml.*;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * A logger for a random local tree variable.
 * It logs only the selected parameters in pairs of node number, variable value
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TreeLogger.java,v 1.25 2006/09/05 13:29:34 rambaut Exp $
 */
public class RLTVLogger extends MCLogger {

    private static final String RANDOM_LOCAL_LOGGER = "randomLocalLogger";
    private static final String FILENAME = "fileName";
    private static final String LOG_EVERY = "logEvery";

    private TreeModel treeModel;
    private RandomLocalTreeVariable randomLocal;

    public RLTVLogger(LogFormatter formatter,
                      int logEvery, TreeModel treeModel,
                      RandomLocalTreeVariable randomLocal) {

        super(formatter, logEvery, false);
        this.treeModel = treeModel;
        this.randomLocal = randomLocal;
    }

    public void startLogging() {
        logLine("State\tRate changes");
    }

    public void log(int state) {

        if (logEvery <= 0 || ((state % logEvery) == 0)) {

            int nodeCount = treeModel.getNodeCount();

            StringBuilder builder = new StringBuilder();
            builder.append(state);
            for (int i = 0; i < nodeCount; i++) {

                NodeRef node = treeModel.getNode(i);

                if (randomLocal.isVariableSelected(treeModel, node)) {
                    builder.append("\t");
                    builder.append(node.getNumber());
                    builder.append("\t");
                    builder.append(randomLocal.getVariable(treeModel, node));
                }
            }

            logLine(builder.toString());
        }
    }

    public void stopLogging() {

        super.stopLogging();
    }

    /**
     * Parses an element from an DOM document into a SpeciationModel. Recognises
     * YuleModel.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return RANDOM_LOCAL_LOGGER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            RandomLocalTreeVariable randomLocal =
                    (RandomLocalTreeVariable) xo.getChild(RandomLocalTreeVariable.class);

            String fileName = xo.getStringAttribute(FILENAME);
            int logEvery = xo.getIntegerAttribute(LOG_EVERY);

            TabDelimitedFormatter formatter =
                    null;
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
    };
}
