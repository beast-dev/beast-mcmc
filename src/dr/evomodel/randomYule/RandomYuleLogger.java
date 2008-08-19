package dr.evomodel.randomYule;

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
 * A logger that the indicator and rate vector for random local clock.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TreeLogger.java,v 1.25 2006/09/05 13:29:34 rambaut Exp $
 */
public class RandomYuleLogger extends MCLogger {

    private static final String RANDOM_YULE_LOGGER = "randomYuleLogger";
    private static final String FILENAME = "fileName";
    private static final String LOG_EVERY = "logEvery";

    private TreeModel treeModel;
    private RandomLocalYuleModel randomYule;

    public RandomYuleLogger(LogFormatter formatter, int logEvery, TreeModel treeModel, RandomLocalYuleModel randomYule) {

        super(formatter, logEvery, false);
        this.treeModel = treeModel;
        this.randomYule = randomYule;
    }

    public void startLogging() {

        super.startLogging();
        logLine("State\tRate changes");
    }

    public void log(int state) {

        if (logEvery <= 0 || ((state % logEvery) == 0)) {

            int nodeCount = treeModel.getNodeCount();

            StringBuilder builder = new StringBuilder();
            builder.append(state);
            for (int i = 0; i < nodeCount; i++) {

                NodeRef node = treeModel.getNode(i);

                if (randomYule.isRateChangeOnBranchAbove(treeModel, node)) {
                    builder.append("\t");
                    builder.append(node.getNumber());
                    builder.append("\t");
                    builder.append(randomYule.getBirthRate(treeModel, node));
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
            return RANDOM_YULE_LOGGER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            RandomLocalYuleModel randomYule = (RandomLocalYuleModel) xo.getChild(RandomLocalYuleModel.class);

            String fileName = xo.getStringAttribute(FILENAME);
            int logEvery = xo.getIntegerAttribute(LOG_EVERY);

            TabDelimitedFormatter formatter =
                    null;
            try {
                formatter = new TabDelimitedFormatter(new PrintWriter(new FileWriter(fileName)));
            } catch (IOException e) {
                e.printStackTrace();
            }

            return new RandomYuleLogger(formatter, logEvery, treeModel, randomYule);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A speciation model of a Yule process whose rate can evolve down the tree.";
        }

        public Class getReturnType() {
            return RandomLocalYuleModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(TreeModel.class),
                new ElementRule(RandomLocalYuleModel.class),
                AttributeRule.newIntegerRule(LOG_EVERY),
                AttributeRule.newStringRule(FILENAME),
        };
    };
}
