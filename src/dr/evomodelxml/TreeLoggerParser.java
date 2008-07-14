package dr.evomodelxml;

import dr.evolution.colouring.BranchColouring;
import dr.evolution.colouring.TreeColouring;
import dr.evolution.colouring.TreeColouringProvider;
import dr.evolution.tree.*;
import dr.evomodel.tree.TreeLogger;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.model.Likelihood;
import dr.xml.*;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class TreeLoggerParser extends LoggerParser {

    public static final String LOG_TREE = "logTree";
    public static final String NEXUS_FORMAT = "nexusFormat";
    public static final String USING_RATES = "usingRates";
    public static final String BRANCH_LENGTHS = "branchLengths";
    public static final String TIME = "time";
    public static final String SUBSTITUTIONS = "substitutions";
    public static final String SORT_TRANSLATION_TABLE = "sortTranslationTable";
    public static final String MAP_NAMES = "mapNamesToNumbers";

    public String getParserName() {
        return LOG_TREE;
    }

    /**
     * @return an object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final Tree tree = (Tree) xo.getChild(Tree.class);

        String title = xo.getAttribute(TITLE, "");

        final boolean nexusFormat = xo.getAttribute(NEXUS_FORMAT, false);

        final boolean sortTranslationTable = xo.getAttribute(SORT_TRANSLATION_TABLE, true);

        boolean substitutions = xo.getAttribute(BRANCH_LENGTHS, "").equals(SUBSTITUTIONS);

        List<TreeAttributeProvider> taps = new ArrayList<TreeAttributeProvider>();
        List<NodeAttributeProvider> naps = new ArrayList<NodeAttributeProvider>();
        List<BranchAttributeProvider> baps = new ArrayList<BranchAttributeProvider>();

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object cxo = xo.getChild(i);
            if (cxo instanceof TreeColouringProvider) {
                final TreeColouringProvider colouringProvider = (TreeColouringProvider) cxo;
                baps.add(new BranchAttributeProvider() {

                    public String getBranchAttributeLabel() {
                        return "deme";
                    }

                    public String getAttributeForBranch(Tree tree, NodeRef node) {
                        TreeColouring colouring = colouringProvider.getTreeColouring(tree);
                        BranchColouring bcol = colouring.getBranchColouring(node);
                        StringBuilder buffer = new StringBuilder();
                        if (bcol != null) {
                            buffer.append("{");
                            buffer.append(bcol.getChildColour());
                            for (int i = 1; i <= bcol.getNumEvents(); i++) {
                                buffer.append(",");
                                buffer.append(bcol.getBackwardTime(i));
                                buffer.append(",");
                                buffer.append(bcol.getBackwardColourAbove(i));
                            }
                            buffer.append("}");
                        }
                        return buffer.toString();
                    }
                });

            } else if (cxo instanceof Likelihood) {
                final Likelihood likelihood = (Likelihood) cxo;
                taps.add(new TreeAttributeProvider() {

                    public String getTreeAttributeLabel() {
                        return "lnP";
                    }

                    public String getAttributeForTree(Tree tree) {
                        return Double.toString(likelihood.getLogLikelihood());
                    }
                });

            } //else {
            if (cxo instanceof TreeAttributeProvider) {
                taps.add((TreeAttributeProvider) cxo);
            }
            if (cxo instanceof NodeAttributeProvider) {
                naps.add((NodeAttributeProvider) cxo);
            }
            if (cxo instanceof BranchAttributeProvider) {
                baps.add((BranchAttributeProvider) cxo);
            }
            //}
        }
        BranchRateController branchRateProvider = null;
        if (substitutions) {
            branchRateProvider = (BranchRateController) xo.getChild(BranchRateController.class);
        }
        if (substitutions && branchRateProvider == null) {
            throw new XMLParseException("To log trees in units of substitutions a BranchRateModel must be provided");
        }

        // logEvery of zero only displays at the end
        int logEvery = xo.getAttribute(LOG_EVERY, 1);

        PrintWriter pw = getLogFile(xo, getParserName());

        LogFormatter formatter = new TabDelimitedFormatter(pw);

        TreeAttributeProvider[] treeAttributeProviders = new TreeAttributeProvider[taps.size()];
        taps.toArray(treeAttributeProviders);
        NodeAttributeProvider[] nodeAttributeProviders = new NodeAttributeProvider[naps.size()];
        naps.toArray(nodeAttributeProviders);
        BranchAttributeProvider[] branchAttributeProviders = new BranchAttributeProvider[baps.size()];
        baps.toArray(branchAttributeProviders);

        // I think the default should be to have names rather than numbers, thus the false default - AJD
        boolean mapNames = xo.getAttribute(MAP_NAMES, false);

        TreeLogger logger =
                new TreeLogger(tree, branchRateProvider,
                        treeAttributeProviders, nodeAttributeProviders, branchAttributeProviders,
                        formatter, logEvery, nexusFormat, sortTranslationTable, mapNames);

        if (title != null) {
            logger.setTitle(title);
        }

        return logger;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newIntegerRule(LOG_EVERY),
            new StringAttributeRule(FILE_NAME,
                    "The name of the file to send log output to. " +
                            "If no file name is specified then log is sent to standard output", true),
            new StringAttributeRule(TITLE, "The title of the log", true),
            AttributeRule.newBooleanRule(NEXUS_FORMAT, true,
                    "Whether to use the NEXUS format for the tree log"),
            AttributeRule.newBooleanRule(SORT_TRANSLATION_TABLE, true,
                    "Whether the translation table is sorted."),
            new StringAttributeRule(BRANCH_LENGTHS, "What units should the branch lengths be in",
                    new String[]{TIME, SUBSTITUTIONS}, true),
            new ElementRule(Tree.class, "The tree which is to be logged"),
            new ElementRule(BranchRateController.class, true),
            new ElementRule(TreeColouringProvider.class, true),
            new ElementRule(Likelihood.class, true),
            new ElementRule(TreeAttributeProvider.class, 0, Integer.MAX_VALUE),
            new ElementRule(NodeAttributeProvider.class, 0, Integer.MAX_VALUE),
            new ElementRule(BranchAttributeProvider.class, 0, Integer.MAX_VALUE),
            AttributeRule.newBooleanRule(MAP_NAMES, true)
    };

    public String getParserDescription() {
        return "Logs a tree to a file";
    }

    public String getExample() {
        final String name = getParserName();
        return
                "<!-- The " + name + " element takes a treeModel to be logged -->\n" +
                        "<" + name + " " + LOG_EVERY + "=\"100\" " + FILE_NAME + "=\"log.trees\" "
                        + NEXUS_FORMAT + "=\"true\">\n" +
                        "	<treeModel idref=\"treeModel1\"/>\n" +
                        "</" + name + ">\n";
    }

    public Class getReturnType() {
        return TreeLogger.class;
    }
}
