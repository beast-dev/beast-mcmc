package dr.evomodel.epidemiology.casetocase;

import dr.evolution.tree.*;
import dr.evomodel.tree.TreeLogger;
import dr.evomodelxml.tree.TreeLoggerParser;
import dr.inference.loggers.LogFormatter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObjectParser;

import java.text.NumberFormat;

/**
 * Logs a partitioned tree. Breaks branches in two at infection times and adds a node of degree 2.
 * TreeAnnotator won't work.
 *
 * @author Matthew Hall
 */
public class PartitionedTreeLogger extends TreeLogger {

    CaseToCaseTreeLikelihood c2cTL;
    Tree originalTree;
    private LogUpon condition = null;
    private TreeTraitProvider.Helper treeTraitHelper;

    public PartitionedTreeLogger(CaseToCaseTreeLikelihood c2cTL, Tree tree, BranchRates branchRates,
                      TreeAttributeProvider[] treeAttributeProviders,
                      TreeTraitProvider[] treeTraitProviders,
                      LogFormatter formatter, int logEvery, boolean nexusFormat,
                      boolean sortTranslationTable, boolean mapNames, NumberFormat format,
                      TreeLogger.LogUpon condition) {
        super(tree, branchRates, treeAttributeProviders, treeTraitProviders, formatter, logEvery, nexusFormat,
                sortTranslationTable, mapNames, format, condition);
        this.c2cTL = c2cTL;
        this.originalTree = tree;
        this.condition = condition;
    }

    public void log(long state){

        final boolean doIt = condition != null ? condition.logNow(state) :
                (logEvery < 0 || ((state % logEvery) == 0));
        if(!doIt)
            return;

        setTree(c2cTL.rewireTree(originalTree));

        super.log(state);

    }



}
