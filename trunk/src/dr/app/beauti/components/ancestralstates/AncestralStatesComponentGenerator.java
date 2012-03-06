package dr.app.beauti.components.ancestralstates;

import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.AbstractPartitionData;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.util.XMLWriter;
import dr.util.Attribute;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class AncestralStatesComponentGenerator extends BaseComponentGenerator {

    private final static boolean DEBUG = true;

    private final static String DNDS_LOG_SUFFIX = ".dNdS.log";
    private final static String STATE_LOG_SUFFIX = ".states.log";

    public AncestralStatesComponentGenerator(final BeautiOptions options) {
        super(options);
    }

    public boolean usesInsertionPoint(InsertionPoint point) {

        AncestralStatesComponentOptions component = (AncestralStatesComponentOptions) options
                .getComponentOptions(AncestralStatesComponentOptions.class);

        boolean reconstructAtNodes = false;
        boolean reconstructAtMRCA = false;
        boolean countingStates = false;
        boolean dNdSRobustCounting = false;

        for (AbstractPartitionData partition : options.getDataPartitions()) {
            if (component.reconstructAtNodes(partition)) reconstructAtNodes = true;
            if (component.reconstructAtMRCA(partition)) reconstructAtMRCA = true;
            if (component.isCountingStates(partition)) countingStates = true;
            if (countingStates && component.dNdSRobustCounting(partition)) dNdSRobustCounting = true;
        }

        if (!reconstructAtNodes && !reconstructAtMRCA && !countingStates && !dNdSRobustCounting) {
            return false;
        }

        switch (point) {
            case IN_FILE_LOG_PARAMETERS:
                return countingStates || dNdSRobustCounting;

            case IN_TREE_LIKELIHOOD:
                return countingStates;

            case IN_TREES_LOG:
                return reconstructAtNodes || countingStates || dNdSRobustCounting;

            case IN_OPERATORS:
                return dNdSRobustCounting;

            case AFTER_TREES_LOG:
                return reconstructAtMRCA || dNdSRobustCounting;

            case AFTER_MCMC:
                return dNdSRobustCounting;

            default:
                return false;
        }

    }// END: usesInsertionPoint

    protected void generate(InsertionPoint point, Object item, XMLWriter writer) {

        AncestralStatesComponentOptions component = (AncestralStatesComponentOptions) options
                .getComponentOptions(AncestralStatesComponentOptions.class);

        switch (point) {
            case IN_OPERATORS:
                writeCodonPartitionedRobustCounting(writer, component);
                break;
            case IN_TREE_LIKELIHOOD:
                writeCountingParameter(writer, (AbstractPartitionData)item);

            case IN_FILE_LOG_PARAMETERS:
                writeLogs(writer, component);
                break;
            case IN_TREES_LOG:
                writeTreeLogs(writer, component);
                break;
            case AFTER_TREES_LOG:
                writeAncestralStateLoggers(writer, component);
                break;
            case AFTER_MCMC:
                writeDNdSPerSiteAnalysisReport(writer, component);
                break;
            default:
                throw new IllegalArgumentException(
                        "This insertion point is not implemented for "
                                + this.getClass().getName());
        }

    }// END: generate

    private void writeCountingParameter(XMLWriter writer, AbstractPartitionData partition) {
        StringBuilder matrix = new StringBuilder();

        for (int i = 0; i < partition.getDataType().getStateCount(); i++) {
            for (int j = 0; j < partition.getDataType().getStateCount(); j++) {
                if (i == j) {
                    matrix.append(" 0.0");
                } else {
                    matrix.append(" 1.0");
                }
            }
        }
        writer.writeTag("parameter",
                new Attribute[] {
                        new Attribute.Default<String>("id", partition.getPrefix() + "count"),
                        new Attribute.Default<String>("value", matrix.toString()) },
                true);

  }

    protected String getCommentLabel() {
        return "Ancestral state reconstruction";
    }

    private void writeCodonPartitionedRobustCounting(XMLWriter writer,
                                                     AncestralStatesComponentOptions component) {

        for (AbstractPartitionData partition : options.getDataPartitions()) {

            if (component.dNdSRobustCounting(partition)) {
                writeCodonPartitionedRobustCounting(writer, partition);
            }
        }
    }

    // Called for each model that requires robust counting (can be more than
    // one)
    private void writeCodonPartitionedRobustCounting(XMLWriter writer,
                                                     AbstractPartitionData partition) {

        if (DEBUG) {
            System.err.println("DEBUG: Writing RB for " + partition.getName());
        }

        writer.writeComment("Robust counting for: " + partition.getName());

        // TODO: Hand coding is so 90s
        String prefix = partition.getName() + ".";

        // S operator
        writer.writeOpenTag("codonPartitionedRobustCounting",
                new Attribute[] {
                        new Attribute.Default<String>("id", prefix + "robustCounting1"),
                        new Attribute.Default<String>("labeling", "S"),
                        new Attribute.Default<String>("useUniformization",
                                "true"),
                        new Attribute.Default<String>("unconditionedPerBranch",
                                "true") });

        writer.writeIDref("treeModel", "treeModel");
        writer.writeOpenTag("firstPosition");
        writer.writeIDref("ancestralTreeLikelihood", partition.getPrefix() + "CP1.treeLikelihood");
        writer.writeCloseTag("firstPosition");

        writer.writeOpenTag("secondPosition");
        writer.writeIDref("ancestralTreeLikelihood", partition.getPrefix() + "CP2.treeLikelihood");
        writer.writeCloseTag("secondPosition");

        writer.writeOpenTag("thirdPosition");
        writer.writeIDref("ancestralTreeLikelihood", partition.getPrefix() + "CP3.treeLikelihood");
        writer.writeCloseTag("thirdPosition");

        writer.writeCloseTag("codonPartitionedRobustCounting");

        writer.writeBlankLine();

        // N operator:
        writer.writeOpenTag("codonPartitionedRobustCounting",
                new Attribute[] {
                        new Attribute.Default<String>("id", prefix + "robustCounting2"),
                        new Attribute.Default<String>("labeling", "N"),
                        new Attribute.Default<String>("useUniformization",
                                "true"),
                        new Attribute.Default<String>("unconditionedPerBranch",
                                "true") });

        writer.writeIDref("treeModel", "treeModel");
        writer.writeOpenTag("firstPosition");
        writer.writeIDref("ancestralTreeLikelihood", partition.getPrefix() + "CP1.treeLikelihood");
        writer.writeCloseTag("firstPosition");

        writer.writeOpenTag("secondPosition");
        writer.writeIDref("ancestralTreeLikelihood", partition.getPrefix() + "CP2.treeLikelihood");
        writer.writeCloseTag("secondPosition");

        writer.writeOpenTag("thirdPosition");
        writer.writeIDref("ancestralTreeLikelihood", partition.getPrefix() + "CP3.treeLikelihood");
        writer.writeCloseTag("thirdPosition");

        writer.writeCloseTag("codonPartitionedRobustCounting");

    }// END: writeCodonPartitionedRobustCounting()

    private void writeLogs(XMLWriter writer, AncestralStatesComponentOptions component) {

        for (AbstractPartitionData partition : options.getDataPartitions()) {

            if (component.dNdSRobustCounting(partition)) {
                String prefix = partition.getName() + ".";

                writer.writeIDref("codonPartitionedRobustCounting", prefix + "robustCounting1");
                writer.writeIDref("codonPartitionedRobustCounting", prefix + "robustCounting2");
            } else {

            }
        }
    }

    private void writeTreeLogs(XMLWriter writer, AncestralStatesComponentOptions component) {

        for (AbstractPartitionData partition : options.getDataPartitions()) {

            if (component.dNdSRobustCounting(partition)) {
                String prefix = partition.getName() + ".";

                writer.writeIDref("codonPartitionedRobustCounting", prefix + "robustCounting1");
                writer.writeIDref("codonPartitionedRobustCounting", prefix + "robustCounting2");
            }

            if (component.reconstructAtNodes(partition)) {
                writer.writeOpenTag("trait",
                        new Attribute[] {
                                new Attribute.Default<String>("name", "states"),
                                new Attribute.Default<String>("tag", partition.getPrefix() + "state")
                        }
                );
                writer.writeIDref("ancestralTreeLikelihood", partition.getPrefix() + "treeLikelihood");
                writer.writeCloseTag("trait");
            }

            if (component.isCountingStates(partition)) {
                writer.writeOpenTag("trait",
                        new Attribute[] {
                                new Attribute.Default<String>("name", partition.getPrefix() + "count"),
                                new Attribute.Default<String>("tag", partition.getPrefix() + "count")
                        }
                );
                writer.writeIDref("ancestralTreeLikelihood", partition.getPrefix() + "treeLikelihood");
                writer.writeCloseTag("trait");
            }
        }
    }

    private void writeAncestralStateLoggers(XMLWriter writer,
                                           AncestralStatesComponentOptions component) {

        for (AbstractPartitionData partition : options.getDataPartitions()) {

            if (component.dNdSRobustCounting(partition)) {
                writeDNdSLogger(writer, partition);
            } else if (component.reconstructAtMRCA(partition)) {
                writeStateLogger(writer, partition, component.getMRCATaxonSet(partition));
            }
        }
    }

    private void writeStateLogger(XMLWriter writer, AbstractPartitionData partition, String mrcaId) {

        writer.writeComment("Ancestral state reconstruction for: " + partition.getName());

        writer.writeOpenTag("log", new Attribute[] {
                new Attribute.Default<String>("id", "fileLog_" + partition.getName()),
                new Attribute.Default<String>("logEvery", Integer.toString(options.logEvery)),
                new Attribute.Default<String>("fileName", partition.getName() + STATE_LOG_SUFFIX) });

        writer.writeOpenTag("ancestralTrait",
                        new Attribute[]{new Attribute.Default<String>("name", partition.getName())});
        writer.writeIDref("treeModel", partition.getPartitionTreeModel().getPrefix() + "treeModel");
        writer.writeIDref("ancestralTreeLikelihood", partition.getPrefix() + "treeLikelihood");

        if (mrcaId != null) {
            writer.writeOpenTag("mrca");
            writer.writeIDref("taxa", mrcaId);
            writer.writeCloseTag("mrca");

        }
        writer.writeCloseTag("ancestralTrait");

        writer.writeCloseTag("log");

    }// END: writeLogs

    private void writeDNdSLogger(XMLWriter writer, AbstractPartitionData partition) {

        String prefix = partition.getName() + ".";

        writer.writeComment("Robust counting for: " + partition.getName());

        writer.writeOpenTag("log", new Attribute[] {
                new Attribute.Default<String>("id", "fileLog_dNdS"),
                new Attribute.Default<String>("logEvery", Integer.toString(options.logEvery)),
                new Attribute.Default<String>("fileName", partition.getName() + DNDS_LOG_SUFFIX) });

        writer
                .writeOpenTag("dNdSLogger",
                        new Attribute[]{new Attribute.Default<String>("id",
                                "dNdS")});
        writer.writeIDref("treeModel", "treeModel");
        writer.writeIDref("codonPartitionedRobustCounting", prefix + "robustCounting1");
        writer.writeIDref("codonPartitionedRobustCounting", prefix + "robustCounting2");
        writer.writeCloseTag("dNdSLogger");

        writer.writeCloseTag("log");

    }// END: writeLogs

    private void writeDNdSPerSiteAnalysisReport(XMLWriter writer,
                                                AncestralStatesComponentOptions component) {

        for (AbstractPartitionData partition : options.getDataPartitions()) {

            if (component.dNdSRobustCounting(partition)) {
                writeDNdSPerSiteAnalysisReport(writer, partition);
            }
        }
    }

    private void writeDNdSPerSiteAnalysisReport(XMLWriter writer,
                                                AbstractPartitionData partition) {

        writer.writeComment("Robust counting for: " + partition.getName());

        writer.writeOpenTag("report");

        writer.write("<dNdSPerSiteAnalysis fileName=" + '\"' + partition.getName() + DNDS_LOG_SUFFIX + '\"' + "/> \n");

        writer.writeCloseTag("report");

    }// END: writeDNdSReport

}
