package dr.app.beauti.components.ancestralstates;

import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.*;
import dr.app.beauti.util.XMLWriter;
import dr.util.Attribute;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class AncestralStatesComponentGenerator extends BaseComponentGenerator {

    private final static boolean DEBUG = true;

    private final static String DNDS_LOG_SUFFIX = ".dNdS.log";

    public AncestralStatesComponentGenerator(final BeautiOptions options) {
        super(options);
    }


    public boolean usesInsertionPoint(InsertionPoint point) {

        AncestralStatesComponentOptions component = (AncestralStatesComponentOptions) options
                .getComponentOptions(AncestralStatesComponentOptions.class);

        boolean reconstructAtNodes = false;
        boolean reconstructAtMRCA = false;
        boolean robustCounting = false;
        boolean dNdSRobustCounting = false;

        for (AbstractPartitionData partition : options.getDataPartitions()) {
            if (component.reconstructAtNodes(partition)) reconstructAtNodes = true;
            if (component.reconstructAtMRCA(partition)) reconstructAtMRCA = true;
            if (component.robustCounting(partition)) robustCounting = true;
            if (component.dNdSRobustCounting(partition)) dNdSRobustCounting = true;
        }

        if (!reconstructAtNodes && !reconstructAtMRCA && !robustCounting) {
            return false;
        }

        switch (point) {
            case IN_FILE_LOG_PARAMETERS:
                return true;

            case IN_TREES_LOG:
                return !reconstructAtMRCA;

            case IN_OPERATORS:
            case AFTER_TREES_LOG:
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
                writer.writeIDref("ancestralTreeLikelihood", partition.getPrefix() + "treeLikelihood");
            }
        }
    }

    private void writeTreeLogs(XMLWriter writer, AncestralStatesComponentOptions component) {

        for (AbstractPartitionData partition : options.getDataPartitions()) {

            if (component.dNdSRobustCounting(partition)) {
                String prefix = partition.getName() + ".";

                writer.writeIDref("codonPartitionedRobustCounting", prefix + "robustCounting1");
                writer.writeIDref("codonPartitionedRobustCounting", prefix + "robustCounting2");
            } else {
                writer.writeIDref("ancestralTreeLikelihood", partition.getPrefix() + "treeLikelihood");
            }
        }
    }

    private void writeAncestralStateLoggers(XMLWriter writer,
                                           AncestralStatesComponentOptions component) {

        for (AbstractPartitionData partition : options.getDataPartitions()) {

            if (component.dNdSRobustCounting(partition)) {
                writeDNdSLogger(writer, partition);
            }
        }
    }

    private void writeDNdSLogger(XMLWriter writer, AbstractPartitionData partition) {

        String prefix = partition.getName() + ".";

        writer.writeComment("Robust counting for: " + partition.getName());

        writer.writeOpenTag("log", new Attribute[] {
                new Attribute.Default<String>("id", "fileLog_dNdS"),
                new Attribute.Default<String>("logEvery", "10000"),
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
