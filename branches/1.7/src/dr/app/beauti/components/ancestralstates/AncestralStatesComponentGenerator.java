package dr.app.beauti.components.ancestralstates;

import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.generator.Generator;
import dr.app.beauti.options.AbstractPartitionData;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionData;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.ContinuousDataType;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.datatype.Nucleotides;
import dr.evomodelxml.treelikelihood.AncestralStateTreeLikelihoodParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
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
            if (component.dNdSRobustCounting(partition)) dNdSRobustCounting = true;
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

            case AFTER_OPERATORS:
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
            case AFTER_OPERATORS:
                writeCodonPartitionedRobustCounting(writer, component);
                break;
            case IN_TREE_LIKELIHOOD:
                writeCountingParameter(writer, (AbstractPartitionData)item);
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

    private void writeCountingParameter(XMLWriter writer, AbstractPartitionData partition) {
        AncestralStatesComponentOptions component = (AncestralStatesComponentOptions) options
                .getComponentOptions(AncestralStatesComponentOptions.class);
        if (!component.isCountingStates(partition)) {
            return;
        }

        StringBuilder matrix = new StringBuilder();


        DataType dataType = partition.getDataType();
        int stateCount = dataType.getStateCount();

        if (dataType == GeneralDataType.INSTANCE) {
            PartitionSubstitutionModel substModel = partition.getPartitionSubstitutionModel();
            stateCount = substModel.getDiscreteStateSet().size();
        }

        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                if (i == j) {
                    matrix.append(" 0.0");
                } else {
                    matrix.append(" 1.0");
                }
            }
        }
        writer.writeTag("parameter",
                new Attribute[] {
                        new Attribute.Default<String>("id", partition.getPrefix() + "count"),  // TODO Pass codon partition number, so can construct unique name
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

//        if (DEBUG) {
//            System.err.println("DEBUG: Writing RC for " + partition.getName());
//        }

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

                // is an alignment data partition
                PartitionSubstitutionModel substModel = partition.getPartitionSubstitutionModel();
                int cpCount = partition.getPartitionSubstitutionModel().getCodonPartitionCount();

                if (cpCount > 1) {
                    for (int i = 1; i <= substModel.getCodonPartitionCount(); i++) {
                        String prefix = partition.getPrefix() + substModel.getPrefix(i);
                        String name = partition.getName() + "." + substModel.getPrefix(i);
                        if (name.endsWith(".")) {
                            name = name.substring(0, name.length() - 1);
                        }
                        writeTrait(writer, partition, prefix, AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG, name);
                    }
                } else {
                    writeTrait(writer, partition, partition.getPrefix(), AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG, partition.getName());

                }
            }

            if (component.isCountingStates(partition)) {
                if (partition.getDataType() == ContinuousDataType.INSTANCE)  {
                    throw new RuntimeException("Can't do counting on Continuous data partition");
                }

                PartitionSubstitutionModel substModel = partition.getPartitionSubstitutionModel();
                int cpCount = partition.getPartitionSubstitutionModel().getCodonPartitionCount();

                if (cpCount > 1) {
                    for (int i = 1; i <= substModel.getCodonPartitionCount(); i++) {
                        String prefix = partition.getPrefix() + substModel.getPrefix(i);
                        String name = partition.getName() + "." + substModel.getPrefix(i) + "count";
                        writeTrait(writer, partition, prefix, "count", name);
                    }
                } else {
                    writeTrait(writer, partition, partition.getPrefix(), "count", partition.getName() + ".count");

                }

            }
        }
    }

    private void writeTrait(XMLWriter writer, AbstractPartitionData partition, String prefix, String tag, String name) {
        String traitName = prefix + tag;

        if (partition.getDataType() == ContinuousDataType.INSTANCE)  {
            traitName = partition.getName();
        }

        writer.writeOpenTag("trait",
                new Attribute[] {
                        new Attribute.Default<String>("name", traitName),
                        new Attribute.Default<String>("tag", name)
                }
        );
        if (partition.getDataType() == ContinuousDataType.INSTANCE)  {
            writer.writeIDref("multivariateTraitLikelihood", prefix + "traitLikelihood");
        } else {
            writer.writeIDref("ancestralTreeLikelihood", prefix + "treeLikelihood");
        }
        writer.writeCloseTag("trait");
    }

    private void writeAncestralStateLoggers(XMLWriter writer,
                                            AncestralStatesComponentOptions component) {

        for (AbstractPartitionData partition : options.getDataPartitions()) {

            if (component.dNdSRobustCounting(partition)) {
                writeDNdSLogger(writer, partition);
            }
            if (component.reconstructAtMRCA(partition)) {
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

        PartitionSubstitutionModel substModel = partition.getPartitionSubstitutionModel();
        int cpCount = partition.getPartitionSubstitutionModel().getCodonPartitionCount();

        if (cpCount > 1) {
            for (int i = 1; i <= substModel.getCodonPartitionCount(); i++) {
                String prefix = partition.getPrefix() + substModel.getPrefix(i);
                String name = partition.getName() + "." + substModel.getPrefix(i);
                if (name.endsWith(".")) {
                    name = name.substring(0, name.length() - 1);
                }
                writeAncestralTrait(writer, partition, mrcaId, prefix, name);
            }
        } else {
            writeAncestralTrait(writer, partition, mrcaId, partition.getPrefix(), partition.getName());
        }

        writer.writeCloseTag("log");

    }// END: writeLogs

    private void writeAncestralTrait(XMLWriter writer, AbstractPartitionData partition, String mrcaId, String prefix, String nameString) {
        String traitName = prefix + AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG;
        if (partition.getDataType() == ContinuousDataType.INSTANCE)  {
            traitName = partition.getName();
        }

        writer.writeOpenTag("ancestralTrait",
                new Attribute[]{
                        new Attribute.Default<String>("name", nameString),
                        new Attribute.Default<String>("traitName", traitName)
                }
        );
        writer.writeIDref("treeModel", partition.getPartitionTreeModel().getPrefix() + "treeModel");

        if (partition.getDataType() == ContinuousDataType.INSTANCE)  {
            writer.writeIDref("multivariateTraitLikelihood", prefix + "traitLikelihood");
        } else {
            writer.writeIDref("ancestralTreeLikelihood", prefix + "treeLikelihood");
        }

        if (mrcaId != null) {
            writer.writeOpenTag("mrca");
            writer.writeIDref("taxa", mrcaId);
            writer.writeCloseTag("mrca");
        }
        writer.writeCloseTag("ancestralTrait");
    }

    private void writeDNdSLogger(XMLWriter writer, AbstractPartitionData partition) {

        String prefix = partition.getName() + ".";

        writer.writeComment("Robust counting for: " + partition.getName());

        writer.writeOpenTag("log", new Attribute[] {
                new Attribute.Default<String>("id", "fileLog_dNdS"),
                new Attribute.Default<String>("logEvery", Integer.toString(options.logEvery)),
                new Attribute.Default<String>("fileName", partition.getName() + DNDS_LOG_SUFFIX) });

        writer.writeOpenTag("dNdSLogger", new Attribute[]{new Attribute.Default<String>("id",
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
