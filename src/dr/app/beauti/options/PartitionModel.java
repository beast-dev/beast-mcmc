package dr.app.beauti.options;

import dr.evolution.datatype.DataType;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public class PartitionModel extends AbstractModelOptions {

    public PartitionModel(DataPartition partition) {
        this(partition.getName(), partition.getAlignment().getDataType());
    }

    public String getName() {
        return name;
    }

    public PartitionModel(String name, DataType dataType) {

        this.name = name;
        this.dataType = dataType;

        double substWeights = 1.0;
        double branchWeights = 30.0;
        double treeWeights = 15.0;

        //Substitution model parameters
        createParameter("hky.frequencies", "HKY base frequencies", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("hky1.frequencies", "HKY base frequencies for codon position 1", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("hky2.frequencies", "HKY base frequencies for codon position 2", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("hky3.frequencies", "HKY base frequencies for codon position 3", UNITY_SCALE, 0.25, 0.0, 1.0);

        createScaleParameter("hky.kappa", "HKY transition-transversion parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("hky1.kappa", "HKY transition-transversion parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("hky2.kappa", "HKY transition-transversion parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("hky3.kappa", "HKY transition-transversion parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        createParameter("gtr.frequencies", "GTR base frequencies", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("gtr1.frequencies", "GTR base frequencies for codon position 1", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("gtr2.frequencies", "GTR base frequencies for codon position 2", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("gtr3.frequencies", "GTR base frequencies for codon position 3", UNITY_SCALE, 0.25, 0.0, 1.0);

        createScaleParameter("gtr.ac", "GTR A-C substitution parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr.ag", "GTR A-G substitution parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr.at", "GTR A-T substitution parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr.cg", "GTR C-G substitution parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr.gt", "GTR G-T substitution parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        createScaleParameter("gtr1.ac", "GTR A-C substitution parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr1.ag", "GTR A-G substitution parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr1.at", "GTR A-T substitution parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr1.cg", "GTR C-G substitution parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr1.gt", "GTR G-T substitution parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        createScaleParameter("gtr2.ac", "GTR A-C substitution parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr2.ag", "GTR A-G substitution parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr2.at", "GTR A-T substitution parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr2.cg", "GTR C-G substitution parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr2.gt", "GTR G-T substitution parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        createScaleParameter("gtr3.ac", "GTR A-C substitution parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr3.ag", "GTR A-G substitution parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr3.at", "GTR A-T substitution parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr3.cg", "GTR C-G substitution parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr3.gt", "GTR G-T substitution parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        createParameter("bsimple.frequencies", "Binary Simple frequencies", UNITY_SCALE, 0.5, 0.0, 1.0);

        createParameter("bcov.frequencies", "Binary Covarion frequencies of the visible states", UNITY_SCALE, 0.5, 0.0, 1.0);
        createParameter("bcov.hfrequencies", "Binary Covarion frequencies of the hidden rates", UNITY_SCALE, 0.5, 0.0, 1.0);
        createParameter("bcov.alpha", "Binary Covarion rate of evolution in slow mode", UNITY_SCALE, 0.5, 0.0, 1.0);
        createParameter("bcov.s", "Binary Covarion rate of flipping between slow and fast modes", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, 100.0);

        createParameter("siteModel.alpha", "gamma shape parameter", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter("siteModel1.alpha", "gamma shape parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter("siteModel2.alpha", "gamma shape parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter("siteModel3.alpha", "gamma shape parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);

        createParameter("siteModel.pInv", "proportion of invariant sites parameter", NONE, 0.5, 0.0, 1.0);
        createParameter("siteModel1.pInv", "proportion of invariant sites parameter for codon position 1", NONE, 0.5, 0.0, 1.0);
        createParameter("siteModel2.pInv", "proportion of invariant sites parameter for codon position 2", NONE, 0.5, 0.0, 1.0);
        createParameter("siteModel3.pInv", "proportion of invariant sites parameter for codon position 3", NONE, 0.5, 0.0, 1.0);

        createParameter("siteModel1.mu", "relative rate parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("siteModel2.mu", "relative rate parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("siteModel3.mu", "relative rate parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("allMus", "All the relative rates");

        createOperator("hky.kappa", SCALE, 0.75, substWeights);
        createOperator("hky1.kappa", SCALE, 0.75, substWeights);
        createOperator("hky2.kappa", SCALE, 0.75, substWeights);
        createOperator("hky3.kappa", SCALE, 0.75, substWeights);
        createOperator("hky.frequencies", DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("hky1.frequencies", DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("hky2.frequencies", DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("hky3.frequencies", DELTA_EXCHANGE, 0.01, substWeights);

        createOperator("gtr.ac", SCALE, 0.75, substWeights);
        createOperator("gtr.ag", SCALE, 0.75, substWeights);
        createOperator("gtr.at", SCALE, 0.75, substWeights);
        createOperator("gtr.cg", SCALE, 0.75, substWeights);
        createOperator("gtr.gt", SCALE, 0.75, substWeights);

        createOperator("gtr1.ac", SCALE, 0.75, substWeights);
        createOperator("gtr1.ag", SCALE, 0.75, substWeights);
        createOperator("gtr1.at", SCALE, 0.75, substWeights);
        createOperator("gtr1.cg", SCALE, 0.75, substWeights);
        createOperator("gtr1.gt", SCALE, 0.75, substWeights);

        createOperator("gtr2.ac", SCALE, 0.75, substWeights);
        createOperator("gtr2.ag", SCALE, 0.75, substWeights);
        createOperator("gtr2.at", SCALE, 0.75, substWeights);
        createOperator("gtr2.cg", SCALE, 0.75, substWeights);
        createOperator("gtr2.gt", SCALE, 0.75, substWeights);

        createOperator("gtr3.ac", SCALE, 0.75, substWeights);
        createOperator("gtr3.ag", SCALE, 0.75, substWeights);
        createOperator("gtr3.at", SCALE, 0.75, substWeights);
        createOperator("gtr3.cg", SCALE, 0.75, substWeights);
        createOperator("gtr3.gt", SCALE, 0.75, substWeights);

        createOperator("gtr.frequencies", DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("gtr1.frequencies", DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("gtr2.frequencies", DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("gtr3.frequencies", DELTA_EXCHANGE, 0.01, substWeights);

        createOperator("bcov.alpha", SCALE, 0.75, substWeights);
        createOperator("bcov.s", SCALE, 0.75, substWeights);
        createOperator("bcov.frequencies", DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("bcov.hfrequencies", DELTA_EXCHANGE, 0.01, substWeights);

        createOperator("siteModel.alpha", SCALE, 0.75, substWeights);
        createOperator("siteModel1.alpha", SCALE, 0.75, substWeights);
        createOperator("siteModel2.alpha", SCALE, 0.75, substWeights);
        createOperator("siteModel3.alpha", SCALE, 0.75, substWeights);

        createOperator("siteModel.pInv", SCALE, 0.75, substWeights);
        createOperator("siteModel1.pInv", SCALE, 0.75, substWeights);
        createOperator("siteModel2.pInv", SCALE, 0.75, substWeights);
        createOperator("siteModel3.pInv", SCALE, 0.75, substWeights);

        createOperator("centeredMu", "Relative rates", "Scales codon position rates relative to each other maintaining mean", "allMus", CENTERED_SCALE, 0.75, substWeights);
        createOperator("deltaMu", "Relative rates", "Changes codon position rates relative to each other maintaining mean", "allMus", DELTA_EXCHANGE, 0.75, substWeights);
    }

    private void selectOperators(ArrayList<Operator> ops) {

        switch (dataType.getType()) {
            case DataType.NUCLEOTIDES:

                switch (nucSubstitutionModel) {
                    case HKY:
                        if (codonPartitionCount > 1 && unlinkedSubstitutionModel) {
                            for (int i = 1; i <= codonPartitionCount; i++) {
                                ops.add(getOperator("hky" + i + ".kappa"));
                            }
                        } else {
                            ops.add(getOperator("hky.kappa"));
                        }
                        if (frequencyPolicy == BeautiOptions.ESTIMATED) {
                            if (codonPartitionCount > 1 && unlinkedSubstitutionModel) {
                                for (int i = 1; i <= codonPartitionCount; i++) {
                                    ops.add(getOperator("hky" + i + ".frequencies"));
                                }
                            } else {
                                ops.add(getOperator("hky.frequencies"));
                            }
                        }
                        break;

                    case GTR:
                        //if (frequencyPolicy == BeautiOptions.ESTIMATED || frequencyPolicy == BeautiOptions.EMPIRICAL){
                        if (codonPartitionCount > 1 && unlinkedSubstitutionModel) {
                            for (int i = 1; i <= codonPartitionCount; i++) {
                                ops.add(getOperator("gtr" + i + ".ac"));
                                ops.add(getOperator("gtr" + i + ".ag"));
                                ops.add(getOperator("gtr" + i + ".at"));
                                ops.add(getOperator("gtr" + i + ".cg"));
                                ops.add(getOperator("gtr" + i + ".gt"));
                            }
                        } else {
                            ops.add(getOperator("gtr.ac"));
                            ops.add(getOperator("gtr.ag"));
                            ops.add(getOperator("gtr.at"));
                            ops.add(getOperator("gtr.cg"));
                            ops.add(getOperator("gtr.gt"));
                        }
                        //}

                        if (frequencyPolicy == BeautiOptions.ESTIMATED) {
                            if (codonPartitionCount > 1 && unlinkedSubstitutionModel) {
                                for (int i = 1; i <= codonPartitionCount; i++) {
                                    ops.add(getOperator("gtr" + i + ".frequencies"));
                                }
                            } else {
                                ops.add(getOperator("gtr.frequencies"));
                            }
                        }
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown nucleotides substitution model");
                }

                break;

            case DataType.AMINO_ACIDS:
                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:
                switch (binarySubstitutionModel) {
                    case BIN_SIMPLE:
                        break;

                    case BIN_COVARION:
                        ops.add(getOperator("bcov.alpha"));
                        ops.add(getOperator("bcov.s"));
                        ops.add(getOperator("bcov.frequencies"));
                        ops.add(getOperator("bcov.hfrequencies"));
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown binary substitution model");
                }
                break;

            default:
                throw new IllegalArgumentException("Unknown data type");
        }

        // if gamma do shape move
        if (gammaHetero) {
            if (codonPartitionCount > 1 && unlinkedHeterogeneityModel) {
                for (int i = 1; i <= codonPartitionCount; i++) {
                    ops.add(getOperator("siteModel" + i + ".alpha"));
                }
            } else {
                ops.add(getOperator("siteModel.alpha"));
            }
        }
        // if pinv do pinv move
        if (invarHetero) {
            if (codonPartitionCount > 1 && unlinkedHeterogeneityModel) {
                for (int i = 1; i <= codonPartitionCount; i++) {
                    ops.add(getOperator("siteModel" + i + ".pInv"));
                }
            } else {
                ops.add(getOperator("siteModel.pInv"));
            }
        }

        if (codonPartitionCount > 1) {
            if (!codonHeteroPattern.equals("112")) {
                ops.add(getOperator("centeredMu"));
            }
            ops.add(getOperator("deltaMu"));
        }
    }


    /**
     * return a list of parameters that are required
     */
    List<Parameter> getParameters() {

        List<Parameter> params = new ArrayList<Parameter>();

        if (codonPartitionCount > 1) {
            for (int i = 1; i <= codonPartitionCount; i++) {
                params.add(getParameter("siteModel" + i + ".mu"));
            }
        }
        switch (dataType.getType()) {
            case DataType.NUCLEOTIDES:
                switch (nucSubstitutionModel) {
                    case HKY:
                        if (codonPartitionCount > 1 && unlinkedSubstitutionModel) {
                            for (int i = 1; i <= codonPartitionCount; i++) {
                                params.add(getParameter("hky" + i + ".kappa"));
                            }
                        } else {
                            params.add(getParameter("hky.kappa"));
                        }
                        break;
                    case GTR:
                        if (codonPartitionCount > 1 && unlinkedSubstitutionModel) {
                            for (int i = 1; i <= codonPartitionCount; i++) {
                                params.add(getParameter("gtr" + i + ".ac"));
                                params.add(getParameter("gtr" + i + ".ag"));
                                params.add(getParameter("gtr" + i + ".at"));
                                params.add(getParameter("gtr" + i + ".cg"));
                                params.add(getParameter("gtr" + i + ".gt"));
                            }
                        } else {
                            params.add(getParameter("gtr.ac"));
                            params.add(getParameter("gtr.ag"));
                            params.add(getParameter("gtr.at"));
                            params.add(getParameter("gtr.cg"));
                            params.add(getParameter("gtr.gt"));
                        }
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown nucleotides substitution model");
                }
                break;

            case DataType.AMINO_ACIDS:
                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:
                switch (binarySubstitutionModel) {
                    case BIN_SIMPLE:
                        break;

                    case BIN_COVARION:
                        params.add(getParameter("bcov.alpha"));
                        params.add(getParameter("bcov.s"));
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown binary substitution model");
                }
                break;

            default:
                throw new IllegalArgumentException("Unknown data type");
        }

        // if gamma do shape move
        if (gammaHetero) {
            if (codonPartitionCount > 1 && unlinkedHeterogeneityModel) {
                for (int i = 1; i <= codonPartitionCount; i++) {
                    params.add(getParameter("siteModel" + i + ".alpha"));
                }
            } else {
                params.add(getParameter("siteModel.alpha"));
            }
        }
        // if pinv do pinv move
        if (invarHetero) {
            if (codonPartitionCount > 1 && unlinkedHeterogeneityModel) {
                for (int i = 1; i <= codonPartitionCount; i++) {
                    params.add(getParameter("siteModel" + i + ".pInv"));
                }
            } else {
                params.add(getParameter("siteModel.pInv"));
            }
        }

        return params;
    }

    protected Parameter getParameter(String name) {
        Parameter parameter = parameters.get(name);
        parameter.setPrefix(getName());

        if (parameter == null)
            throw new IllegalArgumentException("Parameter with name, " + parameter.getName() + ", is unknown");
        return parameter;
    }


    public int nucSubstitutionModel = BeautiOptions.HKY;
    public int aaSubstitutionModel = BeautiOptions.BLOSUM_62;
    public int binarySubstitutionModel = BeautiOptions.BIN_SIMPLE;

    public int frequencyPolicy = BeautiOptions.ESTIMATED;
    public boolean gammaHetero = false;
    public int gammaCategories = 4;
    public boolean invarHetero = false;
    public String codonHeteroPattern = null;
    public double meanSubstitutionRate = 1.0;
    public boolean unlinkedSubstitutionModel = false;
    public boolean unlinkedHeterogeneityModel = false;
    public boolean unlinkedFrequencyModel = false;

    public int codonPartitionCount;
    boolean fixedSubstitutionRate = false;

    public DataType dataType;
    public String name;
}
