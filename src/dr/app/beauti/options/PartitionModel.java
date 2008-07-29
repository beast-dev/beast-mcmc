package dr.app.beauti.options;

import dr.evolution.datatype.DataType;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class PartitionModel extends ModelOptions {

    public PartitionModel(BeautiOptions options, DataPartition partition) {
        this(options, partition.getName(), partition.getAlignment().getDataType());
    }

    /**
     * A copy constructor
     *
     * @param options the beauti options
     * @param name    the name of the new model
     * @param source  the source model
     */
    public PartitionModel(BeautiOptions options, String name, PartitionModel source) {
        this(options, name, source.dataType);

        nucSubstitutionModel = source.nucSubstitutionModel;
        aaSubstitutionModel = source.aaSubstitutionModel;
        binarySubstitutionModel = source.binarySubstitutionModel;

        frequencyPolicy = source.frequencyPolicy;
        gammaHetero = source.gammaHetero;
        gammaCategories = source.gammaCategories;
        invarHetero = source.invarHetero;
        codonHeteroPattern = source.codonHeteroPattern;
        unlinkedSubstitutionModel = source.unlinkedSubstitutionModel;
        unlinkedHeterogeneityModel = source.unlinkedHeterogeneityModel;
        unlinkedFrequencyModel = source.unlinkedFrequencyModel;
    }

    public PartitionModel(BeautiOptions options, String name, DataType dataType) {

        this.options = options;
        this.name = name;
        this.dataType = dataType;

        double substWeights = 1.0;

        //Substitution model parameters
        createParameter("frequencies", "HKY base frequencies", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("CP1.frequencies", "HKY base frequencies for codon position 1", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("CP2.frequencies", "HKY base frequencies for codon position 2", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("CP1+2.frequencies", "HKY base frequencies for codon positions 1 & 2", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("CP3.frequencies", "HKY base frequencies for codon position 3", UNITY_SCALE, 0.25, 0.0, 1.0);

        createScaleParameter("kappa", "HKY transition-transversion parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP1.kappa", "HKY transition-transversion parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP2.kappa", "HKY transition-transversion parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP1+2.kappa", "HKY transition-transversion parameter for codon positions 1 & 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP3.kappa", "HKY transition-transversion parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        createParameter("frequencies", "GTR base frequencies", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("CP1.frequencies", "GTR base frequencies for codon position 1", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("CP2.frequencies", "GTR base frequencies for codon position 2", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("CP1+2.frequencies", "GTR base frequencies for codon positions 1 & 2", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("CP3.frequencies", "GTR base frequencies for codon position 3", UNITY_SCALE, 0.25, 0.0, 1.0);

        createScaleParameter("ac", "GTR A-C substitution parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("ag", "GTR A-G substitution parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("at", "GTR A-T substitution parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("cg", "GTR C-G substitution parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gt", "GTR G-T substitution parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        createScaleParameter("CP1.ac", "GTR A-C substitution parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP1.ag", "GTR A-G substitution parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP1.at", "GTR A-T substitution parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP1.cg", "GTR C-G substitution parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP1.gt", "GTR G-T substitution parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        createScaleParameter("CP2.ac", "GTR A-C substitution parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP2.ag", "GTR A-G substitution parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP2.at", "GTR A-T substitution parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP2.cg", "GTR C-G substitution parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP2.gt", "GTR G-T substitution parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        createScaleParameter("CP1+2.ac", "GTR A-C substitution parameter for codon positions 1 & 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP1+2.ag", "GTR A-G substitution parameter for codon positions 1 & 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP1+2.at", "GTR A-T substitution parameter for codon positions 1 & 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP1+2.cg", "GTR C-G substitution parameter for codon positions 1 & 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP1+2.gt", "GTR G-T substitution parameter for codon positions 1 & 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        createScaleParameter("CP3.ac", "GTR A-C substitution parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP3.ag", "GTR A-G substitution parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP3.at", "GTR A-T substitution parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP3.cg", "GTR C-G substitution parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP3.gt", "GTR G-T substitution parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        createParameter("frequencies", "Binary Simple frequencies", UNITY_SCALE, 0.5, 0.0, 1.0);

        createParameter("frequencies", "Binary Covarion frequencies of the visible states", UNITY_SCALE, 0.5, 0.0, 1.0);
        createParameter("hfrequencies", "Binary Covarion frequencies of the hidden rates", UNITY_SCALE, 0.5, 0.0, 1.0);
        createParameter("bcov.alpha", "Binary Covarion rate of evolution in slow mode", UNITY_SCALE, 0.5, 0.0, 1.0);
        createParameter("bcov.s", "Binary Covarion rate of flipping between slow and fast modes", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, 100.0);

        createParameter("alpha", "gamma shape parameter", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter("CP1.alpha", "gamma shape parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter("CP2.alpha", "gamma shape parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter("CP1+2.alpha", "gamma shape parameter for codon positions 1 & 2", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter("CP3.alpha", "gamma shape parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);

        createParameter("pInv", "proportion of invariant sites parameter", NONE, 0.5, 0.0, 1.0);
        createParameter("CP1.pInv", "proportion of invariant sites parameter for codon position 1", NONE, 0.5, 0.0, 1.0);
        createParameter("CP2.pInv", "proportion of invariant sites parameter for codon position 2", NONE, 0.5, 0.0, 1.0);
        createParameter("CP1+2.pInv", "proportion of invariant sites parameter for codon positions 1 & 2", NONE, 0.5, 0.0, 1.0);
        createParameter("CP3.pInv", "proportion of invariant sites parameter for codon position 3", NONE, 0.5, 0.0, 1.0);

        createParameter("mu", "relative rate parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("CP1.mu", "relative rate parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("CP2.mu", "relative rate parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("CP1+2.mu", "relative rate parameter for codon positions 1 & 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("CP3.mu", "relative rate parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

        createScaleOperator("kappa", substWeights);
        createScaleOperator("CP1.kappa", substWeights);
        createScaleOperator("CP2.kappa", substWeights);
        createScaleOperator("CP1+2.kappa", substWeights);
        createScaleOperator("CP3.kappa", substWeights);

        createOperator("frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("CP1.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("CP2.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("CP1+2.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("CP3.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);

        createScaleOperator("ac", substWeights);
        createScaleOperator("ag", substWeights);
        createScaleOperator("at", substWeights);
        createScaleOperator("cg", substWeights);
        createScaleOperator("gt", substWeights);

        for (int i = 1; i <= 3; i++) {
            createScaleOperator("CP" + i + ".ac", substWeights);
            createScaleOperator("CP" + i + ".ag", substWeights);
            createScaleOperator("CP" + i + ".at", substWeights);
            createScaleOperator("CP" + i + ".cg", substWeights);
            createScaleOperator("CP" + i + ".gt", substWeights);
        }

        createScaleOperator("CP1+2.ac", substWeights);
        createScaleOperator("CP1+2.ag", substWeights);
        createScaleOperator("CP1+2.at", substWeights);
        createScaleOperator("CP1+2.cg", substWeights);
        createScaleOperator("CP1+2.gt", substWeights);

//      These already exist, above
//        createOperator("gtr.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
//        createOperator("gtr1.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
//        createOperator("gtr2.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
//        createOperator("gtr3.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
//
        createScaleOperator("alpha", substWeights);
        for (int i = 1; i <= 3; i++) {
            createScaleOperator("CP" + i + ".alpha", substWeights);
        }
        createScaleOperator("CP1+2.alpha", substWeights);

        createScaleOperator("pInv", substWeights);
        for (int i = 1; i <= 3; i++) {
            createScaleOperator("CP" + i + ".pInv", substWeights);
        }
        createScaleOperator("CP1+2.pInv", substWeights);

        createScaleOperator("bcov.alpha", substWeights);
        createScaleOperator("bcov.s", substWeights);
        //createOperator("frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("hfrequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);

    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public List<Operator> getOperators() {
        List<Operator> operators = new ArrayList<Operator>();

        switch (dataType.getType()) {
            case DataType.NUCLEOTIDES:

                if (getCodonPartitionCount() > 1 && unlinkedSubstitutionModel) {
                    if (codonHeteroPattern.equals("123")) {
                        switch (nucSubstitutionModel) {
                            case HKY:
                                operators.add(getOperator("CP1.kappa"));
                                operators.add(getOperator("CP2.kappa"));
                                operators.add(getOperator("CP3.kappa"));
                                break;

                            case GTR:
                                for (int i = 1; i <= getCodonPartitionCount(); i++) {
                                    operators.add(getOperator("CP" + i + ".ac"));
                                    operators.add(getOperator("CP" + i + ".ag"));
                                    operators.add(getOperator("CP" + i + ".at"));
                                    operators.add(getOperator("CP" + i + ".cg"));
                                    operators.add(getOperator("CP" + i + ".gt"));
                                }
                                break;

                            default:
                                throw new IllegalArgumentException("Unknown nucleotides substitution model");
                        }

                        if (frequencyPolicy == FrequencyPolicy.ESTIMATED) {
                            if (getCodonPartitionCount() > 1 && unlinkedSubstitutionModel) {
                                operators.add(getOperator("CP1.frequencies"));
                                operators.add(getOperator("CP2.frequencies"));
                                operators.add(getOperator("CP3.frequencies"));
                            } else {
                                operators.add(getOperator("frequencies"));
                            }
                        }
                    } else if (codonHeteroPattern.equals("112")) {
                        switch (nucSubstitutionModel) {
                            case HKY:
                                operators.add(getOperator("CP1+2.kappa"));
                                operators.add(getOperator("CP3.kappa"));
                                break;

                            case GTR:
                                operators.add(getOperator("CP1+2.ac"));
                                operators.add(getOperator("CP1+2.ag"));
                                operators.add(getOperator("CP1+2.at"));
                                operators.add(getOperator("CP1+2.cg"));
                                operators.add(getOperator("CP1+2.gt"));

                                operators.add(getOperator("CP3.ac"));
                                operators.add(getOperator("CP3.ag"));
                                operators.add(getOperator("CP3.at"));
                                operators.add(getOperator("CP3.cg"));
                                operators.add(getOperator("CP3.gt"));
                                break;

                            default:
                                throw new IllegalArgumentException("Unknown nucleotides substitution model");
                        }
                        if (frequencyPolicy == FrequencyPolicy.ESTIMATED) {
                            if (getCodonPartitionCount() > 1 && unlinkedSubstitutionModel) {
                                operators.add(getOperator("CP1+2.frequencies"));
                                operators.add(getOperator("CP3.frequencies"));
                            } else {
                                operators.add(getOperator("frequencies"));
                            }
                        }

                    } else {
                        throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                    }
                } else { // no codon partitioning
                    switch (nucSubstitutionModel) {
                        case HKY:
                            operators.add(getOperator("kappa"));
                            break;

                        case GTR:
                            operators.add(getOperator("ac"));
                            operators.add(getOperator("ag"));
                            operators.add(getOperator("at"));
                            operators.add(getOperator("cg"));
                            operators.add(getOperator("gt"));
                            break;

                        default:
                            throw new IllegalArgumentException("Unknown nucleotides substitution model");
                    }
                    if (frequencyPolicy == FrequencyPolicy.ESTIMATED) {
                        operators.add(getOperator("frequencies"));
                    }
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
                        operators.add(getOperator("bcov.alpha"));
                        operators.add(getOperator("bcov.s"));
                        operators.add(getOperator("bcov.frequencies"));
                        operators.add(getOperator("bcov.hfrequencies"));
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
            if (getCodonPartitionCount() > 1 && unlinkedHeterogeneityModel) {
                if (codonHeteroPattern.equals("123")) {
                    operators.add(getOperator("CP1.alpha"));
                    operators.add(getOperator("CP2.alpha"));
                    operators.add(getOperator("CP3.alpha"));
                } else if (codonHeteroPattern.equals("112")) {
                    operators.add(getOperator("CP1+2.alpha"));
                    operators.add(getOperator("CP3.alpha"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                operators.add(getOperator("alpha"));
            }
        }
        // if pinv do pinv move
        if (invarHetero) {
            if (getCodonPartitionCount() > 1 && unlinkedHeterogeneityModel) {
                if (codonHeteroPattern.equals("123")) {
                    operators.add(getOperator("CP1.pInv"));
                    operators.add(getOperator("CP2.alpha"));
                    operators.add(getOperator("CP3.pInv"));
                } else if (codonHeteroPattern.equals("112")) {
                    operators.add(getOperator("CP1+2.pInv"));
                    operators.add(getOperator("CP3.pInv"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                operators.add(getOperator("pInv"));
            }
        }

        return operators;
    }


    /**
     * @param includeRelativeRates true if relative rate parameters should be added
     * @return a list of parameters that are required
     */
    List<Parameter> getParameters(boolean includeRelativeRates) {

        List<Parameter> params = new ArrayList<Parameter>();

        switch (dataType.getType()) {
            case DataType.NUCLEOTIDES:
                if (getCodonPartitionCount() > 1 && unlinkedSubstitutionModel) {
                    if (codonHeteroPattern.equals("123")) {
                        switch (nucSubstitutionModel) {
                            case HKY:
                                params.add(getParameter("CP1.kappa"));
                                params.add(getParameter("CP2.kappa"));
                                params.add(getParameter("CP3.kappa"));
                                break;
                            case GTR:
                                for (int i = 1; i <= getCodonPartitionCount(); i++) {
                                    params.add(getParameter("CP" + i + ".ac"));
                                    params.add(getParameter("CP" + i + ".ag"));
                                    params.add(getParameter("CP" + i + ".at"));
                                    params.add(getParameter("CP" + i + ".cg"));
                                    params.add(getParameter("CP" + i + ".gt"));
                                }
                                break;

                            default:
                                throw new IllegalArgumentException("Unknown nucleotides substitution model");
                        }
                        params.add(getParameter("CP1.mu"));
                        params.add(getParameter("CP2.mu"));
                        params.add(getParameter("CP3.mu"));
                    } else if (codonHeteroPattern.equals("112")) {
                        switch (nucSubstitutionModel) {
                            case HKY:
                                params.add(getParameter("CP1+2.kappa"));
                                params.add(getParameter("CP3.kappa"));
                                break;
                            case GTR:
                                params.add(getParameter("CP1+2.ac"));
                                params.add(getParameter("CP1+2.ag"));
                                params.add(getParameter("CP1+2.at"));
                                params.add(getParameter("CP1+2.cg"));
                                params.add(getParameter("CP1+2.gt"));

                                params.add(getParameter("CP3.ac"));
                                params.add(getParameter("CP3.ag"));
                                params.add(getParameter("CP3.at"));
                                params.add(getParameter("CP3.cg"));
                                params.add(getParameter("CP3.gt"));
                                break;

                            default:
                                throw new IllegalArgumentException("Unknown nucleotides substitution model");
                        }
                        params.add(getParameter("CP1+2.mu"));
                        params.add(getParameter("CP3.mu"));
                    } else {
                        throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                    }
                } else { // no codon partitioning
                    switch (nucSubstitutionModel) {
                        case HKY:
                            params.add(getParameter("kappa"));
                            break;
                        case GTR:
                            params.add(getParameter("ac"));
                            params.add(getParameter("ag"));
                            params.add(getParameter("at"));
                            params.add(getParameter("cg"));
                            params.add(getParameter("gt"));
                            break;

                        default:
                            throw new IllegalArgumentException("Unknown nucleotides substitution model");
                    }
                    if (includeRelativeRates) {
                        params.add(getParameter("mu"));
                    }
                }
                break;

            case DataType.AMINO_ACIDS:
                if (includeRelativeRates) {
                    params.add(getParameter("mu"));
                }
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
                if (includeRelativeRates) {
                    params.add(getParameter("mu"));
                }
                break;

            default:
                throw new IllegalArgumentException("Unknown data type");
        }

        // if gamma do shape move
        if (gammaHetero) {
            if (getCodonPartitionCount() > 1 && unlinkedHeterogeneityModel) {
                if (codonHeteroPattern.equals("123")) {
                    params.add(getParameter("CP1.alpha"));
                    params.add(getParameter("CP2.alpha"));
                    params.add(getParameter("CP3.alpha"));
                } else if (codonHeteroPattern.equals("112")) {
                    params.add(getParameter("CP1+2.alpha"));
                    params.add(getParameter("CP3.alpha"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                params.add(getParameter("alpha"));
            }
        }
        // if pinv do pinv move
        if (invarHetero) {
            if (getCodonPartitionCount() > 1 && unlinkedHeterogeneityModel) {
                if (codonHeteroPattern.equals("123")) {
                    params.add(getParameter("CP1.pInv"));
                    params.add(getParameter("CP2.alpha"));
                    params.add(getParameter("CP3.pInv"));
                } else if (codonHeteroPattern.equals("112")) {
                    params.add(getParameter("CP1+2.pInv"));
                    params.add(getParameter("CP3.pInv"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                params.add(getParameter("pInv"));
            }
        }

        return params;
    }

    public Parameter getParameter(String name) {

        if (name.startsWith(getName())) {
            name = name.substring(getName().length() + 1);
        }
        Parameter parameter = parameters.get(name);

        if (parameter == null) {
            throw new IllegalArgumentException("Parameter with name, " + name + ", is unknown");
        }

        parameter.setPrefix(getPrefix());

        return parameter;
    }

    Operator getOperator(String name) {
        Operator operator = operators.get(name);
        if (operator == null) throw new IllegalArgumentException("Operator with name, " + name + ", is unknown");
        operator.setPrefix(getName());
        return operator;
    }

    public int getCodonPartitionCount() {
        if (codonHeteroPattern == null || codonHeteroPattern.equals("111")) {
            return 1;
        }
        if (codonHeteroPattern.equals("123")) {
            return 3;
        }
        if (codonHeteroPattern.equals("112")) {
            return 2;
        }
        throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
    }

    public void addWeightsForPartition(DataPartition partition, int[] weights, int offset) {
        int n = partition.getSiteCount();
        int codonCount = n / 3;
        int remainder = n % 3;
        if (codonHeteroPattern == null || codonHeteroPattern.equals("111")) {
            weights[offset] += n;
            return;
        }
        if (codonHeteroPattern.equals("123")) {
            weights[offset] += codonCount + (remainder > 0 ? 1 : 0);
            weights[offset + 1] += codonCount + (remainder > 1 ? 1 : 0);
            weights[offset + 2] += codonCount;
            return;
        }
        if (codonHeteroPattern.equals("112")) {
            weights[offset] += codonCount * 2 + remainder; // positions 1 + 2
            weights[offset + 1] += codonCount; // position 3
            return;
        }
        throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
    }

    public String toString() {
        return getName();
    }

    public NucModelType getNucSubstitutionModel() {
        return nucSubstitutionModel;
    }

    public void setNucSubstitutionModel(NucModelType nucSubstitutionModel) {
        this.nucSubstitutionModel = nucSubstitutionModel;
    }

    public int getAaSubstitutionModel() {
        return aaSubstitutionModel;
    }

    public void setAaSubstitutionModel(int aaSubstitutionModel) {
        this.aaSubstitutionModel = aaSubstitutionModel;
    }

    public int getBinarySubstitutionModel() {
        return binarySubstitutionModel;
    }

    public void setBinarySubstitutionModel(int binarySubstitutionModel) {
        this.binarySubstitutionModel = binarySubstitutionModel;
    }

    public FrequencyPolicy getFrequencyPolicy() {
        return frequencyPolicy;
    }

    public void setFrequencyPolicy(FrequencyPolicy frequencyPolicy) {
        this.frequencyPolicy = frequencyPolicy;
    }

    public boolean isGammaHetero() {
        return gammaHetero;
    }

    public void setGammaHetero(boolean gammaHetero) {
        this.gammaHetero = gammaHetero;
    }

    public int getGammaCategories() {
        return gammaCategories;
    }

    public void setGammaCategories(int gammaCategories) {
        this.gammaCategories = gammaCategories;
    }

    public boolean isInvarHetero() {
        return invarHetero;
    }

    public void setInvarHetero(boolean invarHetero) {
        this.invarHetero = invarHetero;
    }

    public String getCodonHeteroPattern() {
        return codonHeteroPattern;
    }

    public void setCodonHeteroPattern(String codonHeteroPattern) {
        this.codonHeteroPattern = codonHeteroPattern;
    }

    /**
     * @return true if the rate matrix parameters are unlinked across codon positions
     */
    public boolean isUnlinkedSubstitutionModel() {
        return unlinkedSubstitutionModel;
    }

    public void setUnlinkedSubstitutionModel(boolean unlinkedSubstitutionModel) {
        this.unlinkedSubstitutionModel = unlinkedSubstitutionModel;
    }

    public boolean isUnlinkedHeterogeneityModel() {
        return unlinkedHeterogeneityModel;
    }

    public void setUnlinkedHeterogeneityModel(boolean unlinkedHeterogeneityModel) {
        this.unlinkedHeterogeneityModel = unlinkedHeterogeneityModel;
    }

    public boolean isUnlinkedFrequencyModel() {
        return unlinkedFrequencyModel;
    }

    public void setUnlinkedFrequencyModel(boolean unlinkedFrequencyModel) {
        this.unlinkedFrequencyModel = unlinkedFrequencyModel;
    }

    public DataType getDataType() {
        return dataType;
    }

    public void setDataType(DataType dataType) {
        this.dataType = dataType;
    }

    public String getPrefix() {
        String prefix = "";
        if (options.getActiveModels().size() > 1) {
            // There is more than one active partition model
            prefix += getName() + ".";
        }
        return prefix;
    }

    public String getPrefix(int codonPartitionNumber) {
        String prefix = "";
        if (options.getActiveModels().size() > 1) {
            // There is more than one active partition model
            prefix += getName() + ".";
        }
        if (getCodonPartitionCount() > 1 && codonPartitionNumber > 0) {
            if (getCodonHeteroPattern().equals("123")) {
                prefix += "CP" + codonPartitionNumber + ".";
            } else if (getCodonHeteroPattern().equals("112")) {
                if (codonPartitionNumber == 1) {
                    prefix += "CP1+2.";
                } else {
                    prefix += "CP3.";
                }
            } else {
                throw new IllegalArgumentException("unsupported codon hetero pattern");
            }

        }
        return prefix;
    }

    private final BeautiOptions options;

    private NucModelType nucSubstitutionModel = NucModelType.HKY;
    private int aaSubstitutionModel = BeautiOptions.BLOSUM_62;
    private int binarySubstitutionModel = BeautiOptions.BIN_SIMPLE;

    private FrequencyPolicy frequencyPolicy = FrequencyPolicy.ESTIMATED;
    private boolean gammaHetero = false;
    private int gammaCategories = 4;
    private boolean invarHetero = false;
    private String codonHeteroPattern = null;
    private boolean unlinkedSubstitutionModel = false;
    private boolean unlinkedHeterogeneityModel = false;

    private boolean unlinkedFrequencyModel = false;

    public DataType dataType;
    public String name;

}
