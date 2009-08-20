/*
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.options;

import dr.evomodel.substmodel.AminoAcidModelType;
import dr.app.beauti.enumTypes.FrequencyPolicyType;
import dr.evomodel.substmodel.NucModelType;
import dr.app.beauti.enumTypes.OperatorType;
import dr.app.beauti.enumTypes.PriorScaleType;
import dr.evolution.datatype.DataType;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class PartitionSubstitutionModel extends ModelOptions {

    // Instance variables

    public static final String[] GTR_RATE_NAMES = {"ac", "ag", "at", "cg", "gt"};
    private static final String[] GTR_TRANSITIONS = {"A-C", "A-G", "A-T", "C-G", "G-T"};

    private final BeautiOptions options;
    private DataType dataType;
    private String name;

    private List<PartitionData> allPartitionData = new ArrayList<PartitionData>();

    private NucModelType nucSubstitutionModel = NucModelType.HKY;
    private AminoAcidModelType aaSubstitutionModel = AminoAcidModelType.BLOSUM_62;
    private int binarySubstitutionModel = BeautiOptions.BIN_SIMPLE;

    private FrequencyPolicyType frequencyPolicy = FrequencyPolicyType.ESTIMATED;
    private boolean gammaHetero = false;
    private int gammaCategories = 4;
    private boolean invarHetero = false;
    private String codonHeteroPattern = null;
    private boolean unlinkedSubstitutionModel = false;
    private boolean unlinkedHeterogeneityModel = false;

    private boolean unlinkedFrequencyModel = false;

    private boolean dolloModel = false;

    public PartitionSubstitutionModel(BeautiOptions options, PartitionData partition) {
        this(options, partition.getName(), partition.getAlignment().getDataType());

        allPartitionData.clear();
        addPartitionData(partition);
    }

    /**
     * A copy constructor
     *
     * @param options the beauti options
     * @param name    the name of the new model
     * @param source  the source model
     */
    public PartitionSubstitutionModel(BeautiOptions options, String name, PartitionSubstitutionModel source) {
        this(options, name, source.dataType);
        
        this.allPartitionData.clear();
        for (PartitionData partition: source.allPartitionData) {
        	this.allPartitionData.add(partition);			
		}        

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

    public PartitionSubstitutionModel(BeautiOptions options, String name, DataType dataType) {

        this.options = options;
        this.name = name;
        this.dataType = dataType;

        initSubstModelParaAndOpers();
    }


    // only init in PartitionSubstitutionModel
    private void initSubstModelParaAndOpers() {
        double substWeights = 1.0;

        //Substitution model parameters
        createParameter("frequencies", "base frequencies", PriorScaleType.UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("CP1.frequencies", "base frequencies for codon position 1", PriorScaleType.UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("CP2.frequencies", "base frequencies for codon position 2", PriorScaleType.UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("CP1+2.frequencies", "base frequencies for codon positions 1 & 2", PriorScaleType.UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("CP3.frequencies", "base frequencies for codon position 3", PriorScaleType.UNITY_SCALE, 0.25, 0.0, 1.0);

        createScaleParameter("kappa", "HKY transition-transversion parameter", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP1.kappa", "HKY transition-transversion parameter for codon position 1", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP2.kappa", "HKY transition-transversion parameter for codon position 2", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP1+2.kappa", "HKY transition-transversion parameter for codon positions 1 & 2", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP3.kappa", "HKY transition-transversion parameter for codon position 3", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

//        createParameter("frequencies", "GTR base frequencies", UNITY_SCALE, 0.25, 0.0, 1.0);
//        createParameter("CP1.frequencies", "GTR base frequencies for codon position 1", UNITY_SCALE, 0.25, 0.0, 1.0);
//        createParameter("CP2.frequencies", "GTR base frequencies for codon position 2", UNITY_SCALE, 0.25, 0.0, 1.0);
//        createParameter("CP1+2.frequencies", "GTR base frequencies for codon positions 1 & 2", UNITY_SCALE, 0.25, 0.0, 1.0);
//        createParameter("CP3.frequencies", "GTR base frequencies for codon position 3", UNITY_SCALE, 0.25, 0.0, 1.0);

        // create the relative rate parameters for the GTR rate matrix
        for (int j = 0; j < 5; j++) {
            createScaleParameter(GTR_RATE_NAMES[j], "GTR " + GTR_TRANSITIONS[j] + " substitution parameter",
            		PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
            for (int i = 1; i <= 3; i++) {

                createScaleParameter(
                        "CP" + i + "." + GTR_RATE_NAMES[j],
                        "GTR " + GTR_TRANSITIONS[j] + " substitution parameter for codon position " + i,
                        PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
            }
            createScaleParameter("CP1+2." + GTR_RATE_NAMES[j],
                    "GTR " + GTR_TRANSITIONS[j] + " substitution parameter for codon positions 1 & 2",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        }

//        createParameter("frequencies", "Binary Simple frequencies", UNITY_SCALE, 0.5, 0.0, 1.0);
//
//        createParameter("frequencies", "Binary Covarion frequencies of the visible states", UNITY_SCALE, 0.5, 0.0, 1.0);
        createParameter("hfrequencies", "Binary Covarion frequencies of the hidden rates", PriorScaleType.UNITY_SCALE, 0.5, 0.0, 1.0);
        createParameter("bcov.alpha", "Binary Covarion rate of evolution in slow mode", PriorScaleType.UNITY_SCALE, 0.5, 0.0, 1.0);
        createParameter("bcov.s", "Binary Covarion rate of flipping between slow and fast modes", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, 100.0);

        createParameter("alpha", "gamma shape parameter", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, 1000.0);
        createParameter("CP1.alpha", "gamma shape parameter for codon position 1", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, 1000.0);
        createParameter("CP2.alpha", "gamma shape parameter for codon position 2", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, 1000.0);
        createParameter("CP1+2.alpha", "gamma shape parameter for codon positions 1 & 2", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, 1000.0);
        createParameter("CP3.alpha", "gamma shape parameter for codon position 3", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, 1000.0);

        createParameter("pInv", "proportion of invariant sites parameter", PriorScaleType.NONE, 0.5, 0.0, 1.0);
        createParameter("CP1.pInv", "proportion of invariant sites parameter for codon position 1", PriorScaleType.NONE, 0.5, 0.0, 1.0);
        createParameter("CP2.pInv", "proportion of invariant sites parameter for codon position 2", PriorScaleType.NONE, 0.5, 0.0, 1.0);
        createParameter("CP1+2.pInv", "proportion of invariant sites parameter for codon positions 1 & 2", PriorScaleType.NONE, 0.5, 0.0, 1.0);
        createParameter("CP3.pInv", "proportion of invariant sites parameter for codon position 3", PriorScaleType.NONE, 0.5, 0.0, 1.0);

        createParameter("mu", "relative rate parameter", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("CP1.mu", "relative rate parameter for codon position 1", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("CP2.mu", "relative rate parameter for codon position 2", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("CP1+2.mu", "relative rate parameter for codon positions 1 & 2", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("CP3.mu", "relative rate parameter for codon position 3", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

        createScaleOperator("kappa", demoTuning, substWeights);
        createScaleOperator("CP1.kappa", demoTuning, substWeights);
        createScaleOperator("CP2.kappa", demoTuning, substWeights);
        createScaleOperator("CP1+2.kappa", demoTuning, substWeights);
        createScaleOperator("CP3.kappa", demoTuning, substWeights);

        createOperator("frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("CP1.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("CP2.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("CP1+2.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("CP3.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);

        for (String rateName : GTR_RATE_NAMES) {
            createScaleOperator(rateName, demoTuning, substWeights);
            for (int j = 1; j <= 3; j++) {
                createScaleOperator("CP" + j + "." + rateName, demoTuning, substWeights);
            }
            createScaleOperator("CP1+2." + rateName, demoTuning, substWeights);
        }

        createScaleOperator("alpha", demoTuning, substWeights);
        for (int i = 1; i <= 3; i++) {
            createScaleOperator("CP" + i + ".alpha", demoTuning, substWeights);
        }
        createScaleOperator("CP1+2.alpha", demoTuning, substWeights);

        createScaleOperator("pInv", demoTuning, substWeights);
        for (int i = 1; i <= 3; i++) {
            createScaleOperator("CP" + i + ".pInv", demoTuning, substWeights);
        }
        createScaleOperator("CP1+2.pInv", demoTuning, substWeights);

        createScaleOperator("bcov.alpha", demoTuning, substWeights);
        createScaleOperator("bcov.s", demoTuning, substWeights);
        //createOperator("frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("hfrequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
    }

    ////////////////////////////////////////////////////////////////

    /**
     * @param includeRelativeRates true if relative rate parameters should be added
     * @return a list of parameters that are required
     */
    public void selectParameters(List<Parameter> params) {
    	boolean includeRelativeRates = getCodonPartitionCount() > 1;//TODO check 

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
                                    for (String rateName : GTR_RATE_NAMES) {
                                        params.add(getParameter("CP" + i + "." + rateName));
                                    }
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
                                for (String rateName : GTR_RATE_NAMES) {
                                    params.add(getParameter("CP1+2." + rateName));
                                }
                                for (String rateName : GTR_RATE_NAMES) {
                                    params.add(getParameter("CP3." + rateName));
                                }
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
                            for (String rateName : GTR_RATE_NAMES) {
                                params.add(getParameter(rateName));
                            }
                            break;

                        default:
                            throw new IllegalArgumentException("Unknown nucleotides substitution model");
                    }
                    
                    if (includeRelativeRates) {//TODO fix mean
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
                    params.add(getParameter("CP2.pInv"));
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
        if (frequencyPolicy == FrequencyPolicyType.ESTIMATED) {
            if (getCodonPartitionCount() > 1 && unlinkedHeterogeneityModel) {
                if (codonHeteroPattern.equals("123")) {
                    params.add(getParameter("CP1.frequencies"));
                    params.add(getParameter("CP2.frequencies"));
                    params.add(getParameter("CP3.frequencies"));
                } else if (codonHeteroPattern.equals("112")) {
                    params.add(getParameter("CP1+2.frequencies"));
                    params.add(getParameter("CP3.frequencies"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                params.add(getParameter("frequencies"));
            }

        }        
    }

    public void selectOperators(List<Operator> ops) {
        
        switch (dataType.getType()) {
            case DataType.NUCLEOTIDES:

                if (getCodonPartitionCount() > 1 && unlinkedSubstitutionModel) {
                    if (codonHeteroPattern.equals("123")) {
                        switch (nucSubstitutionModel) {
                            case HKY:
                                ops.add(getOperator("CP1.kappa"));
                                ops.add(getOperator("CP2.kappa"));
                                ops.add(getOperator("CP3.kappa"));
                                break;

                            case GTR:
                                for (int i = 1; i <= 3; i++) {
                                    for (String rateName : GTR_RATE_NAMES) {
                                        ops.add(getOperator("CP" + i + "." + rateName));
                                    }
                                }
                                break;

                            default:
                                throw new IllegalArgumentException("Unknown nucleotides substitution model");
                        }

                        if (frequencyPolicy == FrequencyPolicyType.ESTIMATED) {
                            if (getCodonPartitionCount() > 1 && unlinkedSubstitutionModel) {
                                ops.add(getOperator("CP1.frequencies"));
                                ops.add(getOperator("CP2.frequencies"));
                                ops.add(getOperator("CP3.frequencies"));
                            } else {
                                ops.add(getOperator("frequencies"));
                            }
                        }
                    } else if (codonHeteroPattern.equals("112")) {
                        switch (nucSubstitutionModel) {
                            case HKY:
                                ops.add(getOperator("CP1+2.kappa"));
                                ops.add(getOperator("CP3.kappa"));
                                break;

                            case GTR:
                                for (String rateName : GTR_RATE_NAMES) {
                                    ops.add(getOperator("CP1+2." + rateName));
                                }
                                for (String rateName : GTR_RATE_NAMES) {
                                    ops.add(getOperator("CP3." + rateName));
                                }
                                break;

                            default:
                                throw new IllegalArgumentException("Unknown nucleotides substitution model");
                        }
                        if (frequencyPolicy == FrequencyPolicyType.ESTIMATED) {
                            if (getCodonPartitionCount() > 1 && unlinkedSubstitutionModel) {
                                ops.add(getOperator("CP1+2.frequencies"));
                                ops.add(getOperator("CP3.frequencies"));
                            } else {
                                ops.add(getOperator("frequencies"));
                            }
                        }

                    } else {
                        throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                    }
                } else { // no codon partitioning
                    switch (nucSubstitutionModel) {
                        case HKY:
                            ops.add(getOperator("kappa"));
                            break;

                        case GTR:
                            for (String rateName : GTR_RATE_NAMES) {
                                ops.add(getOperator(rateName));
                            }
                            break;

                        default:
                            throw new IllegalArgumentException("Unknown nucleotides substitution model");
                    }
                    if (frequencyPolicy == FrequencyPolicyType.ESTIMATED) {
                        ops.add(getOperator("frequencies"));
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
            if (getCodonPartitionCount() > 1 && unlinkedHeterogeneityModel) {
                if (codonHeteroPattern.equals("123")) {
                    ops.add(getOperator("CP1.alpha"));
                    ops.add(getOperator("CP2.alpha"));
                    ops.add(getOperator("CP3.alpha"));
                } else if (codonHeteroPattern.equals("112")) {
                    ops.add(getOperator("CP1+2.alpha"));
                    ops.add(getOperator("CP3.alpha"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                ops.add(getOperator("alpha"));
            }
        }
        // if pinv do pinv move
        if (invarHetero) {
            if (getCodonPartitionCount() > 1 && unlinkedHeterogeneityModel) {
                if (codonHeteroPattern.equals("123")) {
                    ops.add(getOperator("CP1.pInv"));
                    ops.add(getOperator("CP2.pInv"));
                    ops.add(getOperator("CP3.pInv"));
                } else if (codonHeteroPattern.equals("112")) {
                    ops.add(getOperator("CP1+2.pInv"));
                    ops.add(getOperator("CP3.pInv"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                ops.add(getOperator("pInv"));
            }
        }

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

    public void addWeightsForPartition(PartitionData partition, int[] weights, int offset) {
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
    
    ///////////////////////////////////////////////////////

    public Parameter getParameter(String name) {

        Parameter parameter = parameters.get(name);

        if (parameter == null) {
            throw new IllegalArgumentException("Parameter with name, " + name + ", is unknown");
        }

        parameter.setPrefix(getPrefix());

        return parameter;
    }

    public Operator getOperator(String name) {

        Operator operator = operators.get(name);

        if (operator == null) throw new IllegalArgumentException("Operator with name, " + name + ", is unknown");

        operator.setPrefix(getPrefix());

        return operator;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String toString() {
        return getName();
    }

    public List<PartitionData> getAllPartitionData() {
        return allPartitionData;
    }

    public void clearAllPartitionData() {
        this.allPartitionData.clear();
    }

    public void addPartitionData(PartitionData partition) {
        allPartitionData.add(partition);
    }

    public boolean removePartitionData(PartitionData partition) {
        return allPartitionData.remove(partition);
    }

    public NucModelType getNucSubstitutionModel() {
        return nucSubstitutionModel;
    }

    public void setNucSubstitutionModel(NucModelType nucSubstitutionModel) {
        this.nucSubstitutionModel = nucSubstitutionModel;
    }

    public AminoAcidModelType getAaSubstitutionModel() {
        return aaSubstitutionModel;
    }

    public void setAaSubstitutionModel(AminoAcidModelType aaSubstitutionModel) {
        this.aaSubstitutionModel = aaSubstitutionModel;
    }

    public int getBinarySubstitutionModel() {
        return binarySubstitutionModel;
    }

    public void setBinarySubstitutionModel(int binarySubstitutionModel) {
        this.binarySubstitutionModel = binarySubstitutionModel;
    }

    public FrequencyPolicyType getFrequencyPolicy() {
        return frequencyPolicy;
    }

    public void setFrequencyPolicy(FrequencyPolicyType frequencyPolicy) {
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

    public boolean isDolloModel() {
        return dolloModel;
    }

    public void setDolloModel(boolean dolloModel) {
        this.dolloModel = dolloModel;
    }

    public String getPrefix() {
        String prefix = "";
        if (options.getPartitionSubstitutionModels().size() > 1) { // || options.isSpeciesAnalysis()) {
            // There is more than one active partition model, or doing species analysis
            prefix += getName() + ".";
        }
        return prefix;
    }

    public String getPrefix(int codonPartitionNumber) {
        String prefix = "";
        if (options.getPartitionSubstitutionModels().size() > 1) { //|| options.isSpeciesAnalysis()) {
            // There is more than one active partition model, or doing species analysis
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
}