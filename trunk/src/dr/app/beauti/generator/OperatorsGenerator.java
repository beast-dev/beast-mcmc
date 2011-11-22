package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.*;
import dr.app.beauti.types.ClockType;
import dr.app.beauti.types.FixRateType;
import dr.app.beauti.types.RelativeRatesType;
import dr.app.beauti.types.TreePriorType;
import dr.app.beauti.util.XMLWriter;
import dr.evomodel.operators.BitFlipInSubstitutionModelOperator;
import dr.evomodel.substmodel.AbstractSubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.GMRFSkyrideLikelihoodParser;
import dr.evomodelxml.coalescent.VariableDemographicModelParser;
import dr.evomodelxml.coalescent.operators.GMRFSkyrideBlockUpdateOperatorParser;
import dr.evomodelxml.coalescent.operators.SampleNonActiveGibbsOperatorParser;
import dr.evomodelxml.operators.*;
import dr.evomodelxml.speciation.BirthDeathModelParser;
import dr.evomodelxml.speciation.SpeciesTreeModelParser;
import dr.evomodelxml.speciation.YuleModelParser;
import dr.evomodelxml.substmodel.GeneralSubstitutionModelParser;
import dr.inference.model.ParameterParser;
import dr.inference.operators.RateBitExchangeOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inferencexml.model.CompoundParameterParser;
import dr.inferencexml.operators.*;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class OperatorsGenerator extends Generator {

    public OperatorsGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * Write the operator schedule XML block.
     *
     * @param operators the list of operators
     * @param writer    the writer
     */
    public void writeOperatorSchedule(List<Operator> operators, XMLWriter writer) {
        Attribute[] operatorAttributes;
//		switch (options.coolingSchedule) {
//			case SimpleOperatorSchedule.LOG_SCHEDULE:
//        if (options.nodeHeightPrior == TreePriorType.GMRF_SKYRIDE) {
        // TODO: multi-prior, currently simplify to share same prior case
        if (options.isShareSameTreePrior() && options.getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.GMRF_SKYRIDE) {
            operatorAttributes = new Attribute[2];
            operatorAttributes[1] = new Attribute.Default<String>(SimpleOperatorScheduleParser.OPTIMIZATION_SCHEDULE, SimpleOperatorSchedule.LOG_STRING);
        } else {
//				break;
//			default:
            operatorAttributes = new Attribute[1];
        }
        operatorAttributes[0] = new Attribute.Default<String>(XMLParser.ID, "operators");

        writer.writeComment("Define operators");
        writer.writeOpenTag(
                SimpleOperatorScheduleParser.OPERATOR_SCHEDULE,
                operatorAttributes
//				new Attribute[]{new Attribute.Default<String>(XMLParser.ID, "operators")}
        );

        for (Operator operator : operators) {
            if (operator.weight > 0. && operator.inUse) {
            	setModelPrefix(operator.getPrefix());

            	writeOperator(operator, writer);
            }
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_OPERATORS, writer); // Added for special operators

        writer.writeCloseTag(SimpleOperatorScheduleParser.OPERATOR_SCHEDULE);
    }

    private void writeOperator(Operator operator, XMLWriter writer) {

        switch (operator.operatorType) {

            case SCALE:
                writeScaleOperator(operator, writer);
                break;
            case RANDOM_WALK:
                writeRandomWalkOperator(operator, writer);
                break;
            case RANDOM_WALK_ABSORBING:
                writeRandomWalkOperator(operator, false, writer);
                break;
            case RANDOM_WALK_INT:
                writeRandomWalkIntegerOperator(operator, writer);
                break;
            case RANDOM_WALK_REFLECTING:
                writeRandomWalkOperator(operator, true, writer);
                break;
            case INTEGER_RANDOM_WALK:
                writeIntegerRandomWalkOperator(operator, writer);
                break;
            case UP_DOWN:
                writeUpDownOperator(UpDownOperatorParser.UP_DOWN_OPERATOR, operator, writer);
                break;
            case MICROSAT_UP_DOWN:
                writeUpDownOperator(MicrosatUpDownOperatorParser.MICROSAT_UP_DOWN_OPERATOR, operator, writer);
                break;
            case UP_DOWN_ALL_RATES_HEIGHTS:
            	writeUpDownOperatorAllRatesTrees(operator, writer);
                break;
            case SCALE_ALL:
                writeScaleAllOperator(operator, writer);
                break;
            case SCALE_INDEPENDENTLY:
                writeScaleOperator(operator, writer, true);
                break;
            case CENTERED_SCALE:
                writeCenteredOperator(operator, writer);
                break;
            case DELTA_EXCHANGE:
                writeDeltaOperator(operator, writer);
                break;
            case INTEGER_DELTA_EXCHANGE:
                writeIntegerDeltaOperator(operator, writer);
                break;
            case SWAP:
                writeSwapOperator(operator, writer);
                break;
            case BITFLIP:
                writeBitFlipOperator(operator, writer);
                break;
            case BITFIP_IN_SUBST:
                writeBitFlipInSubstOperator(operator, writer);
                break;
            case RATE_BIT_EXCHANGE:
                writeRateBitExchangeOperator(operator, writer);
                break;
            case TREE_BIT_MOVE:
                writeTreeBitMoveOperator(operator, writer);
                break;
            case UNIFORM:
                writeUniformOperator(operator, writer);
                break;
            case INTEGER_UNIFORM:
                writeIntegerUniformOperator(operator, writer);
                break;
            case SUBTREE_SLIDE:
                writeSubtreeSlideOperator(operator, writer);
                break;
            case NARROW_EXCHANGE:
                writeNarrowExchangeOperator(operator, writer);
                break;
            case WIDE_EXCHANGE:
                writeWideExchangeOperator(operator, writer);
                break;
            case WILSON_BALDING:
                writeWilsonBaldingOperator(operator, writer);
                break;
            case SAMPLE_NONACTIVE:
                writeSampleNonActiveOperator(operator, writer);
                break;
            case SCALE_WITH_INDICATORS:
                writeScaleWithIndicatorsOperator(operator, writer);
                break;
            case GMRF_GIBBS_OPERATOR:
                writeGMRFGibbsOperator(operator, writer);
                break;
            case NODE_REHIGHT:
            	writeSpeciesTreeOperator(operator, writer);
            	break;
            default:
                throw new IllegalArgumentException("Unknown operator type");
        }
    }

    private void writeParameter1Ref(XMLWriter writer, Operator operator) {
        writer.writeIDref(ParameterParser.PARAMETER, operator.parameter1.getName());
    }

    private void writeParameter2Ref(XMLWriter writer, Operator operator) {
        writer.writeIDref(ParameterParser.PARAMETER, operator.parameter2.getName());
    }

    private void writeOperatorRef(XMLWriter writer, Operator operator) {
        writer.writeIDref(ParameterParser.PARAMETER, operator.getName());
    }

    private void writeScaleOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(
                ScaleOperatorParser.SCALE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperatorParser.SCALE_FACTOR, operator.tuning),
                        getWeightAttribute(operator.weight)
                });
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(ScaleOperatorParser.SCALE_OPERATOR);
    }

    private void writeScaleOperator(Operator operator, XMLWriter writer, boolean indepedently) {
        writer.writeOpenTag(
                ScaleOperatorParser.SCALE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperatorParser.SCALE_FACTOR, operator.tuning),
                        getWeightAttribute(operator.weight),
                        new Attribute.Default<String>(ScaleOperatorParser.SCALE_ALL_IND, indepedently ? "true" : "false")
                });
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(ScaleOperatorParser.SCALE_OPERATOR);
    }

    private void writeRandomWalkOperator(Operator operator, XMLWriter writer) {
        final String name = RandomWalkOperatorParser.RANDOM_WALK_OPERATOR;
        writer.writeOpenTag(
                name,
                new Attribute[]{
                        new Attribute.Default<Double>("windowSize", operator.tuning),
                        getWeightAttribute(operator.weight)
                });
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(name);
    }

    private void writeRandomWalkOperator(Operator operator, boolean reflecting, XMLWriter writer) {
        final String name = RandomWalkOperatorParser.RANDOM_WALK_OPERATOR;
        writer.writeOpenTag(
                name,
                new Attribute[]{
                        new Attribute.Default<Double>("windowSize", operator.tuning),
                        getWeightAttribute(operator.weight),
                        new Attribute.Default<String>("boundaryCondition",
                                (reflecting ? "reflecting" : "absorbing"))
                });
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(name);
    }

    private void writeRandomWalkIntegerOperator(Operator operator, XMLWriter writer) {
        final String name = RandomWalkIntegerOperatorParser.RANDOM_WALK_INTEGER_OPERATOR;
        writer.writeOpenTag(
                name,
                new Attribute[]{
                        new Attribute.Default<Double>("windowSize", operator.tuning),
                        getWeightAttribute(operator.weight)
                });
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(name);
    }

    private void writeIntegerRandomWalkOperator(Operator operator, XMLWriter writer) {

        int windowSize = (int) Math.round(operator.tuning);
        if (windowSize < 1) windowSize = 1;
        final String name = RandomWalkIntegerOperatorParser.RANDOM_WALK_INTEGER_OPERATOR;
        writer.writeOpenTag(
                name,
                new Attribute[]{
                        new Attribute.Default<Integer>("windowSize", windowSize),
                        getWeightAttribute(operator.weight)
                });
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(name);
    }

    private void writeScaleAllOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(
                ScaleOperatorParser.SCALE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperatorParser.SCALE_FACTOR, operator.tuning),
                        new Attribute.Default<String>(ScaleOperatorParser.SCALE_ALL, "true"),
                        getWeightAttribute(operator.weight)
                });

        if (operator.parameter2 == null) {
            writeParameter1Ref(writer, operator);
        } else {
            writer.writeOpenTag(CompoundParameterParser.COMPOUND_PARAMETER);
            writeParameter1Ref(writer, operator);
//            writer.writeIDref(ParameterParser.PARAMETER, operator.parameter2.getName());
            writeParameter2Ref(writer, operator);
            writer.writeCloseTag(CompoundParameterParser.COMPOUND_PARAMETER);
        }

        writer.writeCloseTag(ScaleOperatorParser.SCALE_OPERATOR);
    }

    private void writeCenteredOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(CenteredScaleOperatorParser.CENTERED_SCALE,
                new Attribute[]{
                        new Attribute.Default<Double>(CenteredScaleOperatorParser.SCALE_FACTOR, operator.tuning),
                        getWeightAttribute(operator.weight)
                }
        );
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(CenteredScaleOperatorParser.CENTERED_SCALE);
    }

    private void writeDeltaOperator(Operator operator, XMLWriter writer) {

        if (operator.getBaseName().startsWith(RelativeRatesType.MU_RELATIVE_RATES.toString())) {

            int[] parameterWeights = ((PartitionSubstitutionModel) operator.parameter1.getOptions()).getPartitionCodonWeights();

            if (parameterWeights != null && parameterWeights.length > 1) {
                String pw = "" + parameterWeights[0];
                for (int i = 1; i < parameterWeights.length; i++) {
                    pw += " " + parameterWeights[i];
                }
                writer.writeOpenTag(DeltaExchangeOperatorParser.DELTA_EXCHANGE,
                        new Attribute[]{
                                new Attribute.Default<Double>(DeltaExchangeOperatorParser.DELTA, operator.tuning),
                                new Attribute.Default<String>(DeltaExchangeOperatorParser.PARAMETER_WEIGHTS, pw),
                                getWeightAttribute(operator.weight)
                        }
                );
            }

        } else if (operator.getBaseName().startsWith(RelativeRatesType.CLOCK_RELATIVE_RATES.toString())) {

        	int[] parameterWeights = options.clockModelOptions.getPartitionClockWeights(operator.getClockModelGroup());

            if (parameterWeights != null && parameterWeights.length > 1) {
                String pw = "" + parameterWeights[0];
                for (int i = 1; i < parameterWeights.length; i++) {
                    pw += " " + parameterWeights[i];
                }
                writer.writeOpenTag(DeltaExchangeOperatorParser.DELTA_EXCHANGE,
                        new Attribute[]{
                                new Attribute.Default<Double>(DeltaExchangeOperatorParser.DELTA, operator.tuning),
                                new Attribute.Default<String>(DeltaExchangeOperatorParser.PARAMETER_WEIGHTS, pw),
                                getWeightAttribute(operator.weight)
                        }
                );
            }

        } else {
            writer.writeOpenTag(DeltaExchangeOperatorParser.DELTA_EXCHANGE,
                    new Attribute[]{
                            new Attribute.Default<Double>(DeltaExchangeOperatorParser.DELTA, operator.tuning),
                            getWeightAttribute(operator.weight)
                    }
            );
        }

        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(DeltaExchangeOperatorParser.DELTA_EXCHANGE);
    }

    private void writeIntegerDeltaOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(DeltaExchangeOperatorParser.DELTA_EXCHANGE,
                new Attribute[]{
                        new Attribute.Default<String>(DeltaExchangeOperatorParser.DELTA, Integer.toString((int) operator.tuning)),
                        new Attribute.Default<String>("integer", "true"),
                        getWeightAttribute(operator.weight),
                        new Attribute.Default<String>("autoOptimize", "false")
                }
        );
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(DeltaExchangeOperatorParser.DELTA_EXCHANGE);
    }

    private void writeSwapOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(SwapOperatorParser.SWAP_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<String>("size", Integer.toString((int) operator.tuning)),
                        getWeightAttribute(operator.weight),
                        new Attribute.Default<String>("autoOptimize", "false")
                }
        );
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(SwapOperatorParser.SWAP_OPERATOR);
    }

    private void writeBitFlipOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(BitFlipOperatorParser.BIT_FLIP_OPERATOR,
                getWeightAttribute(operator.weight));
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(BitFlipOperatorParser.BIT_FLIP_OPERATOR);
    }

    private void writeBitFlipInSubstOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(BitFlipInSubstitutionModelOperator.BIT_FLIP_OPERATOR, new Attribute[]{
                new Attribute.Default<Double>(ScaleOperatorParser.SCALE_FACTOR, operator.tuning),
                getWeightAttribute(operator.weight)});
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        PartitionSubstitutionModel model = (PartitionSubstitutionModel) operator.getOptions();
        writer.writeIDref(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, model.getPrefix() + AbstractSubstitutionModel.MODEL);
        // <svsGeneralSubstitutionModel idref="originModel"/>
        writer.writeCloseTag(BitFlipInSubstitutionModelOperator.BIT_FLIP_OPERATOR);
    }

    private void writeRateBitExchangeOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(RateBitExchangeOperator.OPERATOR_NAME,
                getWeightAttribute(operator.weight));

        writer.writeOpenTag(RateBitExchangeOperator.BITS);
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(RateBitExchangeOperator.BITS);

        writer.writeOpenTag(RateBitExchangeOperator.RATES);
        writeParameter2Ref(writer, operator);
        writer.writeCloseTag(RateBitExchangeOperator.RATES);

        writer.writeCloseTag(RateBitExchangeOperator.OPERATOR_NAME);
    }

    private void writeTreeBitMoveOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(TreeBitMoveOperatorParser.BIT_MOVE_OPERATOR,
                        getWeightAttribute(operator.weight));
        writer.writeIDref(TreeModel.TREE_MODEL,  modelPrefix + TreeModel.TREE_MODEL);
        writer.writeCloseTag(TreeBitMoveOperatorParser.BIT_MOVE_OPERATOR);
    }

    private void writeUniformOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag("uniformOperator",
                getWeightAttribute(operator.weight));
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag("uniformOperator");
    }

    private void writeIntegerUniformOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(UniformIntegerOperatorParser.UNIFORM_INTEGER_OPERATOR,
                getWeightAttribute(operator.weight));
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(UniformIntegerOperatorParser.UNIFORM_INTEGER_OPERATOR);
    }

    private void writeNarrowExchangeOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(ExchangeOperatorParser.NARROW_EXCHANGE,
                getWeightAttribute(operator.weight));
        writer.writeIDref(TreeModel.TREE_MODEL,  modelPrefix + TreeModel.TREE_MODEL);
        writer.writeCloseTag(ExchangeOperatorParser.NARROW_EXCHANGE);
    }

    private void writeWideExchangeOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(ExchangeOperatorParser.WIDE_EXCHANGE,
                getWeightAttribute(operator.weight));
        writer.writeIDref(TreeModel.TREE_MODEL,  modelPrefix + TreeModel.TREE_MODEL);
        writer.writeCloseTag(ExchangeOperatorParser.WIDE_EXCHANGE);
    }

    private void writeWilsonBaldingOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(WilsonBaldingParser.WILSON_BALDING,
                getWeightAttribute(operator.weight));
        writer.writeIDref(TreeModel.TREE_MODEL,  modelPrefix + TreeModel.TREE_MODEL);
        // not supported anymore. probably never worked. (todo) get it out of GUI too
//        if (options.nodeHeightPrior == TreePriorType.CONSTANT) {
//            treePriorGenerator.writeNodeHeightPriorModelRef(writer);
//        }
        writer.writeCloseTag(WilsonBaldingParser.WILSON_BALDING);
    }

    private void writeSampleNonActiveOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(SampleNonActiveGibbsOperatorParser.SAMPLE_NONACTIVE_GIBBS_OPERATOR,
                getWeightAttribute(operator.weight));

        writer.writeOpenTag(SampleNonActiveGibbsOperatorParser.DISTRIBUTION);
        writeOperatorRef(writer, operator);
        writer.writeCloseTag(SampleNonActiveGibbsOperatorParser.DISTRIBUTION);

        writer.writeOpenTag(SampleNonActiveGibbsOperatorParser.DATA_PARAMETER);
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(SampleNonActiveGibbsOperatorParser.DATA_PARAMETER);

        writer.writeOpenTag(SampleNonActiveGibbsOperatorParser.INDICATOR_PARAMETER);
        writeParameter2Ref(writer, operator);
        writer.writeCloseTag(SampleNonActiveGibbsOperatorParser.INDICATOR_PARAMETER);

        writer.writeCloseTag(SampleNonActiveGibbsOperatorParser.SAMPLE_NONACTIVE_GIBBS_OPERATOR);
    }

    private void writeGMRFGibbsOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(
                GMRFSkyrideBlockUpdateOperatorParser.BLOCK_UPDATE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(GMRFSkyrideBlockUpdateOperatorParser.SCALE_FACTOR, operator.tuning),
                        getWeightAttribute(operator.weight)
                }
        );
        writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD,  modelPrefix + "skyride");
        writer.writeCloseTag(GMRFSkyrideBlockUpdateOperatorParser.BLOCK_UPDATE_OPERATOR);
    }

    private void writeScaleWithIndicatorsOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(
                ScaleOperatorParser.SCALE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperatorParser.SCALE_FACTOR, operator.tuning),
                        getWeightAttribute(operator.weight)
                });
        writeParameter1Ref(writer, operator);
        writer.writeOpenTag(ScaleOperatorParser.INDICATORS, new Attribute.Default<String>(ScaleOperatorParser.PICKONEPROB, "1.0"));
        writeParameter2Ref(writer, operator);
        writer.writeCloseTag(ScaleOperatorParser.INDICATORS);
        writer.writeCloseTag(ScaleOperatorParser.SCALE_OPERATOR);
    }

    private void writeSubtreeSlideOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(SubtreeSlideOperatorParser.SUBTREE_SLIDE,
                new Attribute[]{
                        new Attribute.Default<Double>("size", operator.tuning),
                        new Attribute.Default<String>("gaussian", "true"),
                        getWeightAttribute(operator.weight)
                }
        );
        writer.writeIDref(TreeModel.TREE_MODEL,  modelPrefix + TreeModel.TREE_MODEL);
        writer.writeCloseTag(SubtreeSlideOperatorParser.SUBTREE_SLIDE);
    }

    private void writeSpeciesTreeOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(TreeNodeSlideParser.TREE_NODE_REHEIGHT,
                new Attribute[]{ getWeightAttribute(operator.weight) }
        );
        writer.writeIDref(TraitData.TRAIT_SPECIES,  TraitData.TRAIT_SPECIES);
        writer.writeIDref(SpeciesTreeModelParser.SPECIES_TREE,  Generator.SP_TREE);
        writer.writeCloseTag(TreeNodeSlideParser.TREE_NODE_REHEIGHT);
    }

    private void writeUpDownOperator(String opTag, Operator operator, XMLWriter writer) {
        writer.writeOpenTag(opTag,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperatorParser.SCALE_FACTOR, operator.tuning),
                        getWeightAttribute(operator.weight)
                }
        );

        writer.writeOpenTag(UpDownOperatorParser.UP);
        // for isEstimatedRate() = false, write nothing on up part of upDownOp
        if (!operator.parameter1.isFixed && operator.getClockModelGroup().getRateTypeOption() != FixRateType.FIX_MEAN) {
        	writeParameter1Ref(writer, operator);
        }
        writer.writeCloseTag(UpDownOperatorParser.UP);

        writer.writeOpenTag(UpDownOperatorParser.DOWN);
        if (operator.tag == null) {
//	        writer.writeIDref(ParameterParser.PARAMETER,  operator.parameter2.getName());
            writeParameter2Ref(writer, operator);
        } else {
        	writer.writeIDref(operator.tag,  operator.idref);
        }
        writer.writeCloseTag(UpDownOperatorParser.DOWN);

        writer.writeCloseTag(opTag);
    }

    private void writeUpDownOperatorAllRatesTrees(Operator operator, XMLWriter writer) {
    	writer.writeOpenTag(UpDownOperatorParser.UP_DOWN_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperatorParser.SCALE_FACTOR, operator.tuning),
                        getWeightAttribute(operator.weight)
                }
        );

        writer.writeOpenTag(UpDownOperatorParser.UP);

        for (PartitionClockModel model : options.getPartitionClockModels()) {
			if (model.isEstimatedRate()) {
				switch (model.getClockType()) {
	            case STRICT_CLOCK:
	            case RANDOM_LOCAL_CLOCK:
	            	writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "clock.rate");
	                break;

	            case UNCORRELATED:
                    switch (model.getClockDistributionType()) {
                        case LOGNORMAL:
                            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCLD_MEAN);
                            break;
                        case GAMMA:
                            throw new UnsupportedOperationException("Uncorrelated gamma relaxed clock model not implemented yet");
//                            break;
                        case CAUCHY:
                            throw new UnsupportedOperationException("Uncorrelated Cauchy relaxed clock model not implemented yet");
//                            break;
                        case EXPONENTIAL:
                            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCED_MEAN);
                            break;
                    }
	            	break;

	            case AUTOCORRELATED:
                    throw new UnsupportedOperationException("Autocorrelated relaxed clock model not implemented yet");
//	                break;

	            default:
	                throw new IllegalArgumentException("Unknown clock model");
				}
			}
        }
        if (options.useStarBEAST) {
	        if (options.getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_BIRTH_DEATH) {
	        	writer.writeIDref(ParameterParser.PARAMETER, TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME);
	        } else if (options.getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_YULE) {
	        	writer.writeIDref(ParameterParser.PARAMETER, TraitData.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE);
	        }
        }// nothing for EBSP

        writer.writeCloseTag(UpDownOperatorParser.UP);

        writer.writeOpenTag(UpDownOperatorParser.DOWN);

        if (options.useStarBEAST) {
	        writer.writeIDref(SpeciesTreeModelParser.SPECIES_TREE, SP_TREE); // <speciesTree idref="sptree" /> has to be the 1st always
	        writer.writeIDref(ParameterParser.PARAMETER, TraitData.TRAIT_SPECIES + "." + options.starBEASTOptions.POP_MEAN);
	        writer.writeIDref(ParameterParser.PARAMETER, SpeciesTreeModelParser.SPECIES_TREE + "." + SPLIT_POPS);
        } else if (options.isEBSPSharingSamePrior()) {
        	writer.writeIDref(ParameterParser.PARAMETER, VariableDemographicModelParser.demoElementName + ".populationMean");
	        writer.writeIDref(ParameterParser.PARAMETER, VariableDemographicModelParser.demoElementName + ".popSize");
        }

        for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
        	writer.writeIDref(ParameterParser.PARAMETER, tree.getPrefix() + "treeModel.allInternalNodeHeights");
        }

        writer.writeCloseTag(UpDownOperatorParser.DOWN);

        writer.writeCloseTag(UpDownOperatorParser.UP_DOWN_OPERATOR);
    }

    private Attribute getWeightAttribute(double weight) {
        if (weight == (int)weight) {
            return new Attribute.Default<Integer>("weight", (int)weight);
        } else {
            return new Attribute.Default<Double>("weight", weight);
        }
    }
}
