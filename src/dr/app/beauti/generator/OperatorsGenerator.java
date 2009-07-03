package dr.app.beauti.generator;

import dr.app.beauti.util.XMLWriter;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.Operator;
import dr.app.beauti.options.TraitGuesser;
import dr.app.beauti.options.TreePrior;
import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.evomodel.coalescent.operators.GMRFSkyrideBlockUpdateOperator;
import dr.evomodel.coalescent.operators.SampleNonActiveGibbsOperator;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.operators.TreeBitMoveOperator;
import dr.evomodel.operators.TreeNodeSlide;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.speciation.SpeciesTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.ParameterParser;
import dr.inference.operators.*;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.List;

/**
 * @author Alexei Drummond
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
//        if (options.nodeHeightPrior == TreePrior.GMRF_SKYRIDE) {
        // TODO: multi-prior, currently simplify to share same prior case
        if (options.shareSameTreePrior && options.activedSameTreePrior.getNodeHeightPrior() == TreePrior.GMRF_SKYRIDE) {
            operatorAttributes = new Attribute[2];
            operatorAttributes[1] = new Attribute.Default<String>(SimpleOperatorSchedule.OPTIMIZATION_SCHEDULE, SimpleOperatorSchedule.LOG_STRING);
        } else {
//				break;
//			default:
            operatorAttributes = new Attribute[1];
        }
        operatorAttributes[0] = new Attribute.Default<String>(XMLParser.ID, "operators");
        
        writer.writeComment("Define operators");
        writer.writeOpenTag(
                SimpleOperatorSchedule.OPERATOR_SCHEDULE,
                operatorAttributes
//				new Attribute[]{new Attribute.Default<String>(XMLParser.ID, "operators")}
        );

        for (Operator operator : operators) {
            if (operator.weight > 0. && operator.inUse)
                writeOperator(operator, writer);
            
            if (options.isSpeciesAnalysis()) { // TODO
            	this.setModelPrefix(operator.getPrefix() + ".");
            } else {
            	this.setModelPrefix("");
            }
        }

        writer.writeCloseTag(SimpleOperatorSchedule.OPERATOR_SCHEDULE);
    }

    private void writeOperator(Operator operator, XMLWriter writer) {

        switch (operator.type) {

            case SCALE:
                writeScaleOperator(operator, writer);
                break;
            case RANDOM_WALK:
                writeRandomWalkOperator(operator, writer);
            case RANDOM_WALK_ABSORBING:
                writeRandomWalkOperator(operator, false, writer);
                break;
            case RANDOM_WALK_REFLECTING:
                writeRandomWalkOperator(operator, true, writer);
                break;
            case INTEGER_RANDOM_WALK:
                writeIntegerRandomWalkOperator(operator, writer);
                break;
            case UP_DOWN:
                writeUpDownOperator(operator, writer);
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

    private Attribute getRef(String name) {
        return new Attribute.Default<String>(XMLParser.IDREF, name);
    }

    private void writeParameterRefByName(XMLWriter writer, String name) {
        writer.writeTag(ParameterParser.PARAMETER, getRef(name), true);
    }

    private void writeParameter1Ref(XMLWriter writer, Operator operator) {
        writeParameterRefByName(writer, operator.parameter1.getName());
    }

    private void writeScaleOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(
                ScaleOperator.SCALE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperator.SCALE_FACTOR, operator.tuning),
                        getWeightAttribute(operator.weight)
                });
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(ScaleOperator.SCALE_OPERATOR);
    }

    private void writeScaleOperator(Operator operator, XMLWriter writer, boolean indepedently) {
        writer.writeOpenTag(
                ScaleOperator.SCALE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperator.SCALE_FACTOR, operator.tuning),
                        getWeightAttribute(operator.weight),
                        new Attribute.Default<String>(ScaleOperator.SCALE_ALL_IND, indepedently ? "true" : "false")
                });
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(ScaleOperator.SCALE_OPERATOR);
    }

    private void writeRandomWalkOperator(Operator operator, XMLWriter writer) {
        final String name = RandomWalkOperator.PARSER.getParserName();
        writer.writeOpenTag(
                name,
                new Attribute[]{
                        new Attribute.Default<Double>("windowSize", operator.tuning),
                        getWeightAttribute(operator.weight)
                });
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(name);
    }

    private void writeRandomWalkOperator(Operator operator, boolean reflecting, XMLWriter writer) {
        final String name = RandomWalkOperator.PARSER.getParserName();
        writer.writeOpenTag(
                name,
                new Attribute[]{
                        new Attribute.Default<Double>("windowSize", operator.tuning),
                        getWeightAttribute(operator.weight),
                        new Attribute.Default<String>("boundaryCondition",
                                (reflecting ? "reflecting" : "absorbing"))
                });
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(name);
    }

    private void writeIntegerRandomWalkOperator(Operator operator, XMLWriter writer) {

        int windowSize = (int) Math.round(operator.tuning);
        if (windowSize < 1) windowSize = 1;
        final String name = RandomWalkIntegerOperator.PARSER.getParserName();
        writer.writeOpenTag(
                name,
                new Attribute[]{
                        new Attribute.Default<Integer>("windowSize", windowSize),
                        getWeightAttribute(operator.weight)
                });
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(name);
    }

    private void writeScaleAllOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(
                ScaleOperator.SCALE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperator.SCALE_FACTOR, operator.tuning),
                        new Attribute.Default<String>(ScaleOperator.SCALE_ALL, "true"),
                        getWeightAttribute(operator.weight)
                });

        if (operator.parameter2 == null) {
            writeParameter1Ref(writer, operator);
        } else {
            writer.writeOpenTag(CompoundParameter.COMPOUND_PARAMETER);
            writeParameter1Ref(writer, operator);
            writer.writeIDref(ParameterParser.PARAMETER, operator.parameter2.getName());
            writer.writeCloseTag(CompoundParameter.COMPOUND_PARAMETER);
        }

        writer.writeCloseTag(ScaleOperator.SCALE_OPERATOR);
    }

    private void writeUpDownOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(UpDownOperator.UP_DOWN_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperator.SCALE_FACTOR, operator.tuning),
                        getWeightAttribute(operator.weight)
                }
        );

        writer.writeOpenTag(UpDownOperator.UP);
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(UpDownOperator.UP);

        writer.writeOpenTag(UpDownOperator.DOWN);
        writer.writeIDref(ParameterParser.PARAMETER,  operator.parameter2.getName());
        writer.writeCloseTag(UpDownOperator.DOWN);

        writer.writeCloseTag(UpDownOperator.UP_DOWN_OPERATOR);
    }

    private void writeCenteredOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(CenteredScaleOperator.CENTERED_SCALE,
                new Attribute[]{
                        new Attribute.Default<Double>(CenteredScaleOperator.SCALE_FACTOR, operator.tuning),
                        getWeightAttribute(operator.weight)
                }
        );
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(CenteredScaleOperator.CENTERED_SCALE);
    }

    private void writeDeltaOperator(Operator operator, XMLWriter writer) {

        if (operator.getName().equals("Relative rates")) {

            int[] parameterWeights = options.getPartitionWeights();

            if (parameterWeights != null && parameterWeights.length > 1) {
                String pw = "" + parameterWeights[0];
                for (int i = 1; i < parameterWeights.length; i++) {
                    pw += " " + parameterWeights[i];
                }
                writer.writeOpenTag(DeltaExchangeOperator.DELTA_EXCHANGE,
                        new Attribute[]{
                                new Attribute.Default<Double>(DeltaExchangeOperator.DELTA, operator.tuning),
                                new Attribute.Default<String>(DeltaExchangeOperator.PARAMETER_WEIGHTS, pw),
                                getWeightAttribute(operator.weight)
                        }
                );
            }
        } else {
            writer.writeOpenTag(DeltaExchangeOperator.DELTA_EXCHANGE,
                    new Attribute[]{
                            new Attribute.Default<Double>(DeltaExchangeOperator.DELTA, operator.tuning),
                            getWeightAttribute(operator.weight)
                    }
            );
        }

        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(DeltaExchangeOperator.DELTA_EXCHANGE);
    }

    private void writeIntegerDeltaOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(DeltaExchangeOperator.DELTA_EXCHANGE,
                new Attribute[]{
                        new Attribute.Default<String>(DeltaExchangeOperator.DELTA, Integer.toString((int) operator.tuning)),
                        new Attribute.Default<String>("integer", "true"),
                        getWeightAttribute(operator.weight),
                        new Attribute.Default<String>("autoOptimize", "false")
                }
        );
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(DeltaExchangeOperator.DELTA_EXCHANGE);
    }

    private void writeSwapOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(SwapOperator.SWAP_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<String>("size", Integer.toString((int) operator.tuning)),
                        getWeightAttribute(operator.weight),
                        new Attribute.Default<String>("autoOptimize", "false")
                }
        );
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(SwapOperator.SWAP_OPERATOR);
    }

    private void writeBitFlipOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(BitFlipOperator.BIT_FLIP_OPERATOR,
                getWeightAttribute(operator.weight));
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(BitFlipOperator.BIT_FLIP_OPERATOR);
    }

    private void writeTreeBitMoveOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(TreeBitMoveOperator.BIT_MOVE_OPERATOR,
                        getWeightAttribute(operator.weight));
        writer.writeIDref(TreeModel.TREE_MODEL,  modelPrefix + TreeModel.TREE_MODEL);
        writer.writeCloseTag(TreeBitMoveOperator.BIT_MOVE_OPERATOR);
    }

    private void writeUniformOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag("uniformOperator",
                getWeightAttribute(operator.weight));
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag("uniformOperator");
    }

    private void writeIntegerUniformOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag("uniformIntegerOperator",
                getWeightAttribute(operator.weight));
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag("uniformIntegerOperator");
    }

    private void writeNarrowExchangeOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(ExchangeOperator.NARROW_EXCHANGE,
                getWeightAttribute(operator.weight));
        writer.writeIDref(TreeModel.TREE_MODEL,  modelPrefix + TreeModel.TREE_MODEL);
        writer.writeCloseTag(ExchangeOperator.NARROW_EXCHANGE);
    }

    private void writeWideExchangeOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(ExchangeOperator.WIDE_EXCHANGE,
                getWeightAttribute(operator.weight));
        writer.writeIDref(TreeModel.TREE_MODEL,  modelPrefix + TreeModel.TREE_MODEL);
        writer.writeCloseTag(ExchangeOperator.WIDE_EXCHANGE);
    }

    private void writeWilsonBaldingOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(WilsonBalding.WILSON_BALDING,
                getWeightAttribute(operator.weight));
        writer.writeIDref(TreeModel.TREE_MODEL,  modelPrefix + TreeModel.TREE_MODEL);
        // not supported anymore. probably never worked. (todo) get it out of GUI too
//        if (options.nodeHeightPrior == TreePrior.CONSTANT) {
//            treePriorGenerator.writeNodeHeightPriorModelRef(writer);
//        }
        writer.writeCloseTag(WilsonBalding.WILSON_BALDING);
    }

    private void writeSampleNonActiveOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(SampleNonActiveGibbsOperator.SAMPLE_NONACTIVE_GIBBS_OPERATOR,
                getWeightAttribute(operator.weight));

        writer.writeOpenTag(SampleNonActiveGibbsOperator.DISTRIBUTION);
        writeParameterRefByName(writer, operator.getName());
        writer.writeCloseTag(SampleNonActiveGibbsOperator.DISTRIBUTION);

        writer.writeOpenTag(SampleNonActiveGibbsOperator.DATA_PARAMETER);
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(SampleNonActiveGibbsOperator.DATA_PARAMETER);

        writer.writeOpenTag(SampleNonActiveGibbsOperator.INDICATOR_PARAMETER);
        writeParameterRefByName(writer, operator.parameter2.getName());
        writer.writeCloseTag(SampleNonActiveGibbsOperator.INDICATOR_PARAMETER);

        writer.writeCloseTag(SampleNonActiveGibbsOperator.SAMPLE_NONACTIVE_GIBBS_OPERATOR);
    }

    private void writeGMRFGibbsOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(
                GMRFSkyrideBlockUpdateOperator.BLOCK_UPDATE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(GMRFSkyrideBlockUpdateOperator.SCALE_FACTOR, operator.tuning),
                        getWeightAttribute(operator.weight)
                }
        );
        writer.writeIDref(GMRFSkyrideLikelihood.SKYLINE_LIKELIHOOD,  modelPrefix + "skyride");
        writer.writeCloseTag(GMRFSkyrideBlockUpdateOperator.BLOCK_UPDATE_OPERATOR);
    }

    private void writeScaleWithIndicatorsOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(
                ScaleOperator.SCALE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperator.SCALE_FACTOR, operator.tuning),
                        getWeightAttribute(operator.weight)
                });
        writeParameter1Ref(writer, operator);
        writer.writeOpenTag(ScaleOperator.INDICATORS, new Attribute.Default<String>(ScaleOperator.PICKONEPROB, "1.0"));
        writeParameterRefByName(writer, operator.parameter2.getName());
        writer.writeCloseTag(ScaleOperator.INDICATORS);
        writer.writeCloseTag(ScaleOperator.SCALE_OPERATOR);
    }

    private void writeSubtreeSlideOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(SubtreeSlideOperator.SUBTREE_SLIDE,
                new Attribute[]{
                        new Attribute.Default<Double>("size", operator.tuning),
                        new Attribute.Default<String>("gaussian", "true"),
                        getWeightAttribute(operator.weight)
                }
        );
        writer.writeIDref(TreeModel.TREE_MODEL,  modelPrefix + TreeModel.TREE_MODEL);
        writer.writeCloseTag(SubtreeSlideOperator.SUBTREE_SLIDE);
    }
    
    private void writeSpeciesTreeOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(TreeNodeSlide.TREE_NODE_REHEIGHT,
                new Attribute[]{ getWeightAttribute(operator.weight) }
        );
        writer.writeIDref(TraitGuesser.Traits.TRAIT_SPECIES.toString(),  TraitGuesser.Traits.TRAIT_SPECIES.toString());
        writer.writeIDref(SpeciesTreeModel.SPECIES_TREE,  Generator.SP_TREE);
        writer.writeCloseTag(TreeNodeSlide.TREE_NODE_REHEIGHT);
    }

    private Attribute getWeightAttribute(double weight) {
        if (weight == (int)weight) {
            return new Attribute.Default<Integer>("weight", (int)weight);
        } else {
            return new Attribute.Default<Double>("weight", weight);
        }
    }
}
