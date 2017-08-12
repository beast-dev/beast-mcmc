/*
 * Generator.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.*;
import dr.app.beauti.util.XMLWriter;
import dr.evoxml.MergePatternsParser;
import dr.evoxml.SitePatternsParser;
import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.distribution.GammaDistributionModel;
import dr.inference.model.ParameterParser;
import dr.inferencexml.distribution.DistributionModelParser;
import dr.inferencexml.distribution.LogNormalDistributionModelParser;
import dr.inferencexml.distribution.NormalDistributionModelParser;
import dr.inferencexml.distribution.UniformDistributionModelParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public abstract class Generator {

    protected static final String COALESCENT = "coalescent";
    public static final String SP_TREE = "sptree";
    protected static final String SP_START_TREE = "spStartingTree";
    protected static final String SPECIATION_LIKE = "speciation.likelihood";
    public static final String SPLIT_POPS = "splitPopSize";
    protected static final String PDIST = "pdist";
    //	protected static final String STP = "stp";
    protected static final String SPOPS = TraitData.TRAIT_SPECIES + "." + "popSizesLikelihood";

    protected final BeautiOptions options;

    protected Generator(BeautiOptions options) {
        this.options = options;
    }

    public Generator(BeautiOptions options, ComponentFactory[] components) {
        this.options = options;
        if (components != null) {
            for (ComponentFactory component : components) {
                this.components.add(component.createGenerator(options));
            }
        }
    }

    public final void checkComponentOptions() throws GeneratorException {
        for (ComponentGenerator component : components) {
                component.checkOptions();
        }
    }

    /**
     * fix a parameter
     *
     * @param id    the id
     * @param value the value
     */
//    public void fixParameter(Parameter parameter, double value) {
////        dr.app.beauti.options.Parameter parameter = options.getParameter(id);
//        if (parameter == null) {
//            throw new IllegalArgumentException("parameter with name, " + parameter.getName() + ", is unknown");
//        }
//        parameter.isFixed = true;
//        parameter.initial = value;
//    }


    /**
     * write a parameter
     *
     * @param wrapperName wrapperName
     * @param id          the id
     * @param writer      the writer
     */
    public void writeParameterRef(String wrapperName, String id, XMLWriter writer) {
        writer.writeOpenTag(wrapperName);
        writer.writeIDref(ParameterParser.PARAMETER, id);
        writer.writeCloseTag(wrapperName);
    }

    public void writeParameterRef(String id, XMLWriter writer) {
        writer.writeIDref(ParameterParser.PARAMETER, id);
    }

    /**
     * write a parameter
     *
     * @param id      the id
     * @param options PartitionOptions
     * @param writer  the writer
     */
    public void writeParameter(String id, PartitionOptions options, XMLWriter writer) {
        Parameter parameter = options.getParameter(id);
        String prefix = options.getPrefix();

        if (parameter == null) {
            throw new IllegalArgumentException("parameter with name, " + id + ", is unknown; and its prefix is " + options.getPrefix());
        }
        writeParameter(prefix + id, parameter, writer);
    }


    public void writeParameter(int num, String id, PartitionSubstitutionModel model, XMLWriter writer) {
        Parameter parameter = model.getParameter(model.getPrefixCodon(num) + id);
        String prefix = model.getPrefix(num);

        if (parameter == null) {
            throw new IllegalArgumentException("parameter with name, " + id + ", is unknown; and its prefix is " + model.getPrefix());
        }
        writeParameter(prefix + id, parameter, writer);
    }

    /**
     * write a parameter
     *
     * @param wrapperName wrapperName
     * @param id          the id
     * @param dimension   dimension
     * @param writer      the writer
     */
    public void writeParameter(String wrapperName, String id, int dimension, XMLWriter writer) {
        Parameter parameter = options.getParameter(id);
        writer.writeOpenTag(wrapperName);
        writeParameter(parameter, dimension, writer);
        writer.writeCloseTag(wrapperName);
    }

    /**
     * write a parameter
     *
     * @param num         num
     * @param wrapperName wrapperName
     * @param id          the id
     * @param model       PartitionSubstitutionModel
     * @param writer      the writer
     */
    public void writeParameter(int num, String wrapperName, String id, PartitionSubstitutionModel model, XMLWriter writer) {
        writer.writeOpenTag(wrapperName);
        writeParameter(num, id, model, writer);
        writer.writeCloseTag(wrapperName);
    }

    public void writeParameter(String wrapperName, String id, PartitionOptions options, XMLWriter writer) {
        writer.writeOpenTag(wrapperName);
        writeParameter(id, options, writer);
        writer.writeCloseTag(wrapperName);
    }

    /**
     * write a parameter
     *
     * @param parameter the parameter
     * @param dimension the dimension
     * @param writer    the writer
     */
    public void writeParameter(Parameter parameter, int dimension, XMLWriter writer) {
//        Parameter parameter = options.getParameter(id);
        if (parameter == null) {
            throw new IllegalArgumentException("parameter (== null) is unknown");
        }
        if (parameter.isFixed()) { // with prefix
            writeParameter(parameter.getName(), dimension, parameter.getInitial(), Double.NaN, Double.NaN, writer);
        } else {
            double lower = Double.NaN;
            double upper = Double.NaN;
            if (parameter.isNonNegative) {
                lower = 0.0;
            }
            if (parameter.isZeroOne) {
                lower = 0.0;
                upper = 1.0;
            }
            writeParameter(parameter.getName(), dimension, parameter.getInitial(), lower, upper, writer);
        }
    }

    /**
     * write a parameter
     *
     * @param parameterId     the parameter name/id
     * @param parameterColumn the parameter column from which the samples are taken
     * @param fileName        the file from which the samples are taken
     * @param burnin          the number of samples to be discarded
     * @param writer          the writer
     */
    public void writeParameter(String parameterId, String parameterColumn, String fileName, int burnin, XMLWriter writer) {
        ArrayList<Attribute.Default> attributes = new ArrayList<Attribute.Default>();
        attributes.add(new Attribute.Default<String>(XMLParser.ID, parameterId));

        attributes.add(new Attribute.Default<String>("parameterColumn", parameterColumn));
        attributes.add(new Attribute.Default<String>("fileName", fileName));
        attributes.add(new Attribute.Default<String>("burnin", "" + burnin));

        Attribute[] attrArray = new Attribute[attributes.size()];
        for (int i = 0; i < attrArray.length; i++) {
            attrArray[i] = attributes.get(i);
        }

        writer.writeTag(ParameterParser.PARAMETER, attrArray, true);
    }

    public void writeParameter(String id, Parameter parameter, XMLWriter writer) {
        if (parameter.isFixed()) {
            writeParameter(id, 1, parameter.getInitial(), Double.NaN, Double.NaN, writer);
        } else {
            double lower = Double.NaN;
            double upper = Double.NaN;
            if (parameter.isNonNegative) {
                lower = 0.0;
            }
            if (parameter.isZeroOne) {
                lower = 0.0;
                upper = 1.0;
            }
            writeParameter(id, 1, parameter.getInitial(), lower, upper, writer);
        }
    }


    /**
     * write a parameter
     *
     * @param id        the id
     * @param dimension the dimension
     * @param value     the value
     * @param lower     the lower bound
     * @param upper     the upper bound
     * @param writer    the writer
     */
    public void writeParameter(String id, int dimension, double value, double lower, double upper, XMLWriter writer) {
        ArrayList<Attribute.Default> attributes = new ArrayList<Attribute.Default>();
        if (id != null && id.length() > 0) {
            attributes.add(new Attribute.Default<String>(XMLParser.ID, id));
        }
        if (dimension > 1) {
            attributes.add(new Attribute.Default<String>(ParameterParser.DIMENSION, dimension + ""));
        }
        if (!Double.isNaN(value)) {
            attributes.add(new Attribute.Default<String>(ParameterParser.VALUE, multiDimensionValue(dimension, value)));
        }
        if (!Double.isNaN(lower)) {
            attributes.add(new Attribute.Default<String>(ParameterParser.LOWER, multiDimensionValue(dimension, lower)));
        }
        if (!Double.isNaN(upper)) {
            attributes.add(new Attribute.Default<String>(ParameterParser.UPPER, multiDimensionValue(dimension, upper)));
        }

        Attribute[] attrArray = new Attribute[attributes.size()];
        for (int i = 0; i < attrArray.length; i++) {
            attrArray[i] = attributes.get(i);
        }

        writer.writeTag(ParameterParser.PARAMETER, attrArray, true);
    }

    /**
     * write a parameter
     *
     * @param wrapperName wrapperName
     * @param id          the id
     * @param dimension   the dimension
     * @param value       the value
     * @param lower       the lower bound
     * @param upper       the upper bound
     * @param writer      the writer
     */
    public void writeParameter(String wrapperName, String id, int dimension, double value, double lower, double upper, XMLWriter writer) {
        writer.writeOpenTag(wrapperName);
        writeParameter(id, dimension, value, lower, upper, writer);
        writer.writeCloseTag(wrapperName);
    }

    protected void writeCodonPatternsRef(String prefix, int num, int codonPartitionCount, XMLWriter writer) {
        if (codonPartitionCount == 2 && num == 1) { // "11" of "112", num start from 1
            writer.writeIDref(MergePatternsParser.MERGE_PATTERNS, prefix + SitePatternsParser.PATTERNS);
        } else { // "2" of "112" and "123"
            writer.writeIDref(SitePatternsParser.PATTERNS, prefix + SitePatternsParser.PATTERNS);
        }
    }

    private String multiDimensionValue(int dimension, double value) {
        String multi = "";

        multi += value + "";

        // AR: A multidimensional parameter only needs to give initial values for every dimension
        // if they are actually different. A single value will automatically be expanded to every
        // dimension and make for a cleaner looking XML (and more robust to changes in the number
        // of groups/taxa etc.

//        for (int i = 2; i <= dimension; i++)
//            multi += " " + value;

        return multi;
    }


    /**
     * Write the distribution for *DistributionModel
     *
     * @param parameter the parameter
     * @param isRef     only work for uniform dist
     * @param writer    the writer
     */
    protected void writeDistribution(Parameter parameter, boolean isRef, XMLWriter writer) {

        switch (parameter.priorType) {
            case UNIFORM_PRIOR:
                String id = parameter.taxaId + "-uniformDist";
                if (isRef) {
                    writer.writeIDref(UniformDistributionModelParser.UNIFORM_DISTRIBUTION_MODEL, id);
                } else {
                    writer.writeOpenTag(UniformDistributionModelParser.UNIFORM_DISTRIBUTION_MODEL,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, id)
                            });
                    writer.writeOpenTag(UniformDistributionModelParser.LOWER);
                    writer.writeText(Double.toString(parameter.uniformLower));
                    writer.writeCloseTag(UniformDistributionModelParser.LOWER);

                    writer.writeOpenTag(UniformDistributionModelParser.UPPER);
                    writer.writeText(Double.toString(parameter.uniformUpper));
                    writer.writeCloseTag(UniformDistributionModelParser.UPPER);

                    writer.writeCloseTag(UniformDistributionModelParser.UNIFORM_DISTRIBUTION_MODEL);
                }
                break;
            case EXPONENTIAL_PRIOR:
                writer.writeOpenTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);
                writer.writeOpenTag(DistributionModelParser.MEAN);
                writer.writeText(Double.toString(parameter.mean));
                writer.writeCloseTag(DistributionModelParser.MEAN);

                writer.writeOpenTag(DistributionModelParser.OFFSET);
                writer.writeText(Double.toString(parameter.offset));
                writer.writeCloseTag(DistributionModelParser.OFFSET);
                writer.writeCloseTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);
                break;
            case NORMAL_PRIOR:
                writer.writeOpenTag(NormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);
                writer.writeOpenTag(NormalDistributionModelParser.MEAN);
                writer.writeText(Double.toString(parameter.mean));
                writer.writeCloseTag(NormalDistributionModelParser.MEAN);

                writer.writeOpenTag(NormalDistributionModelParser.STDEV);
                writer.writeText(Double.toString(parameter.stdev));
                writer.writeCloseTag(NormalDistributionModelParser.STDEV);
                writer.writeCloseTag(NormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);
                break;
            case LOGNORMAL_PRIOR:
                writer.writeOpenTag(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL,
                        new Attribute[]{
                                new Attribute.Default<Boolean>(LogNormalDistributionModelParser.MEAN_IN_REAL_SPACE, parameter.isInRealSpace()),
                                new Attribute.Default<Boolean>(LogNormalDistributionModelParser.STDEV_IN_REAL_SPACE, parameter.isInRealSpace())
                        });
                writer.writeOpenTag(LogNormalDistributionModelParser.MEAN);
                writer.writeText(Double.toString(parameter.mean));
                writer.writeCloseTag(LogNormalDistributionModelParser.MEAN);

                writer.writeOpenTag(LogNormalDistributionModelParser.STDEV);
                writer.writeText(Double.toString(parameter.stdev));
                writer.writeCloseTag(LogNormalDistributionModelParser.STDEV);

                writer.writeOpenTag(LogNormalDistributionModelParser.OFFSET);
                writer.writeText(Double.toString(parameter.offset));
                writer.writeCloseTag(LogNormalDistributionModelParser.OFFSET);
                writer.writeCloseTag(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL);
                break;
            case GAMMA_PRIOR:
                writer.writeOpenTag(GammaDistributionModel.GAMMA_DISTRIBUTION_MODEL);
                writer.writeOpenTag(DistributionModelParser.SHAPE);
                writer.writeText(Double.toString(parameter.shape));
                writer.writeCloseTag(DistributionModelParser.SHAPE);

                writer.writeOpenTag(DistributionModelParser.SCALE);
                writer.writeText(Double.toString(parameter.scale));
                writer.writeCloseTag(DistributionModelParser.SCALE);

                writer.writeOpenTag(DistributionModelParser.OFFSET);
                writer.writeText(Double.toString(parameter.offset));
                writer.writeCloseTag(DistributionModelParser.OFFSET);
                writer.writeCloseTag(GammaDistributionModel.GAMMA_DISTRIBUTION_MODEL);
                break;
            default:
                throw new IllegalArgumentException("Unknown Distribution Model for " + parameter.getName());
        }
    }

    public void writeReferenceComment(String[] lines, XMLWriter writer) {
        for (String line : lines)
            writer.writeComment(line);
    }

    public void generateInsertionPoint(final ComponentGenerator.InsertionPoint ip, final XMLWriter writer) {
        generateInsertionPoint(ip, null, writer);
    }

    public void generateInsertionPoint(final ComponentGenerator.InsertionPoint ip, final Object item, final XMLWriter writer) {
        for (ComponentGenerator component : components) {
            if (component.usesInsertionPoint(ip)) {
                component.generateAtInsertionPoint(this, ip, item, writer);
            }
        }
    }

    public void generateInsertionPoint(final ComponentGenerator.InsertionPoint ip, final Object item, final String prefix, final XMLWriter writer) {
        for (ComponentGenerator component : components) {
            if (component.usesInsertionPoint(ip)) {
                component.generateAtInsertionPoint(this, ip, item, prefix, writer);
            }
        }
    }

    private final List<ComponentGenerator> components = new ArrayList<ComponentGenerator>();

    public class GeneratorException extends Exception {
        public GeneratorException(String message) {
            super(message);
            switchToPanel = null;
        }

        public GeneratorException(String message, String switchToPanel) {
            super(message);
            this.switchToPanel = switchToPanel;
        }

        public String getSwitchToPanel() {
            return switchToPanel;
        }

        private final String switchToPanel;
    }
}
