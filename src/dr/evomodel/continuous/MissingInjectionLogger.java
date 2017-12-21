/*
 * MissingInjectionLogger.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.continuous;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.xml.*;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class MissingInjectionLogger implements Loggable, Reportable {

    private static final String MISSING_INJECTION_LOGGER = "injectedMissingTraitsLogger";
    private static final String ADJUSTMENT = "adjustment";

    private final TreeTrait traitProvider;
    private final List<MissingInjection.TaxonInformation> taxonInformation;
    private final Tree tree;
    private final Adjustment adjustment;

    public enum Adjustment {
        DIFFERENCE("difference") {
            @Override
            double adjustedValue(double newValue, double oldValue) {
                return newValue - oldValue;
            }
        },
        ORIGINAL("original") {
            @Override
            double adjustedValue(double newValue, double oldValue) {
                return oldValue;
            }
        },
        RAW("raw"){
            @Override
            double adjustedValue(double newValue, double oldValue) {
                return newValue;
            }
        };

        Adjustment(String name) {
            this.name = name;
        }

        private String name;

        public String getName() {
            return name;
        }

        public static Adjustment parse(String name) {
            name = name.toLowerCase();
            for (Adjustment adjustment : Adjustment.values()) {
                if (name.compareTo(adjustment.getName()) == 0) {
                    return adjustment;
                }
            }
            throw new IllegalArgumentException("Unknown adjustment type");
        }

        abstract double adjustedValue(double newValue, double oldValue);
    }

    private MissingInjectionLogger(MissingInjection injector,
                                   TreeTrait traitProvider,
                                   Tree tree,
                                   Adjustment adjustment) {

        this.taxonInformation = injector.getTaxonInformation();
        this.traitProvider = traitProvider;
        this.tree = tree;
        this.adjustment = adjustment;
    }

    private LogColumn[] logColumns = null;

    @Override
    public LogColumn[] getColumns() {

        if (logColumns == null) {
            logColumns = createLogColumns();
        }
        return logColumns;
    }

    private LogColumn[] createLogColumns() {

        int count = getNumberOfMissingValues();

        LogColumn[] columns = new LogColumn[count];

        int index = 0;
        for (final MissingInjection.TaxonInformation info : taxonInformation) {
            for (final MissingInjection.InjectedMissingValue missing : info.injectedMissingValues) {
                columns[index] = new NumberColumn(getColumnName(info, missing)) {
                    @Override
                    public double getDoubleValue() {
                        return adjustment.adjustedValue(
                                getTraitValue(info.index, missing.index),
                                missing.originalValue
                        );
                    }
                };
                ++index;
            }
        }
        
        return columns;
    }

    private double getTraitValue(int taxonIndex, int traitIndex) {
        double[] trait = (double[]) traitProvider.getTrait(tree, tree.getExternalNode(taxonIndex));
        return trait[traitIndex];
    }

    private static String getColumnName(MissingInjection.TaxonInformation info,
                                        MissingInjection.InjectedMissingValue missing) {
        return info.taxon.getId() + "." + (missing.index + 1);
    }

    private int getNumberOfMissingValues() {
        int count = 0;
        for (MissingInjection.TaxonInformation info : taxonInformation) {
            count += info.injectedMissingValues.size();
        }
        return count;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MissingInjection injector = (MissingInjection) xo.getChild(MissingInjection.class);

            TreeDataLikelihood dataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
            Tree tree = dataLikelihood.getTree();
            String traitName = injector.getTraitName();
            TreeTrait treeTrait = dataLikelihood.getTreeTrait(traitName);

            if (treeTrait == null) {
                throw new XMLParseException("Unable to find trait '" + traitName + "'");
            }

            Adjustment adjustment = Adjustment.parse(
                    xo.getAttribute(ADJUSTMENT,
                    Adjustment.DIFFERENCE.getName()));

            return new MissingInjectionLogger(injector, treeTrait, tree, adjustment);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return MissingInjectionLogger.class;
        }

        @Override
        public String getParserName() {
            return MISSING_INJECTION_LOGGER;
        }

        private final XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                new ElementRule(MissingInjection.class),
                new ElementRule(TreeDataLikelihood.class),
                AttributeRule.newStringRule(ADJUSTMENT, true),
        };
    };

    @Override
    public String getReport() {
        LogColumn[] columns = getColumns();

        StringBuilder sb = new StringBuilder();
        for (LogColumn column : columns) {
            sb.append(column.getFormatted()).append("\t");
        }
        return sb.toString();

    }
}
