/*
 * RatioMasker.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.VariableListener;
import dr.util.Transform;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class RatioMasker implements VariableListener {

    private final Parameter ratio;
    private final Parameter mask;
    private final NodeHeightTransform nodeHeightTransform;
    private final double ratioSamllValueThreshold;
    private final double heightDistanceThreshold;
    private boolean updatedByHeight;
    private boolean updatedByRatio;

    public RatioMasker(Parameter ratio,
                       Parameter mask,
                       NodeHeightTransform nodeHeightTransform,
                       double ratioSamllValueThreshold,
                       double heightDistanceThreshold) {
        this.ratio = ratio;
        this.mask = mask;
        this.nodeHeightTransform = nodeHeightTransform;
        this.ratioSamllValueThreshold = ratioSamllValueThreshold;
        this.heightDistanceThreshold = heightDistanceThreshold;
        ratio.addParameterListener(this);
        dimensionCheck();
    }

    private void dimensionCheck() {
        if (ratio.getDimension() != mask.getDimension()) {
            throw new RuntimeException("Ratio and mask parameters should have same dimension.");
        }
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        updateMask();
    }

    private void updateMask() {
        double[] maskByHeight = updateMaskByHeight();
        double[] maskByRatio = updateMaskByRatio();
        for (int i = 0; i < mask.getDimension(); i++) {
            final double currentMaskValue = maskByHeight[i] == 0.0 || maskByRatio[i] == 0.0 ? 0.0:1.0;
            mask.setParameterValueQuietly(i, currentMaskValue);
        }
        mask.fireParameterChangedEvent();
    }

    private double[] updateMaskByHeight() {
        if (!updatedByHeight) {
            updatedByHeight = true;
            return nodeHeightTransform.getNodeHeightTransformDelegate().setMaskByHeightDifference(heightDistanceThreshold);
        }
        return mask.getParameterValues();
    }

    private double[] updateMaskByRatio() {
        if (!updatedByRatio) {
            double[] maskByRatio = new double[ratio.getDimension()];
            for (int i = 0; i < ratio.getDimension(); i++) {
                if (ratio.getParameterValue(i) > ratioSamllValueThreshold && ratio.getParameterValue(i) < 1.0 - ratioSamllValueThreshold) {
                    maskByRatio[i] = 1.0;
                }
            }
            updatedByRatio = true;
            return maskByRatio;
        }
        return mask.getParameterValues();
    }


// **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        final static String RATIO_MASKER = "ratioMasker";
        final static String RATIO_THRESHOLD = "ratioThreshold";
        final static String HEIGHT_THRESHOLD = "heightThreshold";
        final static String MASK = "mask";
        final static String RATIO = "ratio";

        public String getParserName() {
            return RATIO_MASKER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter mask = (Parameter) xo.getChild(MASK).getChild(Parameter.class);
            Parameter ratio = (Parameter) xo.getChild(RATIO).getChild(Parameter.class);

            if (mask.getDimension() == 1) {
                mask.setDimension(ratio.getDimension());
            }
            double ratioThreshold = xo.getAttribute(RATIO_THRESHOLD, 0.0);
            double heightThreshold = xo.getAttribute(HEIGHT_THRESHOLD, 0.0);
            Transform.ComposeMultivariable transform = (Transform.ComposeMultivariable) xo.getChild(Transform.ComposeMultivariable.class);

            return new RatioMasker(ratio, mask, (NodeHeightTransform) transform.getInnerTransform(), ratioThreshold, heightThreshold);
        }

        public String getParserDescription() {
            return "A utility to craft mask for filtering dimensions in the ratio space for nodeHeight transform";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(MASK, Parameter.class),
                new ElementRule(RATIO, Parameter.class),
                AttributeRule.newDoubleRule(RATIO_THRESHOLD),
                AttributeRule.newDoubleRule(HEIGHT_THRESHOLD),
                new ElementRule(Transform.ComposeMultivariable.class)
        };

        public Class getReturnType() {
            return RatioMasker.class;
        }
    };

}
