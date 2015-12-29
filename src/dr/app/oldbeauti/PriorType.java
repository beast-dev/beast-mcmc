/*
 * PriorType.java
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

package dr.app.oldbeauti;

import java.text.NumberFormat;

/**
 * @author Alexei Drummond
 */
public enum PriorType {

    NONE,
    UNIFORM_PRIOR,
    EXPONENTIAL_PRIOR,
    NORMAL_PRIOR,
    LOGNORMAL_PRIOR,
    GAMMA_PRIOR,
    JEFFREYS_PRIOR,
    POISSON_PRIOR,
    TRUNC_NORMAL_PRIOR;

    public String toString() {

        switch (this) {
            case NONE: return "none";
            case UNIFORM_PRIOR: return "Uniform";
            case EXPONENTIAL_PRIOR: return "Exponential";
            case NORMAL_PRIOR: return "Normal";
            case LOGNORMAL_PRIOR: return "Lognormal";
            case GAMMA_PRIOR: return "Gamma";
            case JEFFREYS_PRIOR: return "Jeffreys'";
            case POISSON_PRIOR: return "Poisson";
            case TRUNC_NORMAL_PRIOR: return "Truncated Normal";
            default: return "";
        }
    }

    public String getPriorString(BeastGenerator.Parameter param) {

        NumberFormat formatter = NumberFormat.getNumberInstance();
        StringBuffer buffer = new StringBuffer();

        if (!param.priorEdited) {
            buffer.append("* ");
        }
        switch (param.priorType) {
            case NONE:
                buffer.append("Using Tree Prior");
                break;
            case UNIFORM_PRIOR:
                if (!param.isDiscrete && !param.isStatistic) {
                    buffer.append("Uniform [");
                    buffer.append(formatter.format(param.uniformLower));
                    buffer.append(", ");
                    buffer.append(formatter.format(param.uniformUpper));
                    buffer.append("]");
                } else {
                    buffer.append("Uniform");
                }
                break;
            case EXPONENTIAL_PRIOR:
                buffer.append("Exponential [");
                buffer.append(formatter.format(param.exponentialMean));
                buffer.append("]");
                break;
            case NORMAL_PRIOR:
                buffer.append("Normal [");
                buffer.append(formatter.format(param.normalMean));
                buffer.append(", ");
                buffer.append(formatter.format(param.normalStdev));
                buffer.append("]");
                break;
            case LOGNORMAL_PRIOR:
                buffer.append("LogNormal [");
                buffer.append(formatter.format(param.logNormalMean));
                buffer.append(", ");
                buffer.append(formatter.format(param.logNormalStdev));
                buffer.append("]");
                break;
            case GAMMA_PRIOR:
                buffer.append("Gamma [");
                buffer.append(formatter.format(param.gammaAlpha));
                buffer.append(", ");
                buffer.append(formatter.format(param.gammaBeta));
                buffer.append("]");
                break;
            case JEFFREYS_PRIOR:
                buffer.append("Jeffreys");
                break;
            case POISSON_PRIOR:
                buffer.append("Poisson [");
                buffer.append(formatter.format(param.poissonMean));
                buffer.append("]");
                break;
            case TRUNC_NORMAL_PRIOR:
                buffer.append("Truncated Normal [");
                buffer.append(formatter.format(param.normalMean));
                buffer.append(", ");
                buffer.append(formatter.format(param.normalStdev));
                buffer.append("]");
                buffer.append(" in [");
                buffer.append(formatter.format(param.uniformLower));
                buffer.append(", ");
                buffer.append(formatter.format(param.uniformUpper));
                buffer.append("]");

                break;
            default:
                throw new IllegalArgumentException("Unknown prior type");
        }
        if (param.priorType != PriorType.NONE && !param.isStatistic) {
            buffer.append(", initial=").append(param.initial);
        }

        return buffer.toString();
    }
}
