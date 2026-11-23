/*
 * ChainDiagnosticColumn.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.loggers;

/**
 * A column for a screen log that shows a convergence diagnostic of the chains of a run. There
 * is one value for the run, not one for each chain: its value is supplied by the combined
 * screen log once every chain has reached a state.
 *
 * @author Andrew Rambaut
 */
public class ChainDiagnosticColumn extends NumberColumn implements Loggable {

    public enum Type {
        MIN_ESS("minESS", 0),
        MAX_RHAT("maxRhat", 3),
        ASDSF("ASDSF", 4),
        BURNIN("burnin", 0);

        Type(String label, int decimalPlaces) {
            this.label = label;
            this.decimalPlaces = decimalPlaces;
        }

        public String getLabel() {
            return label;
        }

        private final String label;
        private final int decimalPlaces;
    }

    /**
     * @param showName if true the name of the logged value that gives the diagnostic follows it
     */
    public ChainDiagnosticColumn(Type type, boolean showName) {
        super(type.label);
        this.type = type;
        this.showName = showName;
        setDecimalPlaces(type.decimalPlaces);
    }

    public Type getType() {
        return type;
    }

    /**
     * @return the diagnostic as it should be shown ("-" if it is not available)
     */
    public String format(ChainDiagnostics diagnostics) {
        double value = Double.NaN;
        String name = null;
        if (diagnostics != null) {
            switch (type) {
                case MIN_ESS:
                    value = diagnostics.getMinESS();
                    name = diagnostics.getMinESSLabel();
                    break;
                case MAX_RHAT:
                    value = diagnostics.getMaxRhat();
                    name = diagnostics.getMaxRhatLabel();
                    break;
                case ASDSF:
                    value = diagnostics.getASDSF();
                    break;
                case BURNIN:
                    value = diagnostics.getBurnin();
                    break;
            }
        }
        if (Double.isNaN(value)) {
            return "-";
        }
        String formatted = formatValue(value);
        if (showName && name != null) {
            formatted += " (" + name.trim() + ")";
        }
        return formatted;
    }

    /**
     * @return NaN: a chain does not know the diagnostic, which needs every chain
     */
    @Override
    public double getDoubleValue() {
        return Double.NaN;
    }

    @Override
    public LogColumn[] getColumns() {
        return new LogColumn[]{this};
    }

    private final Type type;
    private final boolean showName;
}
