package dr.evomodel.antigenic;

import dr.inference.model.*;
import dr.math.MathUtils;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.DataTable;
import dr.xml.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public abstract class AntigenicTraitLikelihood extends MultidimensionalScalingLikelihood implements Citable {


    public AntigenicTraitLikelihood(String name) {
        super(name);
    }

    protected void initalizeTable(DataTable<String[]> dataTable, double[][] observationValueTable, ObservationType[][] observationTypeTable, boolean log2Transform) {
        // the largest measured value for any given column of data
        double[] maxColumnValue = new double[dataTable.getColumnCount()];
        // the largest measured value over all
        double maxAssayValue = 0.0;

        for (int i = 0; i < dataTable.getRowCount(); i++) {
            String[] dataRow = dataTable.getRow(i);

            for (int j = 0; j < dataTable.getColumnCount(); j++) {
                Double value = null;
                ObservationType type = null;

                if (dataRow[j].startsWith("<")) {
                    // is a lower bound
                    value = convertString(dataRow[j].substring(1));
                    if (Double.isNaN(value)) {
                        throw new RuntimeException("Illegal value in table as a threshold");
                    }
                    type = ObservationType.LOWER_BOUND;
                } else if (dataRow[j].startsWith(">")) {
                    // is a lower bound
                    value = convertString(dataRow[j].substring(1));
                    if (Double.isNaN(value)) {
                        throw new RuntimeException("Illegal value in table as a threshold");
                    }
                    type = ObservationType.UPPER_BOUND;
                } else {
                    value = convertString(dataRow[j]);
                    type = Double.isNaN(value) ? ObservationType.MISSING : ObservationType.POINT;
                }

                observationValueTable[i][j] = value;
                observationTypeTable[i][j] = type;

                if (!Double.isNaN(value)) {
                    if (value > maxColumnValue[j]) {
                        maxColumnValue[j] = value;
                    }
                    if (value > maxAssayValue) {
                        maxAssayValue = value;
                    }

                }

            }
        }

        if (log2Transform) {
            // transform and normalize the data...
            for (int i = 0; i < dataTable.getRowCount(); i++) {
                for (int j = 0; j < dataTable.getColumnCount(); j++) {
                    observationValueTable[i][j] = transform(observationValueTable[i][j], maxColumnValue[j], 2);
                    // the transformation reverses the bounds
                    if (observationTypeTable[i][j] == ObservationType.UPPER_BOUND) {
                        observationTypeTable[i][j] = ObservationType.LOWER_BOUND;
                    } else if (observationTypeTable[i][j] == ObservationType.LOWER_BOUND) {
                        observationTypeTable[i][j] = ObservationType.UPPER_BOUND;
                    }
                }
            }
        }
    }

    private double convertString(String value) {
        try {
            return java.lang.Double.valueOf(value);
        } catch (NumberFormatException nfe) {
            return java.lang.Double.NaN;
        }
    }

    protected double transform(final double value, final double maxValue, final double base) {
        // log2(maxValue / value)
        return (Math.log(maxValue) - Math.log(value)) / Math.log(base);
    }


}
