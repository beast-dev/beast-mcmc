package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.xml.*;

import java.util.Map;
import java.util.TreeMap;

/**
 * @author Marc A. Suchard
 */
public class LogOrderedMatrix implements Loggable {

    private static final String ORDERED_MATRIX = "orderedMatrix";
    private static final String ORDER_BY = "orderBy";

    enum OrderBy {
        COLUMN {
            @Override
            int getRow(int[] map, int row) {
                return row;
            }

            @Override
            int getCol(int[] map, int col) {
                return map[col];
            }

            @Override
            double getObjective(MatrixParameterInterface matrix, int col) {
                double objective = 0.0;
                for (int row = 0; row < matrix.getRowDimension(); ++row) {
                    double entry = matrix.getParameterValue(row, col);
                    objective += entry * entry;
                }

                return objective;
            }

            @Override
            int getCount(MatrixParameterInterface matrix) {
                return matrix.getColumnDimension();
            }
        },
        ROW {
            @Override
            int getRow(int[] map, int row) {
                return map[row];
            }

            @Override
            int getCol(int[] map, int col) {
                return col;
            }

            @Override
            double getObjective(MatrixParameterInterface matrix, int row) {
                double objective = 0.0;
                for (int col = 0; col < matrix.getColumnDimension(); ++col) {
                    double entry = matrix.getParameterValue(row, col);
                    objective += entry * entry;
                }

                return objective;
            }

            @Override
            int getCount(MatrixParameterInterface matrix) {
                return matrix.getRowDimension();
            }
        };

        abstract int getRow(int[] map, int row);

        abstract int getCol(int[] map, int col);

        abstract double getObjective(MatrixParameterInterface matrix, int index);

        abstract int getCount(MatrixParameterInterface matrix);

        static OrderBy parseOrderBy(XMLObject xo) throws XMLParseException {
            String option = xo.getAttribute(ORDER_BY, COLUMN.name());
            if (option.equalsIgnoreCase(COLUMN.name())) {
                return COLUMN;
            } else if (option.equalsIgnoreCase(ROW.name())) {
                return ROW;
            } else {
                throw new XMLParseException("Unknown ordering");
            }
        }
    }

    private LogOrderedMatrix(MatrixParameterInterface matrix, OrderBy orderBy) {
        this.matrix = matrix;
        this.orderBy = orderBy;
    }

    @Override
    public LogColumn[] getColumns() {

        LogColumn[] logs = new LogColumn[matrix.getRowDimension() * matrix.getColumnDimension()];

        int index = 0;
        for (int col = 0; col < matrix.getColumnDimension(); ++col) {
            for (int row = 0; row < matrix.getRowDimension(); ++row) {
                logs[index] = makeColumn(row, col);
                ++index;
            }
        }
        return logs;
    }

    private NumberColumn makeColumn(final int row, final int col) {
        return new NumberColumn(getName(row, col)) {
            @Override
            public double getDoubleValue() {
                if (row == 0 && col == 0) {
                    reorderMatrix();
                }
                return matrix.getParameterValue(
                        orderBy.getRow(orderMap, row),
                        orderBy.getCol(orderMap, col));
            }
        };
    }

    private void reorderMatrix() {

        Map<Double, Integer> tree = new TreeMap<Double, Integer>();

        for (int index = 0; index < orderBy.getCount(matrix); ++index) {
            tree.put(-orderBy.getObjective(matrix, index), index);
        }
        
        orderMap = new int[tree.size()];
        int index = 0;
        for (Map.Entry<Double, Integer> entry : tree.entrySet()) {
            orderMap[index] = entry.getValue();
            ++index;
        }

        if (DEBUG) {
            System.err.println("Order by: " + orderBy.name());
            int idx = 0;
            for (Map.Entry<Double, Integer> entry : tree.entrySet()) {
                System.err.println("\t" + entry.getKey() + " " + entry.getValue());
                ++idx;
            }
            System.err.print("\t\t");
            for (int i : orderMap) {
                System.err.print(" " + i);
            }
            System.err.println();
        }                    
    }
    
    private String getName(int row, int col) {
        return "ordered." + matrix.getParameterName() + "." + row + "." + col;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

          public String getParserName() {
              return ORDERED_MATRIX;
          }

          public Object parseXMLObject(XMLObject xo) throws XMLParseException {

              MatrixParameterInterface matrix = (MatrixParameterInterface)
                      xo. getChild(MatrixParameterInterface.class);

              OrderBy orderBy = OrderBy.parseOrderBy(xo);

              return new LogOrderedMatrix(matrix, orderBy);
          }

          //************************************************************************
          // AbstractXMLObjectParser implementation
          //************************************************************************

          public String getParserDescription() {
              return "An ordered matrix logger.";
          }

          public XMLSyntaxRule[] getSyntaxRules() {
              return rules;
          }

          private final XMLSyntaxRule[] rules = {
                  new ElementRule(MatrixParameterInterface.class, 1, Integer.MAX_VALUE),
                  AttributeRule.newStringRule(ORDER_BY),
          };

          public Class getReturnType() {
              return LogOrderedMatrix.class;
          }
      };
    final private MatrixParameterInterface matrix;
    final private OrderBy orderBy;

    private int[] orderMap;

    private static final boolean DEBUG = false;
}
