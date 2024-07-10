/*
 * RealNumberField.java
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

package dr.app.gui.components;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.PlainDocument;
import java.awt.*;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;

public class RealNumberField extends JTextField implements FocusListener, DocumentListener {

    public static String NaN = "NaN";
    public static String POSITIVE_INFINITY = "+INF";
    public static String NEGATIVE_INFINITY = "-INF";
    public static String MAX_VALUE = "MAX";
    public static String MIN_VALUE = "MIN";

    protected static char MINUS = '-';
    protected static char PERIOD = '.';
    protected EventListenerList changeListeners = new EventListenerList();
    protected double min;
    protected double max;
    protected final boolean includeMin;
    protected final boolean includeMax;
    protected boolean range_check = false;
    protected boolean range_checked = false;
    protected String label; // make sensible error message

    private boolean isValueValid = true;

    protected boolean allowEmpty = false;

    public RealNumberField() { // no FocusListener
        super();
        setLabel("Value");
        includeMin = true;
        includeMax = true;
    }

    public RealNumberField(double min, double max) {
        this(min, max, "Value");
    }

    public RealNumberField(double min, double max, String label) { // no FocusListener
        this(min, true, max, true, label);
    }

    public RealNumberField(double min, boolean includeMin, double max, boolean includeMax, String label) { // no FocusListener
        super();
        this.min = min;
        this.max = max;
        this.includeMin = includeMin;
        this.includeMax = includeMax;
        setLabel(label);
        range_check = true;
        this.addFocusListener(this);
    }

    public void setAllowEmpty(boolean allowEmpty) {
        this.allowEmpty = allowEmpty;
    }

    public void focusGained(FocusEvent evt) {
    }

    public void focusLost(FocusEvent evt) {
        validateField();
    }

    public void validateField() {
        if (range_check && !range_checked) {
            range_checked = true;
            isValueValid = isValueValid();
            if (!isValueValid) {
                if (isFieldEmpty()) {
                    displayEmptyFieldErrorMessage();
                } else {
                    displayInvalidValueErrorMessage();
                }
                // regain focus for this component
                this.requestFocus();
            }
        }

    }

    public boolean isFieldEmpty() {
        return getText().trim().equals("");
    }

    public boolean isValueValid() {
        if (allowEmpty && isFieldEmpty()) {
            return true;
        }
        if (range_check) {
            try {
                if (getValue() == null) {
                    return false;
                }
                double value = getValue();
                if (value < min || value > max) {
                    return false;
                }
                if (!includeMin && value == min) {
                    return false;
                }
                if (!includeMax && value == max) {
                    return false;
                }
            } catch (NumberFormatException e) {
                return false;
            }
        }
        return true;
    }

    public void setText(Double value) {
        if (value == null && allowEmpty) {
            setText("");
            return;
        }
        if (value.isNaN()) {
            setText(NaN);
        } else if (value == Double.POSITIVE_INFINITY) {
            setText(POSITIVE_INFINITY);
        } else if (value == Double.NEGATIVE_INFINITY) {
            setText(NEGATIVE_INFINITY);
        } else if (value == Double.MAX_VALUE) {
            setText(MAX_VALUE);
        } else if (value == Double.MIN_VALUE) {
            setText(MIN_VALUE);
        } else {
            setText(Double.toString(value));
        }
    }

    public void setText(Integer obj) {
        setText(obj.toString()); // where used?
    }

    public void setText(Long obj) {
        setText(obj.toString()); // where used?
    }

    public String getErrorMessage() {
        String message = "";
        if (min == Double.MIN_VALUE) {
            message = " greater than 0";
        } else if (!Double.isInfinite(min) && min != -Double.MAX_VALUE) {
            message = " greater than " + min;
        }
        if (max == -Double.MIN_VALUE) {
            message = " less than 0";
        } else if (!Double.isInfinite(max) && max != Double.MAX_VALUE) {
            if (message.length() > 0) {
                message += " and";
            }
            message = " less than " + max;
        }

        return label + " must be" + message;
    }

    private void displayInvalidValueErrorMessage() {
        JOptionPane.showMessageDialog(null,
                getErrorMessage(), "Field value is not valid", JOptionPane.ERROR_MESSAGE);
    }

    private void displayEmptyFieldErrorMessage() {
        JOptionPane.showMessageDialog(null,
                getErrorMessage(), "Field is empty", JOptionPane.ERROR_MESSAGE);
    }

    public void setRange(double min, double max) {
        this.min = min;
        this.max = max;
        range_check = true;
    }

    public void setValue(double value) {
        if (range_check) {
            if (value < min || value > max) {
                displayInvalidValueErrorMessage();
                return;
            }
            if (!includeMin && value == min) {
                displayInvalidValueErrorMessage();
                return;
            }
            if (!includeMax && value == max) {
                displayInvalidValueErrorMessage();
                return;
            }
        }
        setText(value);
    }

    public void setLabel(String label) {
        this.label = label;
    }

    public Double getValue() {
        try {
            if (allowEmpty && getText().trim().equals("")) {
                return null;
            } else if (getText().equals(POSITIVE_INFINITY)) {
                return Double.POSITIVE_INFINITY;
            } else if (getText().equals(NEGATIVE_INFINITY)) {
                return Double.NEGATIVE_INFINITY;
            } else if (getText().equals(MAX_VALUE)) {
                return Double.MAX_VALUE;
            } else if (getText().equals(MIN_VALUE)) {
                return Double.MIN_VALUE;
            } else if (getText().equals(NaN)) {
                return Double.NaN;
            } else {
//                System.out.println("=" + getText() + "=");
                return Double.parseDouble(getText());
            }
        } catch (NumberFormatException e) {
            JOptionPane.showMessageDialog(null, "Unable to parse number correctly",
                    "Number Format Exception",
                    JOptionPane.ERROR_MESSAGE);
            isValueValid = false;
            return null;
        }
    }

    protected Document createDefaultModel() {
        Document doc = new RealNumberField.RealNumberFieldDocument();
        doc.addDocumentListener(this);
        return doc;
    }

    public void insertUpdate(DocumentEvent e) {
        fireChanged();
    }

    public void removeUpdate(DocumentEvent e) {
        fireChanged();
    }

    public void changedUpdate(DocumentEvent e) {
        fireChanged();
    }

    static char[] numberSet = {
            '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'
    };

    class RealNumberFieldDocument extends PlainDocument {
        public void insertString(int offs, String str, AttributeSet a)
                throws BadLocationException {

            if (str == "" || str == null) return;
            if (str.equals("+INF") || str.equals("-INF") || str.equals("NaN")
                    || str.equals("MAX_VALUE") || str.equals("MIN_VALUE")) {
                super.insertString(offs, str, a);
                return;
            }

            str = str.trim();

            int length = getLength();
            String buf = getText(0, offs) + str + getText(offs, length - offs);
            buf = buf.trim().toUpperCase();
            char[] array = buf.toCharArray();

            if (array.length > 0) {
                if (array[0] != MINUS && !member(array[0], numberSet) &&
                        array[0] != PERIOD) {
                    Toolkit.getDefaultToolkit().beep();
                    return;
                }
            }

            boolean period_found = (array.length > 0 && array[0] == PERIOD);
            boolean exponent_found = false;
            int exponent_index = -1;
            boolean exponent_sign_found = false;

            for (int i = 1; i < array.length; i++) {
                if (!member(array[i], numberSet)) {
                    if (!period_found && array[i] == PERIOD) {
                        period_found = true;
                    } else if (!exponent_found && array[i] == 'E') {
                        exponent_found = true;
                        exponent_index = i;
                    } else if (exponent_found && i == (exponent_index + 1) && !exponent_sign_found && array[i] == '-') {
                        exponent_sign_found = true;
                    } else {
                        Toolkit.getDefaultToolkit().beep();
                        return;
                    }
                }
            }
            super.insertString(offs, str, a);
        }
    }

    static boolean member(char item, char[] array) {
        for (char anArray : array) {
            if (anArray == item) {
                return true;
            }
        }
        return false;
    }
    //------------------------------------------------------------------------
    // Event Methods
    //------------------------------------------------------------------------

    public void addChangeListener(ChangeListener x) {
        changeListeners.add(ChangeListener.class, x);
    }

    public void removeChangeListener(ChangeListener x) {
        changeListeners.remove(ChangeListener.class, x);
    }

    protected void fireChanged() {
        range_checked = false;
        isValueValid = true;

        ChangeEvent c = new ChangeEvent(this);
        Object[] listeners = changeListeners.getListenerList();
        for (int i = listeners.length - 2; i >= 0; i -= 2) {
            if (listeners[i] == ChangeListener.class) {
                ChangeListener cl = (ChangeListener) listeners[i + 1];
                cl.stateChanged(c);
            }
        }
    }
}