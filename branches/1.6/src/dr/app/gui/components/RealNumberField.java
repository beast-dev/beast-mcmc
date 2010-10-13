/*
 * RealNumberField.java
 *
 * Copyright (c) 2009 JAM Development Team
 *
 * This package is distributed under the Lesser Gnu Public Licence (LGPL)
 *
 */

package dr.app.gui.components;

import dr.util.NumberFormatter;
import dr.app.beauti.util.NumberUtil;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.PlainDocument;
import java.awt.*;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.text.DecimalFormat;

public class RealNumberField extends JTextField implements FocusListener, DocumentListener {

    public static String NaN = "NaN";
    public static String POSITIVE_INFINITY = "+INF";
    public static String NEGATIVE_INFINITY = "-INF";
    public static String MAX_VALUE = "MAX_VALUE";
    public static String MIN_VALUE = "MIN_VALUE";

    protected static char MINUS = '-';
    protected static char PERIOD = '.';
    protected EventListenerList changeListeners = new EventListenerList();
    protected double min;
    protected double max;
    protected boolean range_check = false;
    protected boolean range_checked = false;

    public RealNumberField() {
        super();
    }

    public RealNumberField(double min, double max) {
        this();
        this.min = min;
        this.max = max;
        range_check = true;
        this.addFocusListener(this);
    }

    public void focusGained(FocusEvent evt) {
    }

    public void focusLost(FocusEvent evt) {
        if (range_check && !range_checked) {
            range_checked = true;
            try {
                double value = getValue();
                if (value < min || value > max) {
                    errorMsg();
                }
            } catch (NumberFormatException e) {
                errorMsg();
            }
        }
    }

    public void setText(double value) {
        if (value == Double.NaN) {
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
            setText(NumberUtil.formatDecimal(value, 10, 6));
        }
    }

    public void setText(Integer obj) {
        setText(obj.toString()); // where used?
    }

    public void setText(Long obj) {
        setText(obj.toString()); // where used?
    }

    protected void errorMsg() {
        JOptionPane.showMessageDialog(this,
                "Illegal entry\nValue must be between " + min + " and " +
                max + " inclusive", "Error", JOptionPane.ERROR_MESSAGE);
    }

    public void setRange(double min, double max) {
        this.min = min;
        this.max = max;
        range_check = true;
    }

    public void setValue(double value) {
        if (range_check) {
            if (value < min || value > max) {
                errorMsg();
                return;
            }
        }
        setText(value);
    }

    public Double getValue() {
        try {
            if (getText().trim() == "" || getText().trim() == null) {
//                System.out.println("null");
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
                return new Double(getText());
            }
        } catch (NumberFormatException e) {
            JOptionPane.showMessageDialog(this, "Unable to parse number correctly",
                        "Number Format Exception",
                        JOptionPane.ERROR_MESSAGE);
            return null;
        }
    }

    protected Document createDefaultModel() {
        Document doc = new RealNumberField.RealNumberFieldDocument();
        doc.addDocumentListener(this);
        return doc;
    }

    public void insertUpdate(DocumentEvent e) {
        range_checked = false;
        fireChanged();
    }

    public void removeUpdate(DocumentEvent e) {
        range_checked = false;
        fireChanged();
    }

    public void changedUpdate(DocumentEvent e) {
        range_checked = false;
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
            boolean exponent_found =  false;
            int exponent_index = -1;
            boolean exponent_sign_found =  false;

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