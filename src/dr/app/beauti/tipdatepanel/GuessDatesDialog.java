/*
 * GuessDatesDialog.java
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

package dr.app.beauti.tipdatepanel;

import dr.app.beauti.options.DateGuesser;
import dr.app.beauti.options.STARBEASTOptions;
import dr.app.beauti.util.TextUtil;
import dr.app.gui.components.RealNumberField;
import jam.mac.Utils;
import jam.panels.OptionsPanel;


import java.io.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.*;
import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.prefs.Preferences;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class GuessDatesDialog {

    public static Preferences PREFS = Preferences.userNodeForPackage(GuessDatesDialog.class);
    public static final String DELIMIT_RADIO_KEY = "delimitRadio";
    public static final String ORDER_COMBO_KEY = "orderCombo";
    public static final String PREFIX_TEXT_KEY = "prefixText";
    public static final String REGEX_TEXT_KEY = "regexText";
    public static final String LOAD_TEXT_KEY = "loadText";
    public static final String PARSE_RADIO_KEY = "parseRadio";
    public static final String OFFSET_CHECK_KEY = "offsetCheck";
    public static final String OFFSET_TEXT_KEY = "offsetText";
    public static final String UNLESS_CHECK_KEY = "unlessCheck";
    public static final String UNLESS_TEXT_KEY = "unlessText";
    public static final String OFFSET2_TEXT_KEY = "offset2Text";
    public static final String DATE_FORMAT_TEXT_KEY = "dateFormatText";

    private JFrame frame;

    private final OptionsPanel optionPanel;

    private final JRadioButton orderRadio = new JRadioButton("Defined just by its order", true);
    private final JComboBox orderCombo = new JComboBox(new String[]{"first", "second", "third",
            "fourth", "fourth from last",
            "third from last", "second from last", "last"});
    private final JLabel orderLabel = new JLabel("Order:");
    private final JLabel prefixLabel = new JLabel("Prefix:");

    private final JRadioButton prefixRadio = new JRadioButton("Defined by a prefix and its order", false);
    private final JTextField prefixText = new JTextField(16);

    private final JRadioButton regexRadio = new JRadioButton("Defined by regular expression (REGEX)", false);
    private final JTextField regexText = new JTextField(16);

    private final JRadioButton numericalRadio = new JRadioButton("Parse as a number", true);
    private final JRadioButton calendarRadio = new JRadioButton("Parse as a calendar date", true);
    private final JRadioButton calendar2Radio = new JRadioButton("Parse calendar dates with variable precision", true);

    private final JCheckBox offsetCheck = new JCheckBox("Add the following value to each: ", false);
    private final RealNumberField offsetText = new RealNumberField();

    private final JCheckBox unlessCheck = new JCheckBox("...unless less than:", false);
    private final RealNumberField unlessText = new RealNumberField();

    private final RealNumberField offset2Text = new RealNumberField();

    private final JTextField dateFormatText = new JTextField(16);
    private String description = "Extract dates from taxon labels";

    private final int defaultDelimitRadioOption;
    private final int defaultOrderCombo;
    private final String defaultPrefixText;
    private final String defaultRegexText;
    private final int defaultParseRadioOption;
    private final boolean defaultOffsetCheckOption;
    private final String defaultOffsetText;
    private final boolean defaultUnlessCheckOption;
    private final String defaultUnlessText;
    private final String defaultOffset2Text;
    private final String defaultDateFormatText;

    public GuessDatesDialog(final JFrame frame) {
        this.frame = frame;
        defaultDelimitRadioOption = PREFS.getInt(DELIMIT_RADIO_KEY, 0);
        defaultOrderCombo = PREFS.getInt(ORDER_COMBO_KEY, 0);
        defaultPrefixText = PREFS.get(PREFIX_TEXT_KEY, "");
        defaultRegexText = PREFS.get(REGEX_TEXT_KEY, "");
        defaultParseRadioOption = PREFS.getInt(PARSE_RADIO_KEY, 0);
        defaultOffsetCheckOption = PREFS.getBoolean(OFFSET_CHECK_KEY, false);
        defaultOffsetText = PREFS.get(OFFSET_TEXT_KEY, "1900");
        defaultUnlessCheckOption = PREFS.getBoolean(UNLESS_CHECK_KEY, false);
        defaultUnlessText = PREFS.get(UNLESS_TEXT_KEY, "16");
        defaultOffset2Text = PREFS.get(OFFSET2_TEXT_KEY, "2000");
        defaultDateFormatText = PREFS.get(DATE_FORMAT_TEXT_KEY, "yyyy-MM-dd");

        optionPanel = new OptionsPanel(12, 12);
    }

    private void setupPanel(boolean parsingFromFile) {

        optionPanel.removeAll();

        if (parsingFromFile) {

        } else {
            optionPanel.addLabel("The date is given by a numerical field in the taxon label that is:");

            optionPanel.addSpanningComponent(orderRadio);
//        optionPanel.addSeparator();

            optionPanel.addSpanningComponent(prefixRadio);

            optionPanel.addComponents(orderLabel, orderCombo);
            optionPanel.addComponents(prefixLabel, prefixText);

            prefixLabel.setEnabled(false);
            prefixText.setEnabled(false);
            regexText.setEnabled(false);

            optionPanel.addComponents(regexRadio, regexText);

            optionPanel.addSeparator();
        }

        optionPanel.addSpanningComponent(numericalRadio);

        offsetText.setValue(1900);
        offsetText.setColumns(16);
        offsetText.setEnabled(false);
        optionPanel.addComponents(offsetCheck, offsetText);

        Calendar calendar = GregorianCalendar.getInstance();

        int year = calendar.get(Calendar.YEAR) - 1999;
        unlessText.setValue(year);
        unlessText.setColumns(16);
        unlessText.setEnabled(false);
        optionPanel.addComponents(unlessCheck, unlessText);

        offset2Text.setValue(2000);
        offset2Text.setColumns(16);
        offset2Text.setEnabled(false);
        final JLabel offset2Label = new JLabel("...in which case add:");
        optionPanel.addComponents(offset2Label, offset2Text);

        optionPanel.addSpanningComponent(calendarRadio);
        final JLabel dateFormatLabel = new JLabel("Date format:");
        final JButton helpButton = new JButton(Utils.isMacOSX() ? "" : "?");
        helpButton.putClientProperty("JButton.buttonType", "help");
        JPanel panel = new JPanel();
        panel.add(dateFormatText);
        panel.add(helpButton);
        panel.setOpaque(false);
        optionPanel.addComponents(dateFormatLabel, panel);
        dateFormatText.setText("yyyy-MM-dd");

        optionPanel.addSpanningComponent(calendar2Radio);

        dateFormatLabel.setEnabled(false);
        dateFormatText.setEnabled(false);

        numericalRadio.setToolTipText("Parse the date field as a decimal number");
        calendarRadio.setToolTipText("Parse the date field using a standard date format specification");
        calendar2Radio.setToolTipText("Parse the date field yyyy[-mm[-dd]] with possibly missing month or day");

        helpButton.addActionListener(new ActionListener() {
            public void actionPerformed(final ActionEvent actionEvent) {
                JScrollPane scrollPane = TextUtil.createHTMLScrollPane(
                        DATE_FORMAT_HELP, new Dimension(560,480));
                JOptionPane.showMessageDialog(frame, scrollPane,
                        "Date format help",
                        JOptionPane.PLAIN_MESSAGE);
            }
        });

        offsetCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                offsetText.setEnabled(offsetCheck.isSelected());
                unlessCheck.setEnabled(offsetCheck.isSelected());
                unlessText.setEnabled(offsetCheck.isSelected() && unlessCheck.isSelected());
                offset2Label.setEnabled(offsetCheck.isSelected() && unlessCheck.isSelected());
                offset2Text.setEnabled(offsetCheck.isSelected() && unlessCheck.isSelected());
            }
        });

        unlessCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                unlessText.setEnabled(unlessCheck.isSelected());
                offset2Label.setEnabled(unlessCheck.isSelected());
                offset2Text.setEnabled(unlessCheck.isSelected());
            }
        });

        if (!parsingFromFile) {
            ButtonGroup group = new ButtonGroup();
            group.add(orderRadio);
            group.add(prefixRadio);
            group.add(regexRadio);
            orderRadio.setSelected(true);
            ItemListener listener = new ItemListener() {
                public void itemStateChanged(ItemEvent e) {
                    orderLabel.setEnabled(!regexRadio.isSelected());
                    orderCombo.setEnabled(!regexRadio.isSelected());
                    prefixLabel.setEnabled(prefixRadio.isSelected());
                    prefixText.setEnabled(prefixRadio.isSelected());
                    regexText.setEnabled(regexRadio.isSelected());
                }
            };
            orderRadio.addItemListener(listener);
            prefixRadio.addItemListener(listener);
            regexRadio.addItemListener(listener);
        }

        ButtonGroup group = new ButtonGroup();
        group.add(numericalRadio);
        group.add(calendarRadio);
        group.add(calendar2Radio);
        ItemListener listener = new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                boolean nrs = numericalRadio.isSelected();
                boolean ocs = offsetCheck.isSelected();
                boolean ucs = unlessCheck.isSelected();
                boolean crs = calendarRadio.isSelected();
                offsetCheck.setEnabled(nrs);
                offsetText.setEnabled(nrs && ocs);
                unlessCheck.setEnabled(nrs && ocs);
                unlessText.setEnabled(nrs && ocs && ucs);
                offset2Label.setEnabled(nrs && ocs && ucs);
                offset2Text.setEnabled(nrs && ocs && ucs);
                dateFormatLabel.setEnabled(crs);
                dateFormatText.setEnabled(crs);

            }
        };
        numericalRadio.addItemListener(listener);
        calendarRadio.addItemListener(listener);
        calendar2Radio.addItemListener(listener);

        offsetCheck.setSelected(defaultOffsetCheckOption);
        offsetText.setText(defaultOffsetText);
        unlessCheck.setSelected(defaultUnlessCheckOption);
        unlessText.setText(defaultUnlessText);
        offset2Text.setText(defaultOffset2Text);

        dateFormatText.setText(defaultDateFormatText);

        orderCombo.setSelectedIndex(defaultOrderCombo);
        prefixText.setText(defaultPrefixText);
        regexText.setText(defaultRegexText);

        orderCombo.setSelectedIndex(defaultOrderCombo);
        prefixText.setText(defaultPrefixText);
        regexText.setText(defaultRegexText);

        // set from preferences defaults...
        switch (defaultDelimitRadioOption) {
            case 0: orderRadio.setSelected(true); break;
            case 1: prefixRadio.setSelected(true); break;
            case 2: regexRadio.setSelected(true); break;
            //case 3: loadRadio.setSelected(true); break; // Not allow loadRadio in prefs to avoid dialog confusion.
            default: throw new IllegalArgumentException("unknown radio option");
        }

        switch (defaultParseRadioOption) {
            case 0: numericalRadio.setSelected(true); break;
            case 1: calendarRadio.setSelected(true); break;
            case 2: calendar2Radio.setSelected(true); break;
            default: throw new IllegalArgumentException("unknown radio option");
        }

    }

    public int showDialog() {
        return showDialog(false);
    }

    public int showDialog(boolean parsingFromFile) {

        setupPanel(parsingFromFile);

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, description);
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
            setPreferencesFromDialog();
        }

        return result;
    }

    private void setPreferencesFromDialog() {
        PREFS.putInt(DELIMIT_RADIO_KEY,
                (orderRadio.isSelected() ? 0 :
                        (prefixRadio.isSelected() ? 1 :
                                (regexRadio.isSelected() ? 2 : -1))));

        PREFS.putInt(ORDER_COMBO_KEY, orderCombo.getSelectedIndex());
        PREFS.put(PREFIX_TEXT_KEY, prefixText.getText());
        PREFS.put(REGEX_TEXT_KEY, regexText.getText());

        PREFS.putInt(PARSE_RADIO_KEY,
                (numericalRadio.isSelected() ? 0 :
                        (calendarRadio.isSelected() ? 1 :
                                (calendar2Radio.isSelected() ? 2 : -1))));

        PREFS.putBoolean(OFFSET_CHECK_KEY, offsetCheck.isSelected());
        PREFS.put(OFFSET_TEXT_KEY, offsetText.getText());
        PREFS.putBoolean(UNLESS_CHECK_KEY, unlessCheck.isSelected());
        PREFS.put(UNLESS_TEXT_KEY, unlessText.getText());
        PREFS.put(OFFSET2_TEXT_KEY, offset2Text.getText());

        PREFS.put(DATE_FORMAT_TEXT_KEY, dateFormatText.getText());
    }

    public void setupGuesser(DateGuesser guesser) {
        guesser.order = orderCombo.getSelectedIndex();
        guesser.fromLast = false;
        if (guesser.order > 3) {
            guesser.fromLast = true;
            guesser.order = 8 - guesser.order - 1;
        }

        if (orderRadio.isSelected()) {
            guesser.guessType = DateGuesser.GuessType.ORDER;
        } else if (prefixRadio.isSelected()) {
            guesser.guessType = DateGuesser.GuessType.PREFIX;
            guesser.prefix = prefixText.getText();
        } else if (regexRadio.isSelected()) {
            guesser.guessType = DateGuesser.GuessType.REGEX;
            guesser.regex = regexText.getText();
        } else {
            throw new IllegalArgumentException("unknown radio button selected");
        }

        guesser.parseCalendarDatesAndPrecision = calendar2Radio.isSelected();
        guesser.parseCalendarDates = calendarRadio.isSelected();
        guesser.calendarDateFormat = dateFormatText.getText();

        guesser.offset = 0.0;
        guesser.unlessLessThan = 0.0;
        if (offsetCheck.isSelected()) {
            guesser.offset = offsetText.getValue();
            if (unlessCheck.isSelected()) {
                guesser.unlessLessThan = unlessText.getValue();
                guesser.offset2 = offset2Text.getValue();

            }
        }

    }

    public void setDescription(String description) {
        this.description = description;
    }

    private static final String DATE_FORMAT_HELP =
            "<h4>Date and Time Patterns</h4>\n" +
                    " <p>\n" +
                    " Date and time formats are specified by <em>date and time pattern</em>\n" +
                    " strings.\n" +
                    " Within date and time pattern strings, unquoted letters from\n" +
                    " <code>'A'</code> to <code>'Z'</code> and from <code>'a'</code> to\n" +
                    " <code>'z'</code> are interpreted as pattern letters representing the\n" +
                    " components of a date or time string.\n" +
                    " Text can be quoted using single quotes (<code>'</code>) to avoid\n" +
                    " interpretation.\n" +
                    " <code>\"''\"</code> represents a single quote.\n" +
                    " All other characters are not interpreted; they're simply copied into the\n" +
                    " output string during formatting or matched against the input string\n" +
                    " during parsing.\n" +
                    " <p>\n" +
                    " The following pattern letters are defined (all other characters from\n" +
                    " <code>'A'</code> to <code>'Z'</code> and from <code>'a'</code> to\n" +
                    " <code>'z'</code> are reserved):\n" +
                    " <blockquote>\n" +
                    " <table border=0 cellspacing=3 cellpadding=0 summary=\"Chart shows pattern letters, date/time component, presentation, and examples.\">\n" +
                    "     <tr bgcolor=\"#ccccff\">\n" +
                    "         <th align=left>Letter\n" +
                    "         <th align=left>Date or Time Component\n" +
                    "         <th align=left>Presentation\n" +
                    "         <th align=left>Examples\n" +
                    "     <tr>\n" +
                    "         <td><code>G</code>\n" +
                    "         <td>Era designator\n" +
                    "         <td><a href=\"#text\">Text</a>\n" +
                    "         <td><code>AD</code>\n" +
                    "     <tr bgcolor=\"#eeeeff\">\n" +
                    "         <td><code>y</code>\n" +
                    "         <td>Year\n" +
                    "         <td><a href=\"#year\">Year</a>\n" +
                    "         <td><code>1996</code>; <code>96</code>\n" +
                    "     <tr>\n" +
                    "         <td><code>M</code>\n" +
                    "         <td>Month in year\n" +
                    "         <td><a href=\"#month\">Month</a>\n" +
                    "         <td><code>July</code>; <code>Jul</code>; <code>07</code>\n" +
                    "     <tr bgcolor=\"#eeeeff\">\n" +
                    "         <td><code>w</code>\n" +
                    "         <td>Week in year\n" +
                    "         <td><a href=\"#number\">Number</a>\n" +
                    "         <td><code>27</code>\n" +
                    "     <tr>\n" +
                    "         <td><code>W</code>\n" +
                    "         <td>Week in month\n" +
                    "         <td><a href=\"#number\">Number</a>\n" +
                    "         <td><code>2</code>\n" +
                    "     <tr bgcolor=\"#eeeeff\">\n" +
                    "         <td><code>D</code>\n" +
                    "         <td>Day in year\n" +
                    "         <td><a href=\"#number\">Number</a>\n" +
                    "         <td><code>189</code>\n" +
                    "     <tr>\n" +
                    "         <td><code>d</code>\n" +
                    "         <td>Day in month\n" +
                    "         <td><a href=\"#number\">Number</a>\n" +
                    "         <td><code>10</code>\n" +
                    "     <tr bgcolor=\"#eeeeff\">\n" +
                    "         <td><code>F</code>\n" +
                    "         <td>Day of week in month\n" +
                    "         <td><a href=\"#number\">Number</a>\n" +
                    "         <td><code>2</code>\n" +
                    "     <tr>\n" +
                    "         <td><code>E</code>\n" +
                    "         <td>Day in week\n" +
                    "         <td><a href=\"#text\">Text</a>\n" +
                    "         <td><code>Tuesday</code>; <code>Tue</code>\n" +
                    "     <tr bgcolor=\"#eeeeff\">\n" +
                    "         <td><code>a</code>\n" +
                    "         <td>Am/pm marker\n" +
                    "         <td><a href=\"#text\">Text</a>\n" +
                    "         <td><code>PM</code>\n" +
                    "     <tr>\n" +
                    "         <td><code>H</code>\n" +
                    "         <td>Hour in day (0-23)\n" +
                    "         <td><a href=\"#number\">Number</a>\n" +
                    "         <td><code>0</code>\n" +
                    "     <tr bgcolor=\"#eeeeff\">\n" +
                    "         <td><code>k</code>\n" +
                    "         <td>Hour in day (1-24)\n" +
                    "         <td><a href=\"#number\">Number</a>\n" +
                    "         <td><code>24</code>\n" +
                    "     <tr>\n" +
                    "         <td><code>K</code>\n" +
                    "         <td>Hour in am/pm (0-11)\n" +
                    "         <td><a href=\"#number\">Number</a>\n" +
                    "         <td><code>0</code>\n" +
                    "     <tr bgcolor=\"#eeeeff\">\n" +
                    "         <td><code>h</code>\n" +
                    "         <td>Hour in am/pm (1-12)\n" +
                    "         <td><a href=\"#number\">Number</a>\n" +
                    "         <td><code>12</code>\n" +
                    "     <tr>\n" +
                    "         <td><code>m</code>\n" +
                    "         <td>Minute in hour\n" +
                    "         <td><a href=\"#number\">Number</a>\n" +
                    "         <td><code>30</code>\n" +
                    "     <tr bgcolor=\"#eeeeff\">\n" +
                    "         <td><code>s</code>\n" +
                    "         <td>Second in minute\n" +
                    "         <td><a href=\"#number\">Number</a>\n" +
                    "         <td><code>55</code>\n" +
                    "     <tr>\n" +
                    "         <td><code>S</code>\n" +
                    "         <td>Millisecond\n" +
                    "         <td><a href=\"#number\">Number</a>\n" +
                    "         <td><code>978</code>\n" +
                    "     <tr bgcolor=\"#eeeeff\">\n" +
                    "         <td><code>z</code>\n" +
                    "         <td>Time zone\n" +
                    "         <td><a href=\"#timezone\">General time zone</a>\n" +
                    "         <td><code>Pacific Standard Time</code>; <code>PST</code>; <code>GMT-08:00</code>\n" +
                    "     <tr>\n" +
                    "         <td><code>Z</code>\n" +
                    "         <td>Time zone\n" +
                    "         <td><a href=\"#rfc822timezone\">RFC 822 time zone</a>\n" +
                    "         <td><code>-0800</code>\n" +
                    " </table>\n" +
                    " </blockquote>\n" +
                    " Pattern letters are usually repeated, as their number determines the\n" +
                    " exact presentation:\n" +
                    " <ul>\n" +
                    " <li><strong><a name=\"text\">Text:</a></strong>\n" +
                    "     For formatting, if the number of pattern letters is 4 or more,\n" +
                    "     the full form is used; otherwise a short or abbreviated form\n" +
                    "     is used if available.\n" +
                    "     For parsing, both forms are accepted, independent of the number\n" +
                    "     of pattern letters.\n" +
                    " <li><strong><a name=\"number\">Number:</a></strong>\n" +
                    "     For formatting, the number of pattern letters is the minimum\n" +
                    "     number of digits, and shorter numbers are zero-padded to this amount.\n" +
                    "     For parsing, the number of pattern letters is ignored unless\n" +
                    "     it's needed to separate two adjacent fields.\n" +
                    " <li><strong><a name=\"year\">Year:</a></strong>\n" +
                    "     If the formatter's <A HREF=\"../../java/text/DateFormat.html#getCalendar()\"><CODE>Calendar</CODE></A> is the Gregorian\n" +
                    "     calendar, the following rules are applied.<br>\n" +
                    "     <ul>\n" +
                    "     <li>For formatting, if the number of pattern letters is 2, the year\n" +
                    "         is truncated to 2 digits; otherwise it is interpreted as a\n" +
                    "         <a href=\"#number\">number</a>.\n" +
                    "     <li>For parsing, if the number of pattern letters is more than 2,\n" +
                    "         the year is interpreted literally, regardless of the number of\n" +
                    "         digits. So using the pattern \"MM/dd/yyyy\", \"01/11/12\" parses to\n" +
                    "         Jan 11, 12 A.D.\n" +
                    "     <li>For parsing with the abbreviated year pattern (\"y\" or \"yy\"),\n" +
                    "         <code>SimpleDateFormat</code> must interpret the abbreviated year\n" +
                    "         relative to some century.  It does this by adjusting dates to be\n" +
                    "         within 80 years before and 20 years after the time the <code>SimpleDateFormat</code>\n" +
                    "         instance is created. For example, using a pattern of \"MM/dd/yy\" and a\n" +
                    "         <code>SimpleDateFormat</code> instance created on Jan 1, 1997,  the string\n" +
                    "         \"01/11/12\" would be interpreted as Jan 11, 2012 while the string \"05/04/64\"\n" +
                    "         would be interpreted as May 4, 1964.\n" +
                    "         During parsing, only strings consisting of exactly two digits, as defined by\n" +
                    "         <A HREF=\"../../java/lang/Character.html#isDigit(char)\"><CODE>Character.isDigit(char)</CODE></A>, will be parsed into the default century.\n" +
                    "         Any other numeric string, such as a one digit string, a three or more digit\n" +
                    "         string, or a two digit string that isn't all digits (for example, \"-1\"), is\n" +
                    "         interpreted literally.  So \"01/02/3\" or \"01/02/003\" are parsed, using the\n" +
                    "         same pattern, as Jan 2, 3 AD.  Likewise, \"01/02/-3\" is parsed as Jan 2, 4 BC.\n" +
                    "     </ul>\n" +
                    "     Otherwise, calendar system specific forms are applied.\n" +
                    "     For both formatting and parsing, if the number of pattern\n" +
                    "     letters is 4 or more, a calendar specific <A HREF=\"../../java/util/Calendar.html#LONG\">long form</A> is used. Otherwise, a calendar\n" +
                    "     specific <A HREF=\"../../java/util/Calendar.html#SHORT\">short or abbreviated form</A>\n" +
                    "     is used.\n" +
                    " <li><strong><a name=\"month\">Month:</a></strong>\n" +
                    "     If the number of pattern letters is 3 or more, the month is\n" +
                    "     interpreted as <a href=\"#text\">text</a>; otherwise,\n" +
                    "     it is interpreted as a <a href=\"#number\">number</a>.\n" +
                    " </ul>";

}