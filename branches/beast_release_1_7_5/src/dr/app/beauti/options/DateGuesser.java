/*
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.options;

import dr.app.beauti.tipdatepanel.GuessDatesException;
import dr.evolution.util.Date;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;

import java.io.Serializable;
import java.text.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Andrew Rambaut
 */
public class DateGuesser implements Serializable {
    public enum GuessType {
        ORDER,
        PREFIX,
        REGEX
    }

    public boolean guessDates = false;
    public GuessType guessType = GuessType.ORDER;
    public boolean fromLast = false;
    public int order = 0;
    public String prefix;
    public String regex;
    public double offset = 0.0;
    public double unlessLessThan = 0.0;
    public double offset2 = 0.0;
    public boolean parseCalendarDates = false;
    public String calendarDateFormat = "dd/MM/yyyy";

    private DateFormat dateFormat = new SimpleDateFormat(calendarDateFormat);

    public void guessDates(TaxonList taxonList) {

        dateFormat = new SimpleDateFormat(calendarDateFormat);

        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            java.util.Date origin = new java.util.Date(0);

            double d = 0.0;

            try {
                switch (guessType) {
                    case ORDER:
                        d = guessDateFromOrder(taxonList.getTaxonId(i), order, fromLast);
                        break;
                    case PREFIX:
                        d = guessDateFromPrefix(taxonList.getTaxonId(i), prefix, order, fromLast);
                        break;
                    case REGEX:
                        d = guessDateFromRegex(taxonList.getTaxonId(i), regex);
                        break;
                    default:
                        throw new IllegalArgumentException("unknown GuessType");
                }

            } catch (GuessDatesException gfe) {
                //
            }

            if (offset > 0) {
                if (unlessLessThan > 0) {
                    if (d < unlessLessThan) {
                        d += offset2;
                    } else {
                        d += offset;
                    }
                } else {
                    d += offset;
                }
            }

            Date date = Date.createTimeSinceOrigin(d, Units.Type.YEARS, origin);
            taxonList.getTaxon(i).setAttribute("date", date);
        }
    }

    private double guessDateFromOrder(String label, int order, boolean fromLast) throws GuessDatesException {

        String field;

        if (fromLast) {
            int count = 0;
            int i = label.length() - 1;

            char c = label.charAt(i);

            do {
                // first find a part of a number
                while (!Character.isDigit(c) && c != '.') {
                    i--;
                    if (i < 0) break;
                    c = label.charAt(i);
                }

                if (i < 0) throw new GuessDatesException("Missing number field in taxon label, " + label);

                int j = i + 1;

                // now find the beginning of the number
                while (Character.isDigit(c) || c == '.') {
                    i--;
                    if (i < 0) break;
                    c = label.charAt(i);
                }

                field = label.substring(i + 1, j);

                count++;

            } while (count <= order);

        } else {
            int count = 0;
            int i = 0;

            char c = label.charAt(i);

            do {
                // first find a part of a number
                while (!Character.isDigit(c)) {
                    i++;
                    if (i == label.length()) break;
                    c = label.charAt(i);
                }
                int j = i;

                if (i == label.length()) throw new GuessDatesException("Missing number field in taxon label, " + label);

                // now find the beginning of the number
                while (Character.isDigit(c) || c == '.') {
                    i++;
                    if (i == label.length()) break;
                    c = label.charAt(i);
                }

                field = label.substring(j, i);

                count++;

            } while (count <= order);
        }

        return parseDate(label, field);
    }

    private double guessDateFromPrefix(String label, String prefix, int order, boolean fromLast) throws GuessDatesException {

        String[] fields = label.split(prefix);
        int index;
        if (fromLast) {
            index = fields.length - order - 1;
        } else {
            index = order;
        }
        if (index < 0) {
            index = 0;
        }
        if (index >= fields.length) {
            index = fields.length - 1;
        }

        return parseDate(label, fields[index]);
    }

    private double guessDateFromPrefix(String label, String prefix) throws GuessDatesException {

        int i = label.indexOf(prefix);

        if (i == -1) throw new GuessDatesException("Missing prefix in taxon label, " + label);

        i += prefix.length();
        int j = i;

        // now find the beginning of the number
        char c = label.charAt(i);
        while (i < label.length() - 1 && (Character.isDigit(c) || c == '.')) {
            i++;
            c = label.charAt(i);
        }

        if (i == j) throw new GuessDatesException("Missing field after prefix in taxon label, " + label);

        String field = label.substring(j, i + 1);

        return parseDate(label, field);
    }

    private double guessDateFromRegex(String label, String regex) throws GuessDatesException {
        double d;

        if (!regex.contains("(")) {
            // if user hasn't specified a replace element, assume the whole regex should match
            regex = "(" + regex + ")";
        }

        Pattern pattern = Pattern.compile(regex);
        Matcher matcher = pattern.matcher(label);
        if (!matcher.find()) {
            throw new GuessDatesException("Regular expression doesn't find a match in taxon label, " + label);
        }

        if (matcher.groupCount() < 1) {
            throw new GuessDatesException("Date group not defined in regular expression");
        }

        return parseDate(label, matcher.group(0));
    }

    private double parseDate(String label, String value) throws GuessDatesException {
        double d;
        try {
            d = Double.parseDouble(value);
        } catch (NumberFormatException nfe) {
            try {
                Date date = new Date(dateFormat.parse(value));

                d = date.getTimeValue();
            } catch (ParseException pe) {
                throw new GuessDatesException("Badly formated date in taxon label, " + label);
            }
        }

        return d;
    }
}
