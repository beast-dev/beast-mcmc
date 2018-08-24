/*
 * DateGuesser.java
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

package dr.app.beauti.options;

import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;

import java.io.*;
import java.text.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Andrew Rambaut
 * @author Tommy Lam
 */
public class DateGuesser implements Serializable {
    private static final long serialVersionUID = -9106689400887615213L;

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
    public File loadFile;
    public HashMap<String, String> load;
    public double offset = 0.0;
    public double unlessLessThan = 0.0;
    public double offset2 = 0.0;
    public boolean parseCalendarDates = false;
    public boolean parseCalendarDatesAndPrecision = false;
    public String calendarDateFormat = "yyyy-MM-dd";

    private DateFormat dateFormat;

    public void guessDates(TaxonList taxonList) {
        // To avoid duplicating code, add all the taxa into a list and
        // pass it to guessDates(List<Taxon> taxonList)
        List<Taxon> taxa = new ArrayList<Taxon>();
        for (Taxon taxon : taxonList) {
            taxa.add(taxon);
        }

        guessDates(taxa);
    }

    public void guessDates(TaxonList taxonList, Map<Taxon, String> taxonDateMap) {
        // To avoid duplicating code, add all the taxa into a list and
        // pass it to guessDates(List<Taxon> taxonList)
        List<Taxon> taxa = new ArrayList<Taxon>();
        for (Taxon taxon : taxonList) {
            taxa.add(taxon);
        }

        guessDates(taxa, taxonDateMap);
    }

    public void guessDates(List<Taxon> taxonList) {
        guessDates(taxonList, null);
    }

    public void guessDates(List<Taxon> taxonList, Map<Taxon, String> taxonDateMap) {

        dateFormat = new SimpleDateFormat(calendarDateFormat);
        dateFormat.setTimeZone(TimeZone.getTimeZone("GMT"));

        for (int i = 0; i < taxonList.size(); i++) {
            Taxon taxon = taxonList.get(i);

            // Allocates a Date object and initializes it to represent the specified number of milliseconds since the
            // standard base time known as "the epoch", namely January 1, 1970, 00:00:00 GMT
            java.util.Date origin = new java.util.Date(0);

            double[] values = new double[2];

            try {

                if (taxonDateMap != null) {
                    String dateString = taxonDateMap.get(taxon);
                    if (dateString == null) {
                        continue;
                    }

                    parseDate(taxon.getId(), dateString, values);

                } else {
                    switch (guessType) {
                        case ORDER:
                            guessDateFromOrder(taxonList.get(i).getId(), order, fromLast, values);
                            break;
                        case PREFIX:
                            guessDateFromPrefix(taxonList.get(i).getId(), prefix, order, fromLast, values);
                            break;
                        case REGEX:
                            guessDateFromRegex(taxonList.get(i).getId(), regex, values);
                            break;
                        default:
                            throw new IllegalArgumentException("unknown GuessType");
                    }
                }

            } catch (GuessDatesException gfe) {
                // @todo catch errors and give to user
            }

            double d = values[0];

            if (!parseCalendarDates && !parseCalendarDatesAndPrecision) {
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
            }

            // @todo if any taxa aren't set then return warning

            Date date = Date.createTimeSinceOrigin(d, Units.Type.YEARS, origin);
            date.setUncertainty(values[1]);
            taxon.setAttribute("date", date);
        }
    }

    public Date parseDate(String value) throws GuessDatesException {
        double[] values = new double[2];
        parseDate("", value, values);

        // Allocates a Date object and initializes it to represent the specified number of milliseconds since the
        // standard base time known as "the epoch", namely January 1, 1970, 00:00:00 GMT
        java.util.Date origin = new java.util.Date(0);

        return Date.createTimeSinceOrigin(values[0], Units.Type.YEARS, origin);

    }

    private void guessDateFromOrder(String label, int order, boolean fromLast, double[] values) throws GuessDatesException {

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

        parseDate(label, field, values);
    }

    private static final String REGEX_CHARACTERS = "|[].*()-^$";

    private void guessDateFromPrefix(String label, String prefix, int order, boolean fromLast, double[] values) throws GuessDatesException {

        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < prefix.length(); i++) {
            if (REGEX_CHARACTERS.contains("" + prefix.charAt(i))) {
                sb.append("\\").append(prefix.charAt(i));
            } else {
                sb.append(prefix.charAt(i));
            }
        }

        String[] fields = label.split(sb.toString());
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

        parseDate(label, fields[index], values);
    }

    private void guessDateFromPrefix(String label, String prefix, double[] values) throws GuessDatesException {

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

        parseDate(label, field, values);
    }

    private void guessDateFromRegex(String label, String regex, double[] values) throws GuessDatesException {

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

        parseDate(label, matcher.group(0), values);
    }


    private void parseDateFromValue(String label, HashMap<String, String> myload, double[] values) throws GuessDatesException {
        String dateStr = "";
        if (myload.containsKey(label)) {
            dateStr = (String)(myload.get(label));
        } else {
            throw new GuessDatesException("The imported table doesn't contain the taxon label, " + label);
        }

        parseDate(label, dateStr, values);
    }


    private DateFormat dateFormat1 = null;
    private DateFormat dateFormat2 = null;
    private DateFormat dateFormat3 = null;

    private void parseDate(String label, String value, double[] values) throws GuessDatesException {
        double d;
        double p = 0.0;

        if (dateFormat1 == null) {
            // set the timezones to GMT so they match the origin date...
            dateFormat1 = new SimpleDateFormat("yyyy-MM-dd");
            dateFormat1.setTimeZone(TimeZone.getTimeZone("GMT"));
            dateFormat2 = new SimpleDateFormat("yyyy-MM");
            dateFormat2.setTimeZone(TimeZone.getTimeZone("GMT"));
            dateFormat3 = new SimpleDateFormat("yyyy");
            dateFormat3.setTimeZone(TimeZone.getTimeZone("GMT"));
        }

        if (parseCalendarDatesAndPrecision) {
            try {
                Date date = new Date(dateFormat1.parse(value));

                d = date.getTimeValue();
                p = 0.0;

            } catch (ParseException pe) {
                try {
                    Date date = new Date(dateFormat2.parse(value));

                    d = date.getTimeValue();
                    p = 1.0 / 12.0;

                } catch (ParseException pe2) {
                    try {
                        Date date = new Date(dateFormat3.parse(value));

                        d = date.getTimeValue();
                        p = 1.0;

                    } catch (ParseException pe3) {
                        throw new GuessDatesException("Badly formatted date for taxon, " + label);
                    }
                }
            }

        } else if (parseCalendarDates) {
            try {

                Date date = new Date(dateFormat.parse(value));

                d = date.getTimeValue();
            } catch (ParseException pe) {
                throw new GuessDatesException("Badly formatted date for taxon, " + label);
            }

        } else {
            try {
                d = Double.parseDouble(value);
            } catch (NumberFormatException nfe) {
                throw new GuessDatesException("Badly formatted date for taxon, " + label);
            }
        }

        values[0] = d;
        values[1] = p;
    }
}
