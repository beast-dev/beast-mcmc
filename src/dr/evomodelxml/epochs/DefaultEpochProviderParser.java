/*
 * DefaultEpochProviderParser.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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

package dr.evomodelxml.epochs;

import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.TimeScale;
import dr.evolution.util.Units;
import dr.evomodel.epochs.DefaultEpochProvider;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.xml.*;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TimeZone;
import java.util.logging.Logger;

/**
 * Epochs with fixed transition times, given either as a number of epochs of the same duration:
 * <pre>
 *     &lt;epochs id="epochs" epochCount="12" duration="1" units="months"/&gt;
 * </pre>
 * or as the transition times themselves:
 * <pre>
 *     &lt;epochs id="epochs" transitionTimes="1 2 4 8 26" units="weeks"/&gt;
 * </pre>
 * The durations or times are in "units" (days, weeks, months or years; years if it is not given) and are
 * converted to years, the units of time of the tree.
 *
 * Height zero is the date of the most recent tip: the most recent of the dates of all the taxa read so far.
 * The transition times are normally counted back from it. "origin" gives a date, on the scale of the dates of
 * the taxa, from which they are counted back instead:
 * <pre>
 *     &lt;epochs id="epochs" epochCount="12" duration="1" units="months" origin="2020.0"/&gt;
 * </pre>
 * puts the transitions one month, two months ... before the start of 2020, whatever the date of the most
 * recent tip. The origin can also be a calendar date, "2020-01-01", when the dates of the taxa are in years.
 *
 * The transitions can be given as calendar dates too, in any order, again when the dates of the taxa are in
 * years:
 * <pre>
 *     &lt;epochs id="epochs" transitionTimes="2020-03-23 2020-07-04 2020-11-05"/&gt;
 * </pre>
 * A calendar date is the start of that day and is converted to a decimal year as BEAUti converts the dates
 * of tips.
 *
 * The epochs can be given as date elements inside the element, with no epochCount, duration, units
 * or transitionTimes:
 * <pre>
 *     &lt;epochs id="epochs"&gt;
 *         &lt;date value="2021-01-01"/&gt;
 *         &lt;date value="2020-11-05"/&gt;
 *         &lt;date value="2020-07-04"/&gt;
 *         &lt;date value="2020-03-23"/&gt;
 *     &lt;/epochs&gt;
 * </pre>
 * The most recent of the dates is the origin and the others are the transitions. With an origin attribute
 * all the dates are transitions. Each transition is at the height of its date, the time before the most
 * recent tip; if no taxon has a date, the origin is taken as the time of the most recent tip.
 *
 *
 * A transitionTimes element holds a parameter for the transition times. If the parameter has values, they
 * are the transition times (heights in years) and nothing else is needed:
 * <pre>
 *     &lt;epochs id="epochs"&gt;
 *         &lt;transitionTimes&gt;
 *             &lt;parameter id="epochs.times" value="0.5 1.0 2.0"/&gt;
 *         &lt;/transitionTimes&gt;
 *     &lt;/epochs&gt;
 * </pre>
 * If it has no values, the epochs are given in one of the ways above and the parameter is set to the
 * transition times, so that operators can be applied to them:
 * <pre>
 *     &lt;epochs id="epochs" epochCount="12" duration="1" units="months"&gt;
 *         &lt;transitionTimes&gt;
 *             &lt;parameter id="epochs.times"/&gt;
 *         &lt;/transitionTimes&gt;
 *     &lt;/epochs&gt;
 * </pre>
 *
 * @author Andrew Rambaut
 */
public class DefaultEpochProviderParser extends AbstractXMLObjectParser {

    public static final String EPOCHS = DefaultEpochProvider.EPOCHS;
    public static final String EPOCH_COUNT = "epochCount";
    public static final String DURATION = "duration";
    public static final String TRANSITION_TIMES = DefaultEpochProvider.TRANSITION_TIMES;
    public static final String UNITS = "units";
    public static final String ORIGIN = "origin";

    private static final String[] TIME_UNIT_NAMES = {"days", "weeks", "months", "years"};
    private static final Units.Type[] TIME_UNIT_TYPES =
            {Units.Type.DAYS, Units.Type.WEEKS, Units.Type.MONTHS, Units.Type.YEARS};

    // the time of a tree is in years
    private static final Units.Type TREE_UNITS = Units.Type.YEARS;

    public String getParserName() {
        return EPOCHS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final StringBuilder notes = new StringBuilder();

        final List<Date> dates = new ArrayList<Date>();
        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof Date) {
                dates.add((Date) xo.getChild(i));
            }
        }

        // a parameter for the transition times, which either has them as its values or is given them
        Parameter parameter = null;
        boolean parameterHasValues = false;
        if (xo.hasChildNamed(TRANSITION_TIMES)) {
            final XMLObject cxo = xo.getChild(TRANSITION_TIMES);
            parameter = (Parameter) cxo.getChild(Parameter.class);
            if (parameter == null) {
                throw new XMLParseException("The " + TRANSITION_TIMES + " element of the " + EPOCHS + " element" +
                        idText(xo) + " should contain a parameter");
            }
            parameterHasValues = hasValues(cxo);
        }

        final DefaultEpochProvider epochs;
        final double[] transitionTimes;
        try {
            if (parameterHasValues) {
                // the values of the parameter are the transition times, as heights
                for (String attribute : new String[]{EPOCH_COUNT, DURATION, TRANSITION_TIMES, UNITS, ORIGIN}) {
                    if (xo.hasAttribute(attribute)) {
                        throw new XMLParseException("The " + EPOCHS + " element" + idText(xo) + " gives its " +
                                "transition times as the values of a parameter, so it should not have the " +
                                attribute + " attribute");
                    }
                }
                if (!dates.isEmpty()) {
                    throw new XMLParseException("The " + EPOCHS + " element" + idText(xo) + " gives its transition " +
                            "times as the values of a parameter, so it should not have dates as well");
                }
                transitionTimes = parameter.getParameterValues();
                epochs = new DefaultEpochProvider(parameter);
            } else {
                transitionTimes = dates.isEmpty() ?
                        getTransitionTimesFromAttributes(xo, notes) :
                        getTransitionTimesFromDates(xo, dates, notes);
                epochs = parameter != null ?
                        new DefaultEpochProvider(parameter, transitionTimes) :
                        new DefaultEpochProvider(transitionTimes);
            }
        } catch (IllegalArgumentException iae) {
            throw new XMLParseException("In the " + EPOCHS + " element" + idText(xo) + ": " + iae.getMessage());
        }

        final Logger logger = Logger.getLogger("dr.evomodel");
        logger.info("\nEpochs" + idText(xo) + ": " + epochs.getEpochCount() +
                (epochs.getEpochCount() == 1 ? " epoch" : " epochs"));
        if (transitionTimes.length > 0) {
            logger.info(notes + "  transition times (years before the most recent tip): " +
                    Arrays.toString(transitionTimes));
        }

        return epochs;
    }

    /**
     * @return whether the parameter in the element was given values (by itself, or where it is defined if it
     * is a reference): if not, it is there to be given them
     */
    private static boolean hasValues(XMLObject cxo) throws XMLParseException {
        for (int i = 0; i < cxo.getChildCount(); i++) {
            if (cxo.getChild(i) instanceof Parameter) {
                Object raw = cxo.getRawChild(i);
                if (raw instanceof Reference) {
                    raw = ((Reference) raw).getReferenceObject();
                }
                return raw instanceof XMLObject && ((XMLObject) raw).hasAttribute(ParameterParser.VALUE);
            }
        }
        return false;
    }

    /**
     * The transition times from the epochCount and duration attributes, or the transitionTimes attribute.
     */
    private static double[] getTransitionTimesFromAttributes(XMLObject xo, StringBuilder notes) throws XMLParseException {

        final boolean hasEqualDurations = xo.hasAttribute(EPOCH_COUNT) || xo.hasAttribute(DURATION);
        final boolean hasTransitionTimes = xo.hasAttribute(TRANSITION_TIMES);
        if (hasEqualDurations == hasTransitionTimes) {
            throw new XMLParseException("The " + EPOCHS + " element" + idText(xo) + " needs either " +
                    EPOCH_COUNT + " and " + DURATION + ", or " + TRANSITION_TIMES + " (and not both)");
        }

        // the transition times are either all numbers or all calendar dates
        final String[] calendarDates = hasTransitionTimes ? getCalendarDates(xo) : null;
        final boolean hasCalendarDates = calendarDates != null;

        double[] transitionTimes;
        try {
            if (hasEqualDurations) {
                if (!xo.hasAttribute(EPOCH_COUNT) || !xo.hasAttribute(DURATION)) {
                    throw new XMLParseException("The " + EPOCHS + " element" + idText(xo) + " needs both " +
                            EPOCH_COUNT + " and " + DURATION);
                }
                transitionTimes = DefaultEpochProvider.equalDurationTransitionTimes(
                        xo.getIntegerAttribute(EPOCH_COUNT), xo.getDoubleAttribute(DURATION));
            } else if (!hasCalendarDates) {
                transitionTimes = xo.getDoubleArrayAttribute(TRANSITION_TIMES);
            } else {
                transitionTimes = null;
            }
        } catch (IllegalArgumentException iae) {
            throw new XMLParseException("In the " + EPOCHS + " element" + idText(xo) + ": " + iae.getMessage());
        }

        if (hasCalendarDates) {
            if (xo.hasAttribute(UNITS) || xo.hasAttribute(ORIGIN)) {
                throw new XMLParseException("The " + EPOCHS + " element" + idText(xo) + " gives its " +
                        TRANSITION_TIMES + " as calendar dates, which need no " + UNITS + " or " + ORIGIN);
            }
            // a date is the time before the most recent tip that it falls
            final Date mostRecent = getMostRecentCalendarDate(xo, "calendar dates for " + TRANSITION_TIMES);
            transitionTimes = new double[calendarDates.length];
            for (int i = 0; i < calendarDates.length; i++) {
                transitionTimes[i] = mostRecent.getTimeValue() - parseCalendarDate(calendarDates[i], xo);
            }
            Arrays.sort(transitionTimes);
            if (transitionTimes.length > 0 && !(transitionTimes[0] > 0.0)) {
                throw new XMLParseException("In the " + EPOCHS + " element" + idText(xo) + ", the latest of the " +
                        TRANSITION_TIMES + " is not before the most recent tip (" + mostRecent.getTimeValue() + ")");
            }
            String[] sortedDates = calendarDates.clone();
            Arrays.sort(sortedDates, java.util.Collections.reverseOrder());
            notes.append("  transitions at ").append(String.join(", ", sortedDates))
                    .append("; the most recent tip is at ").append(mostRecent.getTimeValue()).append("\n");
        } else {
            if (xo.hasAttribute(UNITS)) {
                transitionTimes = convert(transitionTimes, parseTimeUnits(xo.getStringAttribute(UNITS), xo), TREE_UNITS);
            }

            if (xo.hasAttribute(ORIGIN)) {
                // the transition times are counted back from the origin, not from the most recent tip
                final String originAttribute = xo.getStringAttribute(ORIGIN).trim();
                final Date mostRecent;
                final double origin;
                if (isCalendarDate(originAttribute)) {
                    mostRecent = getMostRecentCalendarDate(xo, "a calendar date for the " + ORIGIN);
                    origin = parseCalendarDate(originAttribute, xo);
                } else {
                    mostRecent = getMostRecentDate(xo);
                    try {
                        origin = Double.parseDouble(originAttribute);
                    } catch (NumberFormatException nfe) {
                        throw new XMLParseException("In the " + EPOCHS + " element" + idText(xo) + ", " + ORIGIN +
                                "=\"" + originAttribute + "\" is neither a number nor a calendar date (yyyy-mm-dd)");
                    }
                }
                // dates run forwards in time (2020.5) or backwards (an age of 3.5)
                final double timeBeforeMostRecent = mostRecent.isBackwards() ?
                        origin - mostRecent.getTimeValue() : mostRecent.getTimeValue() - origin;
                final double originHeight =
                        TimeScale.convertTimeUnits(timeBeforeMostRecent, mostRecent.getUnits(), TREE_UNITS);
                for (int i = 0; i < transitionTimes.length; i++) {
                    transitionTimes[i] += originHeight;
                }
                if (transitionTimes.length > 0 && !(transitionTimes[0] > 0.0)) {
                    throw new XMLParseException("In the " + EPOCHS + " element" + idText(xo) + ", the first " +
                            "transition between epochs is not before the most recent tip (" +
                            mostRecent.getTimeValue() + "): with an " + ORIGIN + " of " + originAttribute +
                            " it is " + (-transitionTimes[0]) + " years after it.");
                }
                notes.append("  counted back from ").append(originAttribute)
                        .append(isCalendarDate(originAttribute) ? " (" + origin + ")" : "").append(", ")
                        .append(Math.abs(originHeight)).append(" years ")
                        .append(originHeight >= 0.0 ? "before" : "after")
                        .append(" the most recent tip (").append(mostRecent.getTimeValue()).append(")\n");
            }
        }

        return transitionTimes;
    }

    /**
     * The transition times from date elements inside the element. The origin is the date the epochs are
     * counted back from: the origin attribute, with every date a transition, or without it the most recent
     * of the dates, with the others the transitions. Each transition is at the height of its date, the time
     * before the most recent tip; if no taxon has a date, the origin is taken as the time of the most recent
     * tip.
     */
    private static double[] getTransitionTimesFromDates(XMLObject xo, List<Date> dates, StringBuilder notes)
            throws XMLParseException {

        for (String attribute : new String[]{EPOCH_COUNT, DURATION, TRANSITION_TIMES, UNITS}) {
            if (xo.hasAttribute(attribute)) {
                throw new XMLParseException("The " + EPOCHS + " element" + idText(xo) + " gives its epochs as " +
                        "dates, so it should not have the " + attribute + " attribute");
            }
        }

        final Date first = dates.get(0);
        for (Date date : dates) {
            if (!Arrays.asList(TIME_UNIT_TYPES).contains(date.getUnits())) {
                throw new XMLParseException("In the " + EPOCHS + " element" + idText(xo) + ", the dates must be " +
                        "calendar dates or have units of time (days, weeks, months or years)");
            }
            if (date.getUnits() != first.getUnits() || date.isBackwards() != first.isBackwards()) {
                throw new XMLParseException("In the " + EPOCHS + " element" + idText(xo) + ", the dates are not " +
                        "all in the same units and direction");
            }
        }

        // the origin, on the scale of the dates
        final double origin;
        final List<Date> transitions = new ArrayList<Date>(dates);
        if (xo.hasAttribute(ORIGIN)) {
            final String originAttribute = xo.getStringAttribute(ORIGIN).trim();
            if (isCalendarDate(originAttribute)) {
                if (first.getUnits() != Units.Type.YEARS || first.isBackwards()) {
                    throw new XMLParseException("The " + EPOCHS + " element" + idText(xo) + " has a calendar " +
                            "date for the " + ORIGIN + ", which needs its dates to be calendar dates too, or in " +
                            "years running forwards (such as 2020.5)");
                }
                origin = parseCalendarDate(originAttribute, xo);
            } else {
                try {
                    origin = Double.parseDouble(originAttribute);
                } catch (NumberFormatException nfe) {
                    throw new XMLParseException("In the " + EPOCHS + " element" + idText(xo) + ", " + ORIGIN +
                            "=\"" + originAttribute + "\" is neither a number nor a calendar date (yyyy-mm-dd)");
                }
            }
            notes.append("  origin ").append(originAttribute);
        } else {
            // the most recent of the dates: going backwards it has the smallest value, going forwards the largest
            Date mostRecent = first;
            for (Date date : dates) {
                if (date.isBackwards() ? date.getTimeValue() < mostRecent.getTimeValue() :
                        date.getTimeValue() > mostRecent.getTimeValue()) {
                    mostRecent = date;
                }
            }
            transitions.remove(mostRecent);
            origin = mostRecent.getTimeValue();
            notes.append("  origin ").append(origin).append(", the most recent of the dates");
        }

        // the time from the most recent tip back to the origin
        double originHeight = 0.0;
        final boolean hasDatedTaxa = Taxon.getMostRecentDate() != null;
        if (hasDatedTaxa) {
            final Date mostRecentTip = getMostRecentDate(xo);
            if (mostRecentTip.getUnits() != first.getUnits() || mostRecentTip.isBackwards() != first.isBackwards()) {
                throw new XMLParseException("In the " + EPOCHS + " element" + idText(xo) + ", the dates are not " +
                        "in the same units and direction as the dates of the taxa");
            }
            originHeight = TimeScale.convertTimeUnits(timeBefore(mostRecentTip, origin), first.getUnits(), TREE_UNITS);
            notes.append(", ").append(Math.abs(originHeight)).append(" years ")
                    .append(originHeight >= 0.0 ? "before" : "after")
                    .append(" the most recent tip (").append(mostRecentTip.getTimeValue()).append(")\n");
        } else {
            notes.append(", taken as the time of the most recent tip\n");
        }

        double[] transitionTimes = new double[transitions.size()];
        for (int i = 0; i < transitionTimes.length; i++) {
            // how long before the origin the date is
            final double timeBeforeOrigin = -timeBefore(transitions.get(i), origin);
            transitionTimes[i] = originHeight +
                    TimeScale.convertTimeUnits(timeBeforeOrigin, first.getUnits(), TREE_UNITS);
        }
        Arrays.sort(transitionTimes);

        if (transitionTimes.length > 0 && !(transitionTimes[0] > 0.0)) {
            throw new XMLParseException("In the " + EPOCHS + " element" + idText(xo) + ", the dates of the " +
                    "transitions between epochs must be before " +
                    (hasDatedTaxa ? "the most recent tip" : "the " + ORIGIN) +
                    (transitionTimes[0] == 0.0 ? ": the latest is at the same time as it" :
                            ": the latest is " + (-transitionTimes[0]) + " years after it"));
        }
        return transitionTimes;
    }

    /**
     * @param date  a date
     * @param value a value on the scale of the date
     * @return how long before the date the value is (negative if it is after), in the units of the date.
     * Dates run forwards in time (2020.5) or backwards (an age of 3.5).
     */
    private static double timeBefore(Date date, double value) {
        return date.isBackwards() ? value - date.getTimeValue() : date.getTimeValue() - value;
    }

    /**
     * @return the times converted from one unit of time to another
     */
    public static double[] convert(double[] times, Units.Type fromUnits, Units.Type toUnits) {
        double[] converted = new double[times.length];
        for (int i = 0; i < times.length; i++) {
            converted[i] = TimeScale.convertTimeUnits(times[i], fromUnits, toUnits);
        }
        return converted;
    }

    /**
     * @return the tokens of the transitionTimes attribute if they are calendar dates, or null if they are
     * numbers
     */
    private static String[] getCalendarDates(XMLObject xo) throws XMLParseException {
        final String value = xo.getStringAttribute(TRANSITION_TIMES).trim();
        if (value.isEmpty()) {
            return null;
        }
        final String[] tokens = value.split("[\\s,]+");
        int dateCount = 0;
        for (String token : tokens) {
            if (isCalendarDate(token)) {
                dateCount++;
            }
        }
        if (dateCount == 0) {
            return null;
        }
        if (dateCount < tokens.length) {
            throw new XMLParseException("In the " + EPOCHS + " element" + idText(xo) + ", the " +
                    TRANSITION_TIMES + " must be either all numbers or all calendar dates (yyyy-mm-dd)");
        }
        return tokens;
    }

    private static boolean isCalendarDate(String value) {
        return value.matches("\\d\\d\\d\\d-\\d\\d-\\d\\d");
    }

    /**
     * @return the calendar date as a decimal year, by the conversion BEAUti uses for the dates of tips
     */
    private static double parseCalendarDate(String value, XMLObject xo) throws XMLParseException {
        final SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd");
        format.setLenient(false);
        format.setTimeZone(TimeZone.getTimeZone("GMT"));
        try {
            return new Date(format.parse(value)).getTimeValue();
        } catch (ParseException pe) {
            throw new XMLParseException("In the " + EPOCHS + " element" + idText(xo) + ", \"" + value +
                    "\" is not a calendar date");
        }
    }

    /**
     * @return the date of the most recent tip, which for calendar dates to be placed relative to it must be
     * in years running forwards
     */
    private static Date getMostRecentCalendarDate(XMLObject xo, String purpose) throws XMLParseException {
        final Date mostRecent = getMostRecentDate(xo);
        if (mostRecent.getUnits() != Units.Type.YEARS || mostRecent.isBackwards()) {
            throw new XMLParseException("The " + EPOCHS + " element" + idText(xo) + " has " + purpose +
                    ", which needs the dates of the taxa to be in years running forwards (such as 2020.5)");
        }
        return mostRecent;
    }

    /**
     * @return the date of the most recent tip, which is height zero: the most recent of the dates of all the
     * taxa that have been read (they are kept by Taxon as they are read, not for one set of taxa or tree)
     */
    private static Date getMostRecentDate(XMLObject xo) throws XMLParseException {
        final Date mostRecent = Taxon.getMostRecentDate();
        if (mostRecent == null) {
            throw new XMLParseException("The " + EPOCHS + " element" + idText(xo) + " has dates to place " +
                    "before the most recent tip, but no taxon read before it has a date");
        }
        if (!Arrays.asList(TIME_UNIT_TYPES).contains(mostRecent.getUnits())) {
            throw new XMLParseException("In the " + EPOCHS + " element" + idText(xo) + ", the dates of the " +
                    "taxa are not in units of time, so dates cannot be placed before the most recent tip");
        }
        return mostRecent;
    }

    private static Units.Type parseTimeUnits(String name, XMLObject xo) throws XMLParseException {
        for (int i = 0; i < TIME_UNIT_NAMES.length; i++) {
            if (TIME_UNIT_NAMES[i].equalsIgnoreCase(name.trim())) {
                return TIME_UNIT_TYPES[i];
            }
        }
        throw new XMLParseException("In the " + EPOCHS + " element" + idText(xo) + ", " + UNITS + "=\"" +
                name + "\" is not a unit of time (days, weeks, months or years)");
    }

    private static String idText(XMLObject xo) {
        try {
            return xo.hasId() ? " '" + xo.getId() + "'" : "";
        } catch (XMLParseException e) {
            return "";
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Epochs of time with fixed transition times, given as a number of epochs of the same duration " +
                "or as the transition times.";
    }

    public Class getReturnType() {
        return DefaultEpochProvider.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(EPOCH_COUNT, true,
                    "The number of epochs, counting the last, which has no end"),
            AttributeRule.newDoubleRule(DURATION, true,
                    "The duration of each epoch"),
            new StringAttributeRule(TRANSITION_TIMES,
                    "The times at which one epoch gives way to the next, in increasing order, or their " +
                            "calendar dates (yyyy-mm-dd) in any order", true),
            new StringAttributeRule(UNITS,
                    "The units of the duration or transition times (years if not given)",
                    TIME_UNIT_NAMES, true),
            new StringAttributeRule(ORIGIN,
                    "A date, on the scale of the dates of the taxa or as a calendar date (yyyy-mm-dd), that the " +
                            "transition times are counted back from (by default, the date of the most recent tip)",
                    true),
            new ElementRule(TRANSITION_TIMES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, "A parameter for the transition times. If it has values they are the transition times (as heights, " +
                    "in years) and nothing else is needed; if not, it is given the transition times, so that " +
                    "operators can be applied to them", true),
            new ElementRule(Date.class,
                    "Dates that give the epochs: the transitions between them and, without an origin attribute, " +
                            "the origin (the most recent of the dates)", 0, Integer.MAX_VALUE)
    };
}
