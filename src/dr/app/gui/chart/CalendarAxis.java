/*
 * LinearAxis.java
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

package dr.app.gui.chart;

import java.time.LocalDateTime;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.util.Calendar;

import static dr.evolution.util.TimeScale.DAYS_PER_YEAR;

public class CalendarAxis extends Axis.AbstractAxis {

	private static final DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yy-MM-dd");

	/**
	 * Empty constructor
	 */
	public CalendarAxis() { }

	/**
	 * Axis flag constructor
	 */
	public CalendarAxis(int minAxisFlag, int maxAxisFlag) {
		setAxisFlags(minAxisFlag, maxAxisFlag);
	}

	/**
	 *	Transform a value
	 */
	public double transform(double value) {
		return value;	// a linear transform !
	}

	/**
	*	Untransform a value
	*/
	public double untransform(double value) {
		return value;	// a linear transform !
	}

	@Override
	public String format(double value) {
		LocalDateTime dataTime = convertToDate(value);
		return dataTime.format(formatter);
	}

	public static LocalDateTime convertToDate(double decimalYear) {
		return LocalDateTime.ofEpochSecond(convertToEpochSeconds(decimalYear), 0, ZoneOffset.UTC);
	}

	private static long convertToEpochSeconds(double decimalYear) {
		return convertToEpochMilliseconds(decimalYear) / 1000;
	}

	private static long convertToEpochMilliseconds(double decimalYear) {
		int year = (int)Math.floor(decimalYear);
		long ms = (long)((decimalYear - Math.floor(decimalYear)) * DAYS_PER_YEAR * 24 * 3600 * 1000);
		Calendar calendar = Calendar.getInstance();
		calendar.set(year, 0, 0, 0, 0, 0);
		return (calendar.getTimeInMillis()) + ms;
	}
}
