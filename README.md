# Thermal-Tolerance-of-Acer-Campestre
Thermal tolerance of  Acer campestre (field maple) under heat and drought stress derived from chlorophyll fluorescence

This file contains raw chlorophyll fluorescence measurements collected to assess the thermal tolerance of Acer campestre under controlled heat ramp conditions. The data were used to compute two physiological thermal thresholds:

Tcrit: The temperature at which PSII efficiency (Fv/Fm) begins to sharply decline.

T₅₀: The temperature at which 50% of PSII efficiency is lost.

These values were derived from fluorescence responses across temperature treatments, both pre- and post-hydration.

Data Columns
Column Name	Description
Sample #	Measurement ID
Date	Measurement date (MM/DD)
Time	Time of measurement (24-hr format)
Mod Int	Measuring light intensity (PAM setting)
Mod Wl	Measuring wavelength (e.g., "R" for red)
Det Gain	Detector gain setting
Sat Flash Int	Saturation flash intensity (%)
Sat PW	Saturation pulse width (ms)
Far Red Mode	Far-red light mode status (On/Off)
Far Red Int	Far-red light intensity
Far Red Dur	Duration of far-red light (s)
Fo	Minimum fluorescence (dark-adapted)
Fm	Maximum fluorescence (dark-adapted)
Fv	Variable fluorescence (Fm - Fo)
Fv/m	Maximum quantum yield of PSII (Fv/Fm)
Fv/o	Fluorescence ratio (Fv/Fo)
Temperature	Temperature during the measurement (°C)

Application
Used in generating thermal response curves of PSII.

Key thresholds (Tcrit and T₅₀) extracted using nonlinear regression and breakpoint analysis.

Values were further analyzed using PCA, clustering, and multivariate statistics to assess the effects of hydration and trial conditions.

Data collected using a standard PAM fluorometer under consistent protocols.

Each temperature point typically includes replicate measurements.

Outliers and anomalies (if any) were handled during preprocessing in R.

