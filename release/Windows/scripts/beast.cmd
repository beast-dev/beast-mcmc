@echo off
setlocal
if ""=="%BEAST%" set BEAST=%~dp0%..
set BEAST_LIB=%BEAST%\lib
set BEAST_PLUG=%BEAST%\plugins
rem The heap. A fixed 256MB was not enough for anything but the smallest analysis. Half the
rem machine's RAM keeps pace with whatever BEAST is run on, and is a ceiling rather than a
rem reservation. Set BEAST_XMX to fix a size instead, e.g. BEAST_XMX=8g.
set HEAP=-XX:MaxRAMPercentage=50
if not "%BEAST_XMX%"=="" set HEAP=-Xmx%BEAST_XMX%
java -Xms64m %HEAP% -Dbeast.plugins.dir="%BEAST_PLUG%" -Djava.library.path="%BEAST_LIB%" -jar "%BEAST_LIB%/beast.jar" %*