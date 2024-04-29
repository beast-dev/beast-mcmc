@echo off
setlocal
if ""=="%BEAST%" set BEAST=%~dp0%..
set BEAST_LIB=%BEAST%\lib
set BEAST_PLUG=%BEAST%\plugins
java -Xms64m -Xmx256m -Dbeast.plugins.dir="%BEAST_PLUG%" -Djava.library.path="%BEAST_LIB%" -jar "%BEAST_LIB%/beast.jar" %*
