@echo off
setlocal
if ""=="%BEAST%" set BEAST=%~dp0%..
set BEAST_LIB=%BEAST%\lib
java -Xms64m -Xmx1024m -Djava.library.path="%BEAST_LIB%" -cp "%BEAST_LIB%/beast.jar" dr.app.beast.BeastMC3 %*
