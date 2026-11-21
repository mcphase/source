@echo off
set "_MCTEST_HELPER=%~1"
if not defined _MCTEST_HELPER (
    echo mctest_helper: missing helper name 1>&2
    exit /b 2
)

if /i not "%_MCTEST_HELPER:~-3%"==".pl" set "_MCTEST_HELPER=%_MCTEST_HELPER%.pl"
if not defined MCTEST_PERL set "MCTEST_PERL=perl"

set "_MCTEST_OUT=%TEMP%\mctest_helper_%RANDOM%.tmp"

"%MCTEST_PERL%" "%MCPHASE_DIR%\bin\%_MCTEST_HELPER%" %2 %3 %4 %5 %6 %7 %8 %9 > "%_MCTEST_OUT%"
if errorlevel 1 (
    set "_MCTEST_RC=1"
    del "%_MCTEST_OUT%" 2>nul
    set "_MCTEST_HELPER="
    exit /b 1
)

for /f "usebackq tokens=1,* delims==" %%A in ("%_MCTEST_OUT%") do set "%%A=%%B"

del "%_MCTEST_OUT%" 2>nul
set "_MCTEST_HELPER="
exit /b 0
