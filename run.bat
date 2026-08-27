@echo off
setlocal
cd /d "%~dp0"

set "RSCRIPT="

REM First, look for Rscript on the system PATH
for /f "delims=" %%I in ('where Rscript.exe 2^>nul') do (
    set "RSCRIPT=%%I"
    goto :found
)

REM Otherwise, search the normal R installation folder
for /d %%D in ("%ProgramFiles%\R\R-*") do (
    if exist "%%D\bin\Rscript.exe" (
        set "RSCRIPT=%%D\bin\Rscript.exe"
    )
)

:found
if not defined RSCRIPT (
    echo ERROR: R could not be found.
    echo.
    echo Please install R from:
    echo https://cran.r-project.org/
    echo.
    pause
    exit /b 1
)

echo Using R:
echo %RSCRIPT%
echo.

"%RSCRIPT%" "%~dp0run_app.R"

if errorlevel 1 (
    echo.
    echo The app stopped because of an error.
)

pause