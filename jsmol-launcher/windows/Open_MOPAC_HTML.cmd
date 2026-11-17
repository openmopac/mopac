@echo off
REM Open_MOPAC_HTML.cmd - open a MOPAC-generated HTML file with JSmol.
REM
REM Use it in any of these ways:
REM   * drag an HTML file onto this .cmd file
REM   * right-click an HTML file > Open with > choose this .cmd file
REM   * from a command prompt:  Open_MOPAC_HTML.cmd  "C:\work\molecule.html"
REM
REM All the work is done by Open-MopacHtml.ps1 in the same folder.  Nothing is
REM installed, no administrator rights are needed, and no files are written
REM into your working folders.
if "%~1"=="" (
    echo.
    echo   Drag a MOPAC HTML file onto this file, or use "Open with".
    echo.
    pause
    exit /b 1
)
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0Open-MopacHtml.ps1" "%~1"
