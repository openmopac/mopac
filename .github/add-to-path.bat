@echo off

IF "%~1" == "" GOTO do_nothing

set NEW_PATH=%PATH%;%1
set PATH=%NEW_PATH%
if "%NEW_PATH:~1024%" == "" setx PATH "%NEW_PATH%"
set NEW_PATH=

:do nothing
