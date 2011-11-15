@echo off
for /f "delims= " %%i in ('dir /b .\*~') do copy /y "%%i" ".\backup"
for /f "delims= " %%i in ('dir /b .\*~') do del "%%i"
pause