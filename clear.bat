@echo off
for /f "delims= " %%i in ('dir /b .\*~') do del "%%i"
del *.o

