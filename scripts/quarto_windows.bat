@echo off
REM Some conda-forge Quarto builds retain a build-machine path for this DLL.
REM Resolve the installed runtime from the active environment instead.
if exist "%CONDA_PREFIX%\Library\lib\deno_dom.dll" set "QUARTO_DENO_DOM=%CONDA_PREFIX%\Library\lib\deno_dom.dll"
