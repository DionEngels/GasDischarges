rem Change the next line so the path matches your installation.
rem It is assumed that Dev-Cpp and gnuplot have been installed
rem in sub-directories of this directory.
rem Then call this file before you do any work.

set NPS_DIR=C:\CPPP
set DEVCPP_DIR=%NPS_DIR%\Dev-CppPortable\App\devcpp\
set GNUPLOT_DIR=%NPS_DIR%\gnuplot
set NPS_PATHS=%NPS_DIR%;%DEVCPP_DIR%;%DEVCPP_DIR%\bin;%GNUPLOT_DIR%\bin
set PATH=%NPS_PATHS%;%PATH%

