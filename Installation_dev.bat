:: Check for Python Installation
python --version 3>NUL
if errorlevel 1 goto errorNoPython
echo Python is installed! Continue with installation of packages.
python3 -m ensurepip --default-pip
pip install -r requirements.txt
echo "Installation successfully finished"
pause

:: skip executing the errorNoPython section
goto:eof

:errorNoPython
echo.
echo Error^: Python not installed. Please install python from the Windows store and press any key to continue afterwards.
python3
pause
python3 -m ensurepip --default-pip
pip install -r requirements.txt
echo "Installation successfully finished"
pause