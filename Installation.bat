@echo off
echo -----------------------------------------
echo.
echo "NOTE!: If the Microsoft store opens: Please search for "Python 3.9" and click "install" 
echo.
echo "This specific version is needed for PFAScreen (not the most recent)!"
echo.
pause
python3
pause
python3 -m ensurepip --default-pip
pip install -r requirements.txt
echo "Installation successfully finished!"
pause