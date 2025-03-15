@echo off
:: Ensure the script runs as administrator
net session >nul 2>&1
if %errorLevel% neq 0 (
    echo Requesting administrative privileges...
    powershell -Command "Start-Process '%0' -Verb RunAs"
    exit /b
)

:: Get the script's directory
set "CURRENT_DIR=%~dp0"
set "CURRENT_DIR=%CURRENT_DIR:~0,-1%"
cd /d "%CURRENT_DIR%"

:: Define the _internal directory
set "INTERNAL_DIR=%CURRENT_DIR%\_internal"

:: Unblock all DLLs inside _internal and output the list of affected files
echo Unblocking the following DLLs:
powershell -Command "Get-ChildItem -Path '%INTERNAL_DIR%' -Recurse -Filter *.dll | ForEach-Object { Unblock-File -Path $_.FullName; echo $_.FullName }"

:: Define the path to regasm.exe
set "REGASMPROG=%SystemRoot%\Microsoft.NET\Framework64\v4.0.30319\regasm.exe"

:: Unregister DLLs
%REGASMPROG% "%INTERNAL_DIR%\BaseCommon.dll" /tlb /u
%REGASMPROG% "%INTERNAL_DIR%\BaseDataAccess.dll" /tlb /u
%REGASMPROG% "%INTERNAL_DIR%\MassSpecDataReader.dll" /tlb /u

:: Register DLLs
%REGASMPROG% "%INTERNAL_DIR%\BaseCommon.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\BaseDataAccess.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\MassSpecDataReader.dll" /tlb


%REGASMPROG% "%INTERNAL_DIR%\Clearcore2.Compression.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\Clearcore2.Data.AnalystDataProvider.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\Clearcore2.Data.CommonInterfaces.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\Clearcore2.Data.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\Clearcore2.InternalRawXYProcessing.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\Clearcore2.Muni.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\Clearcore2.ProjectUtilities.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\Clearcore2.RawXYProcessing.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\Clearcore2.StructuredStorage.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\Clearcore2.Utility.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\Clearcore2.Data.WiffReader.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\Sciex.TofTof.T2DFMan.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\BlaisWiff.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\WiffReaderCOM.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\Sciex.ClearCore.FMAN.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\Sciex.Data.Processing.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\SCiex.Data.SimpleTypes.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\Data.XYData.dll" /tlb
%REGASMPROG% "%INTERNAL_DIR%\Sciex.FMan.dll" /tlb


:: Display the directory for user feedback
echo.
echo The batch file was called from: %CURRENT_DIR%
echo.


:: Delete the "Run" file inside _internal (if it exists)
del "%INTERNAL_DIR%\Run" >nul 2>&1


:: Launch GUI_UniDec.exe after the script completes
echo Launching GUI_UniDec...
start "" "%CURRENT_DIR%\GUI_UniDec.exe"

exit

