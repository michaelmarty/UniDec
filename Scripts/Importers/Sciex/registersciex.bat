@echo off
echo Registering WiffReaderCOM and related DLLs...

:: Ensure we are running as Administrator
net session >nul 2>&1
if %errorLevel% neq 0 (
    echo ERROR: Please run this script as Administrator!
    pause
    exit /b
)

:: Set the path to RegAsm (64-bit version)
set REGASM="C:\Windows\Microsoft.NET\Framework64\v4.0.30319\RegAsm"

:: Register .NET assemblies
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\Clearcore2.Compression.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\Clearcore2.Data.AnalystDataProvider.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\Clearcore2.Data.CommonInterfaces.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\Clearcore2.Data.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\Clearcore2.Data.WiffReader.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\Clearcore2.InternalRawXYProcessing.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\Clearcore2.Muni.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\Clearcore2.ProjectUtilities.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\Clearcore2.RawXYProcessing.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\Clearcore2.StructuredStorage.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\Clearcore2.Utility.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\Sciex.TofTof.T2DFMan.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\BlaisWiff.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\WiffReaderCOM.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\Sciex.ClearCore.FMAN.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\Sciex.Data.Processing.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\SCiex.Data.SimpleTypes.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\Data.XYData.dll" /tlb /codebase
%REGASM% "C:\Python\UniDec3\unidec\UniDecImporter\Sciex\Sciex.FMan.dll" /tlb /codebase


echo Registration complete!
pause
