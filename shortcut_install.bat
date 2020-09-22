@echo off
@Set "SCRIPT=%TEMP%\LinkMaker-%RANDOM%-%RANDOM%.vbs"
@(  echo Set oWS = WScript.CreateObject("WScript.Shell"^)
    echo sLinkFile = "%USERPROFILE%\Desktop\tReasure.lnk"
    echo Set oLink = oWS.CreateShortcut(sLinkFile^)
    echo oLink.TargetPath = "C:\Program Files\R\R-4.0.2\Rscript.exe"
    echo oLink.Arguments = """%USERPROFILE%\Documents\tReasure_test.R"""
    echo oLink.IconLocation = "%USERPROFILE%\Downloads\treasure.ico"
    echo oLink.WorkingDirectory = "%USERPROFILE%\Documents"
    echo oLink.Save
)>"%SCRIPT%"
@"%__AppDir__%cscript.exe" //NoLogo "%SCRIPT%"
@Del "%SCRIPT%"