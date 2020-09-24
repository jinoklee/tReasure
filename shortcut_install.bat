@echo off

set src=c:\Program Files
set search=Rscript.exe

for /r "%src%" %%F in (*%search%*) do (
  set full=%%~fF
  set name=%%~nxF
)

@Set "SCRIPT=%TEMP%\LinkMaker-%RANDOM%-%RANDOM%.vbs"
@(  echo Set oWS = WScript.CreateObject("WScript.Shell"^)
    echo sLinkFile = "%USERPROFILE%\Desktop\tReasure.lnk"
    echo Set oLink = oWS.CreateShortcut(sLinkFile^)
    echo oLink.TargetPath = "%full%"
    echo oLink.Arguments = """%USERPROFILE%\Documents\tReasure.R"""
    echo oLink.IconLocation = "%USERPROFILE%\Downloads\treasure.ico"
    echo oLink.WorkingDirectory = "%USERPROFILE%\Documents"
    echo oLink.Save
)>"%SCRIPT%"
@"%__AppDir__%cscript.exe" //NoLogo "%SCRIPT%"
@Del "%SCRIPT%"
