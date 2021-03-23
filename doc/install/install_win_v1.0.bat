:: Require files
:: pkg_lib.R
:: loadR_v1.R
:: tReasure.R
:: tReasure_source.tar.gz


@Echo Off
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: move files on Documents
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

mkdir %USERPROFILE%\Documents\tReasure_v1
copy * %USERPROFILE%\Documents\tReasure_v1

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set the path Rscript.exe 
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
set "src=c:\Program Files"
set "search=Rscript.exe"

for /r "%src%" %%F in (*%search%*) do (
  set "full=%%~fF"
)
echo %full%


:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: install the tReasure packages
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
for /r "%USERPROFILE%\Documents\tReasure_v1" %%F in (*loadR_v1.R*) do (
 echo %%~fF
 set "load=%%~fF"
 "%full%" "%%~fF" %*
)

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set the .lnk on Desktop for tReasure 
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


set lnk=%userprofile%\desktop\tReasure.lnk

set command=^
$objshell = New-object -ComObject WScript.Shell;^
$olnk =$objshell.CreateShortcut('%lnk%');^
$olnk.TargetPath = '%full%';^
$olnk.Arguments = '%USERPROFILE%\Documents\tReasure_v1\tReasure_v1.R';^
$olnk.IconLocation = '%USERPROFILE%\Documents\tReasure_v1\tReasure.ico';^
$olnk.Save();

powershell -command "& {%command%}"




