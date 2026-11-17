@echo off
REM Stops the local web server started by Open_MOPAC_HTML.cmd.
REM (Closing the minimized "MOPAC JSmol server" window does the same thing.)
powershell.exe -NoProfile -ExecutionPolicy Bypass -Command ^
  "$f = Join-Path $env:LOCALAPPDATA 'mopac-jsmol\server.json';" ^
  "if (Test-Path $f) { $s = Get-Content $f -Raw | ConvertFrom-Json;" ^
  "  Stop-Process -Id $s.pid -Force -ErrorAction SilentlyContinue; Remove-Item $f;" ^
  "  Write-Host ('Stopped server for ' + $s.root) } else { Write-Host 'No server is running.' }"
pause
