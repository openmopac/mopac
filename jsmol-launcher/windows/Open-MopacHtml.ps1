<#
Open-MopacHtml.ps1  -  open a MOPAC-generated HTML file with JSmol, no setup needed.

How it works
  1. Starting from the folder that holds the HTML file, walk up the folder tree
     until a folder containing  jsmol\JSmol.min.js  is found.  That folder is the
     "workspace root".  One copy of JSmol anywhere above your working folders is
     enough; nothing is copied or linked.
  2. Start a small local web server (serve.ps1) rooted at that folder, unless one
     for the same root is already running.  The server only listens on
     localhost and needs no administrator rights.
  3. Open  http://localhost:<port>/<path to the file>  in the default browser.

This matches the path convention used by MOPAC 23.2.1 and later, which writes
  <script src="/jsmol/JSmol.min.js">  into its HTML output.  It also works with
older MOPAC output that uses  ../jsmol/  provided the jsmol folder is in the
expected place relative to the HTML file.

Usage
  Open-MopacHtml.ps1  "C:\path\to\molecule.html"
  Normally invoked through Open_MOPAC_HTML.cmd (drag-and-drop or "Open with").
#>
param(
    [Parameter(Mandatory = $true, Position = 0)]
    [string]$HtmlFile,
    [int]$Port = 8000
)

$ErrorActionPreference = 'Stop'

function Fail($msg) {
    Write-Host ""
    Write-Host "  $msg" -ForegroundColor Red
    Write-Host ""
    Read-Host "  Press Enter to close"
    exit 1
}

# ---- 1. resolve the file and find the workspace root -------------------------
try { $file = (Resolve-Path -LiteralPath $HtmlFile).Path }
catch { Fail "File not found: $HtmlFile" }

$dir = Split-Path -Parent $file
$root = $null
$probe = $dir
while ($probe) {
    if (Test-Path -LiteralPath (Join-Path $probe 'jsmol\JSmol.min.js')) { $root = $probe; break }
    $parent = Split-Path -Parent $probe
    if (-not $parent -or $parent -eq $probe) { break }
    $probe = $parent
}
if (-not $root -and $env:MOPAC_JSMOL_ROOT -and
    (Test-Path -LiteralPath (Join-Path $env:MOPAC_JSMOL_ROOT 'jsmol\JSmol.min.js'))) {
    $root = $env:MOPAC_JSMOL_ROOT
}
if (-not $root) {
    Fail ("No 'jsmol' folder was found in $dir or any folder above it.`n" +
          "  Put one copy of the JSmol 'jsmol' folder (the one containing JSmol.min.js)`n" +
          "  in a folder above your MOPAC working folders, for example C:\Users\<you>\MOPAC\jsmol,`n" +
          "  or set the environment variable MOPAC_JSMOL_ROOT to that parent folder.")
}

# ---- 2. make sure a server rooted at $root is running -------------------------
$stateDir  = Join-Path $env:LOCALAPPDATA 'mopac-jsmol'
$stateFile = Join-Path $stateDir 'server.json'
New-Item -ItemType Directory -Force -Path $stateDir | Out-Null

$serverOk = $false
if (Test-Path -LiteralPath $stateFile) {
    try {
        $st = Get-Content -LiteralPath $stateFile -Raw | ConvertFrom-Json
        $proc = Get-Process -Id $st.pid -ErrorAction SilentlyContinue
        if ($proc -and $st.root -eq $root) { $serverOk = $true; $Port = $st.port }
        elseif ($proc) { Stop-Process -Id $st.pid -Force -ErrorAction SilentlyContinue }
    } catch { }
}

if (-not $serverOk) {
    # pick a free port starting at $Port
    $p = $Port
    while ($p -lt $Port + 20) {
        $busy = Get-NetTCPConnection -LocalPort $p -State Listen -ErrorAction SilentlyContinue
        if (-not $busy) { break }
        $p++
    }
    $Port = $p
    $serve = Join-Path $PSScriptRoot 'serve.ps1'
    $psArgs = @('-NoProfile', '-ExecutionPolicy', 'Bypass', '-WindowStyle', 'Minimized',
              '-File', "`"$serve`"", '-Root', "`"$root`"", '-Port', $Port)
    $proc = Start-Process -FilePath 'powershell.exe' -ArgumentList $psArgs -PassThru -WindowStyle Minimized
    @{ pid = $proc.Id; root = $root; port = $Port } | ConvertTo-Json | Set-Content -LiteralPath $stateFile
    # wait until it answers
    $deadline = (Get-Date).AddSeconds(10)
    do {
        Start-Sleep -Milliseconds 200
        try { $null = Invoke-WebRequest -Uri "http://localhost:$Port/jsmol/JSmol.min.js" -Method Head -UseBasicParsing -TimeoutSec 2; $serverOk = $true } catch { }
    } while (-not $serverOk -and (Get-Date) -lt $deadline)
    if (-not $serverOk) { Fail "The local web server did not start (port $Port)." }
}

# ---- 3. open the page ----------------------------------------------------------
$rel = $file.Substring($root.Length).TrimStart('\').Replace('\', '/')
$relEscaped = ($rel -split '/' | ForEach-Object { [uri]::EscapeDataString($_) }) -join '/'
$url = "http://localhost:$Port/$relEscaped"
Write-Host "  JSmol root : $root"
Write-Host "  Opening    : $url"
Start-Process $url
