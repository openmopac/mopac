<#
serve.ps1  -  minimal static file server for MOPAC HTML output and JSmol.

Serves the files under -Root at http://localhost:<Port>/ .  Listens only on
localhost, so nothing outside this computer can reach it.  Needs no
administrator rights and no Python.  Started automatically by
Open-MopacHtml.ps1; stop it with Stop_MOPAC_server.cmd or just close its window.
#>
param(
    [Parameter(Mandatory = $true)] [string]$Root,
    [int]$Port = 8000
)

$Root = (Resolve-Path -LiteralPath $Root).Path.TrimEnd('\')
$host.UI.RawUI.WindowTitle = "MOPAC JSmol server  -  $Root  -  port $Port"

$mime = @{
    '.html' = 'text/html; charset=utf-8'; '.htm' = 'text/html; charset=utf-8'
    '.js'   = 'application/javascript';   '.css' = 'text/css'
    '.json' = 'application/json';         '.png' = 'image/png'
    '.jpg'  = 'image/jpeg';               '.gif' = 'image/gif'
    '.svg'  = 'image/svg+xml';            '.ico' = 'image/x-icon'
    '.txt'  = 'text/plain; charset=utf-8'; '.pdb' = 'text/plain'
    '.xyz'  = 'text/plain';               '.mgf' = 'text/plain'
    '.mop'  = 'text/plain';               '.arc' = 'text/plain'
    '.out'  = 'text/plain';               '.spt' = 'text/plain'
    '.zip'  = 'application/zip';          '.gz'  = 'application/gzip'
    '.wasm' = 'application/wasm';         '.woff' = 'font/woff'; '.woff2' = 'font/woff2'
}

$listener = New-Object System.Net.HttpListener
$listener.Prefixes.Add("http://localhost:$Port/")
$listener.Start()
Write-Host "Serving $Root at http://localhost:$Port/   (close this window to stop)"

while ($listener.IsListening) {
    try { $ctx = $listener.GetContext() } catch { break }
    $req = $ctx.Request; $res = $ctx.Response
    try {
        $relPath = [uri]::UnescapeDataString($req.Url.AbsolutePath).TrimStart('/').Replace('/', '\')
        $full = [System.IO.Path]::GetFullPath((Join-Path $Root $relPath))
        $inside = $full.Equals($Root, [System.StringComparison]::OrdinalIgnoreCase) -or
                  $full.StartsWith($Root + '\', [System.StringComparison]::OrdinalIgnoreCase)
        if (-not $inside) {
            $res.StatusCode = 403
        } elseif (Test-Path -LiteralPath $full -PathType Container) {
            $index = Join-Path $full 'index.html'
            if (Test-Path -LiteralPath $index) { $full = $index } else { $res.StatusCode = 404 }
        } elseif (-not (Test-Path -LiteralPath $full -PathType Leaf)) {
            $res.StatusCode = 404
        }
        if ($res.StatusCode -eq 200) {
            $ext = [System.IO.Path]::GetExtension($full).ToLower()
            $res.ContentType = if ($mime.ContainsKey($ext)) { $mime[$ext] } else { 'application/octet-stream' }
            $bytes = [System.IO.File]::ReadAllBytes($full)
            $res.ContentLength64 = $bytes.Length
            if ($req.HttpMethod -ne 'HEAD') { $res.OutputStream.Write($bytes, 0, $bytes.Length) }
        }
    } catch {
        try { $res.StatusCode = 500 } catch { }
    } finally {
        try { $res.OutputStream.Close() } catch { }
    }
}
