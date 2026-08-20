Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$environmentBefore = @{}
Get-ChildItem Env: | ForEach-Object {
    $environmentBefore[$_.Name] = $_.Value
}

$vswhere = Join-Path ${env:ProgramFiles(x86)} `
    "Microsoft Visual Studio\Installer\vswhere.exe"
$installationPath = & $vswhere -latest -products * `
    -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
    -property installationPath
$installationPath = [string]($installationPath | Select-Object -First 1)
if ([string]::IsNullOrWhiteSpace($installationPath)) {
    throw "Visual Studio C++ tools not found"
}

$developerShell = Join-Path $installationPath `
    "Common7\Tools\Launch-VsDevShell.ps1"
& $developerShell -Arch amd64 -HostArch amd64 -SkipAutomaticLocation

Get-ChildItem Env: | Where-Object {
    -not $environmentBefore.ContainsKey($_.Name) -or
    $environmentBefore[$_.Name] -ne $_.Value
} | ForEach-Object {
    $delimiter = "DSEAMS_$([guid]::NewGuid().ToString('N'))"
    Add-Content -Path $env:GITHUB_ENV -Value "$($_.Name)<<$delimiter"
    Add-Content -Path $env:GITHUB_ENV -Value $_.Value
    Add-Content -Path $env:GITHUB_ENV -Value $delimiter
}
