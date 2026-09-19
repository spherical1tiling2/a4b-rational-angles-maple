param(
    [Parameter(Mandatory=$true)][string]$TexPath,
    [string]$OutputPath = (Join-Path $PSScriptRoot 'appendix_3vertex_solution_candidates.txt')
)
$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'appendix_case_ids.ps1')
$ids = Read-AppendixCaseIds (Join-Path $PSScriptRoot 'appendix_all_3vertex_cases.txt')
$rows = [System.Collections.Generic.List[string]]::new()
$activeIds = @()
$inAppendix = $false
$headings = 0
foreach ($line in (Get-Content -LiteralPath $TexPath)) {
    if ($line -match '^\s*\\section\*\{Appendix') { $inAppendix = $true }
    if (-not $inAppendix) { continue }
    if ($line -match '\\sub(?:sub)?section\*\{(?:Subcase|Case)\.?\s') {
        $activeIds = @()
        if ($line -match '^\s*%' -or $line -match '\{Subcase') { continue }
        $heading = Get-AppendixHeading $line $ids
        if ($null -ne $heading -and $heading.VertexCount -eq 3) {
            $headings++
            $activeIds = @($ids.Values | Where-Object { $_ -eq $heading.Id -or $_.StartsWith($heading.Id+'-OR') } | Sort-Object)
        }
        continue
    }
    if ($activeIds.Count -eq 0 -or $line -match '^\s*%' -or $line -match '\\sout\{') { continue }
    $match = [regex]::Match($line, '^\s*\$f\s*=\s*(\d+)\s*,.*?=\s*\(([0-9,\s]+)\)\s*/\s*(\d+)\s*\$')
    if (-not $match.Success) { continue }
    $nums = $match.Groups[2].Value -replace '\s',''
    if (($nums -split ',').Count -ne 5) { throw 'Expected five angle numerators' }
    foreach ($id in $activeIds) { $rows.Add(($id+'|'+$match.Groups[1].Value+'|'+$nums+'|'+$match.Groups[3].Value)) }
}
if ($headings -ne 230 -or $rows.Count -eq 0) { throw 'Incomplete solution extraction' }
$text = @('# Finite angle tuples printed in the appendix; OR branches are expanded.') + @($rows | Sort-Object -Unique)
[IO.File]::WriteAllText([IO.Path]::GetFullPath($OutputPath), ($text -join "`n")+"`n", [Text.UTF8Encoding]::new($false))
Write-Output "primary_headings=$headings; finite_records=$($text.Count-1)"
