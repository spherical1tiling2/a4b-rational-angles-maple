param([Parameter(Mandatory=$true)][string]$TexPath)
$ErrorActionPreference = 'Stop'
$testRoot = Join-Path ([IO.Path]::GetTempPath()) ('a4b-extraction-' + [guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $testRoot | Out-Null
function Assert-Test([bool]$Condition, [string]$Message) {
    if (-not $Condition) { throw $Message }
}
function Read-Records([string]$Path) {
    return @(Get-Content -LiteralPath $Path | Where-Object { $_ -match '^\s*\["APP-' } | ForEach-Object { $_.Trim().TrimEnd(',') })
}
$matrix = Join-Path $testRoot 'matrix.txt'
$polys = Join-Path $testRoot 'polys.mpl'
& (Join-Path $PSScriptRoot 'extract_appendix_case_matrix.ps1') -TexPath $TexPath -MatrixPath $matrix -IssuePath (Join-Path $testRoot 'matrix-issues.txt')
& (Join-Path $PSScriptRoot 'extract_appendix_reference_polynomials.ps1') -TexPath $TexPath -OutputPath $polys -IssuePath (Join-Path $testRoot 'poly-issues.txt')
$records = Read-Records $matrix
$references = Read-Records $polys
Assert-Test ($records.Count -eq 257) 'Expected 257 extracted vertex combinations'
Assert-Test ($references.Count -eq 230) 'Expected 230 reference polynomials'
Assert-Test (@($records | Where-Object {$_ -match '"APP-158"'}).Count -eq 0) 'Retired ID returned'
Assert-Test (@($records | Where-Object {$_ -match '"APP-177-OR1"'}).Count -eq 1) 'Missing new OR branch'

# Both accepted punctuation forms must yield identical identities and data.
$plainTex = Join-Path $testRoot 'case-without-period.tex'
[IO.File]::WriteAllText($plainTex, (Get-Content -Raw -LiteralPath $TexPath).Replace('{Case. ','{Case '))
$plainMatrix = Join-Path $testRoot 'plain-matrix.txt'
& (Join-Path $PSScriptRoot 'extract_appendix_case_matrix.ps1') -TexPath $plainTex -MatrixPath $plainMatrix -IssuePath (Join-Path $testRoot 'plain-issues.txt')
Assert-Test (@(Compare-Object $records (Read-Records $plainMatrix)).Count -eq 0) 'Case punctuation changed extracted records'

# Duplicate headings must fail before replacing an existing output.
$firstHeading = @(Get-Content -LiteralPath $TexPath | Where-Object { $_ -match '^\s*\\subsubsection\*\{Case\.?\s' })[0]
$badTex = Join-Path $testRoot 'duplicate.tex'
[IO.File]::WriteAllText($badTex, (Get-Content -Raw -LiteralPath $TexPath) + "`n" + $firstHeading)
$before = (Get-FileHash -LiteralPath $matrix -Algorithm SHA256).Hash
$failed = $false
try {
    & (Join-Path $PSScriptRoot 'extract_appendix_case_matrix.ps1') -TexPath $badTex -MatrixPath $matrix -IssuePath (Join-Path $testRoot 'duplicate-issues.txt') | Out-Null
} catch { $failed = $true }
Assert-Test $failed 'Duplicate heading was accepted'
Assert-Test ((Get-FileHash -LiteralPath $matrix -Algorithm SHA256).Hash -eq $before) 'Failed extraction replaced the output'

# An omitted heading must also fail closed, for both extractors.
$missingTex = Join-Path $testRoot 'missing.tex'
[IO.File]::WriteAllText($missingTex, (Get-Content -Raw -LiteralPath $TexPath).Replace($firstHeading, '% ' + $firstHeading))
foreach ($kind in @('matrix','polynomial')) {
    $target = if ($kind -eq 'matrix') { $matrix } else { $polys }
    $before = (Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash
    $failed = $false
    try {
        if ($kind -eq 'matrix') {
            & (Join-Path $PSScriptRoot 'extract_appendix_case_matrix.ps1') -TexPath $missingTex -MatrixPath $target -IssuePath (Join-Path $testRoot 'missing-matrix-issues.txt') | Out-Null
        } else {
            & (Join-Path $PSScriptRoot 'extract_appendix_reference_polynomials.ps1') -TexPath $missingTex -OutputPath $target -IssuePath (Join-Path $testRoot 'missing-poly-issues.txt') | Out-Null
        }
    } catch { $failed = $true }
    Assert-Test $failed "Missing heading was accepted by $kind extractor"
    Assert-Test ((Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash -eq $before) "Failed $kind extraction replaced output"
}
Write-Output 'PASS: counts, stable IDs, Case/Case. compatibility, duplicate/missing rejection, and output preservation.'
Write-Output "Test artifacts: $testRoot"
