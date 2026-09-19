param(
    [Parameter(Mandatory=$true)]
    [string]$TexPath,
    [string]$MatrixPath = (Join-Path $PSScriptRoot 'appendix_all_3vertex_cases.txt'),
    [string]$IssuePath = (Join-Path $PSScriptRoot 'appendix_case_extraction_issues.txt'),
    [string]$IdReferencePath = (Join-Path $PSScriptRoot 'appendix_all_3vertex_cases.txt')
)

$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'appendix_case_ids.ps1')
$caseIds = Read-AppendixCaseIds $IdReferencePath
$seenIds = @{}

function Convert-PhaseVector([string]$Expr) {
    # Convert a heading exponent such as
    #   2\mathrm{i} \delta, \frac{\mathrm{i} \eee}{2}+\frac{2\mathrm{i} \pi}{f}
    # into coefficients of (a,b,c,d,e,q), where angle=pi*variable and q=1/f.
    $s = $Expr.Replace('\mathrm{i}', '')
    $frac = '\\frac\{([^{}]*)\}\{([^{}]*)\}'
    while ($s -match $frac) {
        $s = [regex]::Replace($s, $frac, '($1)/($2)', 1)
    }
    $s = $s.Replace('\aaa','a').Replace('\bbb','b').Replace('\ccc','c')
    $s = $s.Replace('\ddd','d').Replace('\delta','d').Replace('\eee','e')
    $s = $s.Replace('\epsilon','e').Replace('\pi/f','q')
    $s = $s.Replace(' ','').Replace('{','').Replace('}','').Replace('*','')
    $s = $s.Replace('(','').Replace(')','')
    $s = $s.Replace('\pi/f','q').Replace('\pi/3f','q/3')
    $s = $s.Replace('\pi/5f','q/5').Replace('\pi/7f','q/7')
    $numCoef = @(0,0,0,0,0,0)
    $denCoef = @(1,1,1,1,1,1)
    foreach ($term0 in ($s -split '(?=[+-])' | Where-Object { $_ -ne '' })) {
        $term = $term0
        $sign = 1
        if ($term.StartsWith('-')) { $sign=-1; $term=$term.Substring(1) }
        elseif ($term.StartsWith('+')) { $term=$term.Substring(1) }
        $tm = [regex]::Match($term, '^(?<num>[0-9]+)?(?<var>[abcdeq])(?:/(?<den>[0-9]+))?$')
        if (-not $tm.Success) { throw "unparsed phase term '$term0' from '$Expr' -> '$s'" }
        $num = if ($tm.Groups['num'].Success) {[int]$tm.Groups['num'].Value} else {1}
        $den = if ($tm.Groups['den'].Success) {[int]$tm.Groups['den'].Value} else {1}
        $idx = switch ($tm.Groups['var'].Value) { 'a' {0} 'b' {1} 'c' {2} 'd' {3} 'e' {4} 'q' {5} }
        $newNum = $numCoef[$idx]*$den + $sign*$num*$denCoef[$idx]
        $newDen = $denCoef[$idx]*$den
        $g = [math]::Abs([int]$newNum)
        $h = [math]::Abs([int]$newDen)
        while ($h -ne 0) { $tmp=$g % $h; $g=$h; $h=$tmp }
        if ($g -eq 0) { $g=1 }
        $numCoef[$idx] = [int]($newNum/$g)
        $denCoef[$idx] = [int]($newDen/$g)
    }
    $formatted = for ($j=0; $j -lt 6; $j++) {
        if ($denCoef[$j] -eq 1) { [string]$numCoef[$j] }
        else { ('{0}/{1}' -f $numCoef[$j],$denCoef[$j]) }
    }
    return '[' + ($formatted -join ',') + ']'
}

function Remove-TeXWrappers([string]$Text) {
    $result = $Text
    foreach ($marker in @('\redchange{','\sout{','\textcolor{red}')) {
        while ($true) {
            $start = $result.IndexOf($marker)
            if ($start -lt 0) { break }
            $open = if ($marker.EndsWith('{')) { $start + $marker.Length - 1 } else { $start + $marker.Length }
            if ($open -ge $result.Length -or $result[$open] -ne '{') { break }
            $depth = 0
            $close = -1
            for ($k=$open; $k -lt $result.Length; $k++) {
                if ($result[$k] -eq '{') { $depth++ }
                elseif ($result[$k] -eq '}') {
                    $depth--
                    if ($depth -eq 0) { $close = $k; break }
                }
            }
            if ($close -lt 0) { break }
            $inner = $result.Substring($open + 1, $close - $open - 1)
            $result = $result.Substring(0,$start) + $inner + $result.Substring($close + 1)
        }
    }
    return $result
}

function Convert-Monomial([string]$Part) {
    $v = @(0,0,0,0,0)
    $matches = [regex]::Matches($Part, '\\(aaa|bbb|ccc|ddd|eee)(?:\^\{?([0-9]+)\}?)?')
    foreach ($hit in $matches) {
        $idx = switch ($hit.Groups[1].Value) { 'aaa' {0} 'bbb' {1} 'ccc' {2} 'ddd' {3} 'eee' {4} }
        $power = if ($hit.Groups[2].Success) {[int]$hit.Groups[2].Value} else {1}
        $v[$idx] += $power
    }
    if ($matches.Count -eq 0) { throw "unparsed monomial '$Part'" }
    return '[' + ($v -join ',') + ']'
}

$rows = New-Object System.Collections.Generic.List[string]
$issues = New-Object System.Collections.Generic.List[string]
$headingNo = 0
$recordNo = 0
$twoVertexCount = 0
$lines = Get-Content -LiteralPath $TexPath

foreach ($line in $lines) {
    # Lines beginning with %% are intentionally disabled cases in the TeX source.
    if ($line -match '^\s*%') { continue }
    if ($line -notmatch '^\s*\\subsubsection\*\{Case\.?\s') { continue }
    $line = Remove-TeXWrappers $line
    $headingNo++
    try {
        $heading = Get-AppendixHeading $line $caseIds
        if ($heading.VertexCount -eq 2) { $twoVertexCount++; continue }
        $body = $heading.Body

        $phaseStart = $line.IndexOf(', $(')
        if ($phaseStart -lt 0) { throw 'no phase basis' }
        $phase = $line.Substring($phaseStart + 4)
        $pm = [regex]::Match($phase, 'x=\{\\mathrm e\}\^\{(?<x>.*?)\},y=\{\\mathrm e\}\^\{(?<y>.*?)\}\)')
        if (-not $pm.Success) { throw 'phase basis did not match' }
        $xvec = Convert-PhaseVector $pm.Groups['x'].Value
        $yvec = Convert-PhaseVector $pm.Groups['y'].Value

        $orMarker = '\,(\text{or}\,'
        $orPos = $body.IndexOf($orMarker)
        $primary = if ($orPos -ge 0) { $body.Substring(0,$orPos) } else { $body }
        $parts = @($primary -split ',') | Where-Object { $_.Trim().Length -gt 0 }
        if ($parts.Count -ne 3) { throw "primary monomial count=$($parts.Count)" }
        $baseVectors = @($parts | ForEach-Object { Convert-Monomial $_ })

        $variants = New-Object System.Collections.Generic.List[string]
        $variants.Add(('[' + ($baseVectors -join ',') + ']'))
        if ($orPos -ge 0) {
            $orText = $body.Substring($orPos + $orMarker.Length).TrimEnd(')')
            foreach ($alt in (@($orText -split ',') | Where-Object {$_.Trim().Length -gt 0})) {
                $av = Convert-Monomial $alt
                $variants.Add(('[' + $baseVectors[0] + ',' + $baseVectors[1] + ',' + $av + ']'))
            }
        }
        $branch = 0
        foreach ($vv in $variants) {
            $recordNo++
            $key = Get-AppendixVertexKey $vv
            if (-not $caseIds.ContainsKey($key)) { throw "Unknown variant for $($heading.Id); update the reference matrix first" }
            $id = $caseIds[$key]
            if ($seenIds.ContainsKey($id)) { throw "Duplicate heading or variant: $id" }
            $seenIds[$id] = $true
            $rows.Add(('  ["{0}", {1}, {2}, {3}]' -f $id,$vv,$xvec,$yvec))
            $branch++
        }
    }
    catch {
        $issues.Add(("HEADING-{0}`t{1}`t{2}" -f $headingNo,$_.Exception.Message,$line))
    }
}

foreach ($id in $caseIds.Values) {
    if (-not $seenIds.ContainsKey($id)) { $issues.Add("MISSING`t$id") }
}
Set-Content -LiteralPath $IssuePath -Value $issues -Encoding ASCII
if ($issues.Count -gt 0) { throw "Extraction failed with $($issues.Count) issues; matrix was not overwritten. See $IssuePath" }

$header = @(
    '# Auto-extracted from active Appendix case headings.',
    '# Each record is [case_id, three 5D vertex vectors, x phase vector, y phase vector].',
    '# Phase vectors use normalized variables (a,b,c,d,e,q), with q=1/f.',
    '# Vertex identities determine the IDs; APP-133 covers the repeated heading APP-158.',
    'appendix_all_3vertex_cases := ['
)
$footer = @(']:')
$content = $header + ((($rows | Sort-Object) -join ",`n") + "`n") + $footer
# Maple's command-line reader on Windows rejects the UTF-8 BOM emitted by
# Windows PowerShell's UTF8 mode, so these ASCII-only files are BOM-free.
Set-Content -LiteralPath $MatrixPath -Value $content -Encoding ASCII
Set-Content -LiteralPath $IssuePath -Value $issues -Encoding ASCII
Write-Output ("active_headings={0}; three_vertex_headings={1}; two_vertex_headings={2}; expanded_records={3}; issues={4}" -f $headingNo,($headingNo-$twoVertexCount),$twoVertexCount,$rows.Count,$issues.Count)
