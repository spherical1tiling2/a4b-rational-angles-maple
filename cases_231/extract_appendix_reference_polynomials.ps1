param(
    [Parameter(Mandatory=$true)]
    [string]$TexPath,
    [string]$OutputPath = (Join-Path $PSScriptRoot 'appendix_reference_polys_by_heading.mpl'),
    [string]$IssuePath = (Join-Path $PSScriptRoot 'appendix_reference_polynomial_issues.txt')
)

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

$lines = Get-Content -LiteralPath $TexPath
$rows = New-Object System.Collections.Generic.List[string]
$issues = New-Object System.Collections.Generic.List[string]
$headingNo = 0

for ($i=0; $i -lt $lines.Count; $i++) {
    $line = $lines[$i]
    if ($line -match '^\s*%%') { continue }
    if ($line -notmatch '^\\subsubsection\*\{Case ') { continue }
    $headingNo++
    $eqLine = $null
    for ($j=$i+1; $j -lt [math]::Min($i+12,$lines.Count); $j++) {
        $normalized = Remove-TeXWrappers $lines[$j]
        if ($normalized -match '^\s*\$.*=0(\$|\s*$)') { $eqLine=$normalized; break }
        if ($normalized -match '\\begin\{aligned\}') {
            $block = New-Object System.Collections.Generic.List[string]
            $block.Add($lines[$j])
            for ($k=$j+1; $k -lt [math]::Min($i+40,$lines.Count); $k++) {
                $part = $lines[$k]
                $block.Add($part)
                if ($part -match '\\end\{aligned\}') { break }
            }
            $eqLine = Remove-TeXWrappers ($block -join ' ')
            break
        }
    }
    if ($null -eq $eqLine) {
        $issues.Add("HEADING-$headingNo`tno equation line")
        continue
    }
    $expr = $eqLine.Trim()
    $expr = $expr.Replace('\begin{aligned}','').Replace('\end{aligned}','').Replace('&','').Replace('\\','')
    $expr = $expr.Replace('$','')
    $eqPos = $expr.LastIndexOf('=0')
    $expr = $expr.Substring(0, $eqPos)
    $expr = $expr.Trim()
    $expr = $expr.Replace('\,','').Replace('\ ','').Replace('\!','').Replace('\quad','').Replace('\qquad','')
    $expr = $expr.Replace('\left','').Replace('\right','').Replace('\bigl','').Replace('\bigr','').Replace('\Bigl','').Replace('\Bigr','')
    $expr = [regex]::Replace($expr, '\\zeta_\{?([0-9]+)\}?', 'zeta$1')
    $expr = [regex]::Replace($expr, '\^\{([0-9]+)\}', '^$1')
    # TeX writes products by juxtaposition with spaces; Maple needs `*`.
    $expr = [regex]::Replace($expr, '\s*([+\-])\s*', '$1')
    $expr = [regex]::Replace($expr, '\s+', '*')
    $expr = [regex]::Replace($expr, '(x(?:\^[0-9]+)?)y', '$1*y')
    $expr = [regex]::Replace($expr, '(y(?:\^[0-9]+)?)x', '$1*x')
    $expr = [regex]::Replace($expr, '\)\(', ')*(')
    $expr = [regex]::Replace($expr, '\)\^([0-9]+)\(', ')^$1*(')
    $expr = [regex]::Replace($expr, '\)\s*([xy])', ')*$1')
    $expr = [regex]::Replace($expr, '\)(?=zeta)', ')*')
    $expr = [regex]::Replace($expr, '([0-9])\(', '$1*(')
    $expr = [regex]::Replace($expr, '([0-9])([xy])', '$1*$2')
    $expr = [regex]::Replace($expr, '([0-9])(?=zeta)', '$1*')
    # TeX uses juxtaposition for products, which Maple also accepts.
    $id = 'APP-{0:D3}' -f $headingNo
    $rows.Add(('  ["{0}", {1}]' -f $id,$expr))
}

$content = @(
    '# Reference polynomials copied from the active Appendix equations.',
    '# zetaN denotes exp(2*Pi*I/N).',
    'zeta3 := exp(2*Pi*I/3): zeta4 := I:',
    'zeta5 := exp(2*Pi*I/5): zeta6 := exp(2*Pi*I/6):',
    'zeta7 := exp(2*Pi*I/7): zeta8 := exp(2*Pi*I/8):',
    'zeta9 := exp(2*Pi*I/9): zeta12 := exp(2*Pi*I/12):',
    'zeta15 := exp(2*Pi*I/15): zeta20 := exp(2*Pi*I/20):',
    'appendix_reference_polys_by_heading := ['
)
$content += (($rows -join ",`n") + "`n")
$content += ']: '
Set-Content -LiteralPath $OutputPath -Value $content -Encoding ASCII
Set-Content -LiteralPath $IssuePath -Value $issues -Encoding ASCII
Write-Output ("active_headings={0}; reference_polynomials={1}; issues={2}" -f $headingNo,$rows.Count,$issues.Count)
